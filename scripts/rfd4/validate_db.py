#!/usr/bin/env python3
"""Check a built RFD4-Proteina set DB (``lmdb/rfd4/<set>/<set>.lmdb``) against its worklist and sidecar.

For every key (or a random ``--sample``): the value decodes, holds a non-empty assembly and a deposition date, and its
chain / residue / atom counts agree with what the Hub's ``examples.parquet`` says where it says it (a column that is
present in the sidecar). Also reports the keys of the worklist the DB lacks (records that failed to convert) and the
distribution of entity types. Heavy for the large sets: run through sbatch (``validate_db.sbatch``).
"""

from __future__ import annotations

import argparse
import os
import random
from collections import Counter

import lmdb
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from datacooker.lmdb.sharded import open_env

from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.instructions.transforms.rfd4_worklist import COMPONENTS

# the Hub counts a modified / non-standard residue of a peptide apart (the NCAA sets: "other residues"); a BioMol record's
# protein residues are all of them
PROTEIN_EXTRA = ("n_non_standard_amino_acid_residues", "n_other_residues")
# sidecar column -> (what it counts in the BioMol value)
RESIDUE_COLUMNS = {
    "n_standard_amino_acid_residues": "protein",
    "n_rna_residues": "rna",
    "n_dna_residues": "dna",
}


def count(value: dict) -> dict[str, int]:
    """Chains, residues and atoms of the record's first assembly, residues split by the chain's entity type."""
    a = next(iter(value["assembly_dict"].values()))
    chains = a["chains"]["nodes"]
    etype = chains["entity_type"]["value"]
    ptr = a["index_table"]["chain_res_indptr"]
    out = Counter(chains=len(etype), atoms=len(a["atoms"]["nodes"]["id"]["value"]))
    for i, t in enumerate(etype):
        n = int(ptr[i + 1] - ptr[i])
        out["residues"] += n
        out[f"type:{t}"] += 1
        if t == "polypeptide(L)":
            out["protein"] += n
        elif t == "polyribonucleotide":
            out["rna"] += n
        elif t == "polydeoxyribonucleotide":
            out["dna"] += n
    return dict(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("name", choices=list(COMPONENTS))
    ap.add_argument(
        "--sample",
        type=int,
        default=0,
        help="check this many random keys (default: all)",
    )
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    root = os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean")
    db = open_env(  # an ordinary LMDB or a shard collection (OMol)
        f"{root}/lmdb/rfd4/{a.name}/{a.name}.lmdb",
        readonly=True,
        lock=False,
        readahead=False,
    )
    wl = lmdb.open(
        f"{root}/materials/intermediate/rfd4/{a.name}/worklist.lmdb",
        readonly=True,
        lock=False,
    )
    with wl.begin() as t:
        wkeys = [bytes(k).decode() for k in t.cursor().iternext(values=False)]
    with db.begin() as t:
        keys = [bytes(k).decode() for k in t.cursor().iternext(values=False)]
    lacking = set(wkeys) - set(keys)
    extra = set(keys) - set(wkeys)
    print(
        f"{a.name}: worklist {len(wkeys)}, db {len(keys)}, lacking {len(lacking)}, not in worklist {len(extra)}"
    )
    if lacking:
        print("  e.g. lacking:", sorted(lacking)[:5])
    if a.sample and a.sample < len(keys):
        keys = random.Random(a.seed).sample(keys, a.sample)
    side_path = f"{root}/metadata/rfd4_{a.name}_examples.parquet"
    names = pq.ParquetFile(side_path).schema_arrow.names
    have = [c for c in [*list(RESIDUE_COLUMNS), "n_atoms"] if c in names]
    extra_cols = [c for c in PROTEIN_EXTRA if c in names]
    # only the rows of the keys checked (OMol's sidecar has 53M rows)
    side = pq.read_table(side_path, columns=["key", *have, *extra_cols])
    side = side.filter(pc.is_in(side.column("key"), value_set=pa.array(keys)))
    first: dict[str, dict] = {}
    cols = {c: side.column(c).to_pylist() for c in have + extra_cols}
    for i, k in enumerate(side.column("key").to_pylist()):
        first.setdefault(k, {c: cols[c][i] for c in cols})
    bad = Counter()
    mism = Counter()
    types = Counter()
    ex: dict[str, str] = {}
    n_ok = 0
    with db.begin() as t:
        for k in keys:
            try:
                v = from_bytes(bytes(t.get(k.encode())))
                assert v["assembly_dict"], "empty assembly_dict"
                assert "deposition_date" in v["metadata_dict"], "no deposition_date"
                c = count(v)
            except Exception as e:  # noqa: BLE001
                bad[f"{type(e).__name__}: {str(e)[:80]}"] += 1
                continue
            n_ok += 1
            for kk, vv in c.items():
                if kk.startswith("type:"):
                    types[kk] += vv
            row = first.get(k, {})
            for col in have:
                want = row.get(col)
                if col == "n_standard_amino_acid_residues" and want is not None:
                    # the record's protein residues include the non-standard and unknown ones the Hub counts apart
                    want += sum(row.get(c) or 0 for c in extra_cols)
                got = c.get(RESIDUE_COLUMNS.get(col, "atoms"))
                if want is None or got is None:
                    continue  # the Hub left the count empty (null means not audited), or the record has no such chains
                # residue counts must agree; the record may hold more atoms than the Hub counted (the BioMol residue
                # lists the CCD's atoms, unresolved ones included), never fewer
                bad_count = got < int(want) if col == "n_atoms" else int(want) != got
                if bad_count:
                    mism[col] += 1
                    ex.setdefault(col, f"{k}: sidecar {want} vs record {got}")
    print(f"  decoded and sane: {n_ok}/{len(keys)}; failures: {dict(bad)}")
    print(f"  entity types (chains): {dict(types)}")
    for col in have:
        print(
            f"  {col}: {mism[col]} of {n_ok} {'short' if col == 'n_atoms' else 'differ'}{('  e.g. ' + ex[col]) if col in ex else ''}"
        )


if __name__ == "__main__":
    main()
