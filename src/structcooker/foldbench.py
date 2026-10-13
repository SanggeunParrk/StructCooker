"""FoldBench inference DBs: the slice of the PDB-wide BioMol LMDBs that a FoldBench target set needs.

FoldBench (BEAM-Labs, 9 tasks, 1,522 assemblies) is predicted from three kinds of lookups, all by key:

=============  ==========================  ==========================================  ==============================
asset          key                         value                                       source (full-release LMDB)
=============  ==========================  ==========================================  ==============================
protein MSA    ``P`` + 20 digits           ``{msa_dict}`` (schema E)                   ``lmdb/pdb/msa/a3m.lmdb``
RNA MSA        ``R`` + 20 digits           ``{msa_dict}`` (schema E, a3m_type=rna)     ``lmdb/pdb/msa/msa_rna.lmdb``
template       ``{pdbid}_{chain}``         ``{template_mols}`` (schema D)              ``lmdb/chain/template/template.lmdb``
CCD            component code              inference CCD record (``atoms``/``residues``/``chains``)  ``--ccd`` (see below)
=============  ==========================  ==========================================  ==============================

The sequences were already searched and the templates already resolved for the whole PDB, so a FoldBench DB is a
**key-list subset**, not a new search. This module derives the key lists from a FoldBench package
(``meta/seqids_*.txt``, ``meta/chain_seqid_map.tsv``, ``inputs/<target>/data.yaml``), copies exactly those records
byte-for-byte (the values are never decoded), and verifies that every lookup a ``data.yaml`` makes resolves.

The few sequences the PDB-wide DBs do not cover (FoldBench's 3 chromophore-abbreviated entries and a 99-mer) were
searched from scratch once with :func:`structcooker.instructions.transforms.msa.run_msa_search`; ``--fallback`` points
at a DB that already holds them, so a rebuild does not repeat that search.

Sources are the release tree (``BIOMOL_ROOT``, e.g. ``/data/shared/cssb_data/BioMol_clean``), at the
paths its db configs write (``RELEASE_LAYOUT``; a test keeps the two in step). The CCD is not the release's
``ccd/biomol_CCD.lmdb`` (schema G, ``{chem_comp_dict}``): inference reads CCDMol graphs
(``atoms``/``residues``/``chains``), reshaped from it into ``foldbench/ccd_inference.lmdb``
(:mod:`structcooker.instructions.transforms.ccdmol`, scripts/staging/build_inference_ccd.py). The four
from-scratch edge MSAs are kept in ``foldbench/msa_edge.lmdb``; ``structcooker foldbench build`` uses
both by default.

Run as ``structcooker foldbench keys|build|verify`` (see :mod:`structcooker.cli`) or
``python -m structcooker.foldbench build PACKAGE OUT``.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import lmdb

KEY_WIDTH = 21  # a seq_id key: one letter (P protein, R RNA, D/Q/N other nucleic acids) + 20 zero-padded digits
# ``(ACE)`` / ``(NAG)`` inline residues and ligands. The glycosidic positions ``|(1,2)`` contain a comma, so they never match.
_INLINE_CODE = re.compile(r"\(([A-Za-z0-9]{1,5})\)")
# value-less key prefixes the lab's PDB-wide MSA DBs do not carry (DNA has no MSA, FoldBench matches AF3 there)
_NO_MSA_PREFIXES = ("D", "N", "Q")


def lmdb_key(seq_id: str) -> str:
    """``R0021166`` -> ``R00000000000000021166``: the 21-character form the MSA LMDBs are keyed by."""
    s = seq_id.strip()
    if len(s) == KEY_WIDTH:
        return s
    return s[0] + s[1:].zfill(KEY_WIDTH - 1)


@dataclass
class FoldBenchKeys:
    """The keys one FoldBench package looks up, by asset."""

    protein: list[str] = field(default_factory=list)  # includes the from-scratch (edge) sequences: ``--fallback`` supplies them
    rna: list[str] = field(default_factory=list)
    other: list[str] = field(default_factory=list)  # D / N / Q keys: DNA and short oligos, no MSA by design (not copied)
    template: list[str] = field(default_factory=list)
    ccd: list[str] = field(default_factory=list)

    def counts(self) -> dict[str, int]:
        return {k: len(v) for k, v in asdict(self).items()}


def _read_ids(path: Path) -> list[str]:
    if not path.is_file():
        return []
    return sorted({lmdb_key(ln) for ln in path.read_text().splitlines() if ln.strip()})


def load_keys(package: Path) -> FoldBenchKeys:
    """Derive the key lists from ``package`` = the FoldBench directory (``meta/``, optionally ``inputs/``)."""
    meta = package / "meta"
    keys = FoldBenchKeys(
        protein=_read_ids(meta / "seqids_protein.txt"),
        rna=_read_ids(meta / "seqids_rna.txt"),
    )
    codes: set[str] = set()
    chain_map = meta / "chain_seqid_map.tsv"
    if chain_map.is_file():
        with chain_map.open() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                codes.update(_INLINE_CODE.findall(row.get("sequence", "")))
    keys.ccd = sorted(codes)

    template: set[str] = set()
    msa: dict[str, set[str]] = {"P": set(), "R": set(), "other": set()}
    index = package / "inputs" / "index.txt"
    if index.is_file():
        import yaml

        for target in (ln.strip() for ln in index.read_text().splitlines() if ln.strip()):
            spec = package / "inputs" / target / "data.yaml"
            if not spec.is_file():
                continue
            doc = yaml.safe_load(spec.read_text())
            template.update(str(v) for v in (doc.get("template") or {}).values() if v)
            for v in (doc.get("a3m_key") or {}).values():
                msa[str(v)[0] if str(v)[0] in "PR" else "other"].add(lmdb_key(str(v)))
        # the lookups a ``data.yaml`` actually makes are authoritative; they also name the from-scratch (edge) sequences
        keys.protein, keys.rna, keys.other = sorted(msa["P"]), sorted(msa["R"]), sorted(msa["other"])
    keys.template = sorted(template)
    return keys


@dataclass
class CopyReport:
    """What one subset copy did."""

    name: str
    target: str
    sources: list[str]
    requested: int
    written: int = 0
    from_fallback: int = 0
    missing: list[str] = field(default_factory=list)
    value_bytes: int = 0
    seconds: float = 0.0


def _open_read(path: Path):
    from datacooker.lmdb.sharded import open_env  # an ordinary LMDB or an immutable shard collection

    return open_env(path, readonly=True, lock=False, readahead=False)


def copy_keys(dst: Path, keys: Iterable[str], sources: Sequence[Path], *, name: str = "") -> CopyReport:
    """Write ``keys`` out of the first of ``sources`` that holds each one into a new LMDB at ``dst`` (values byte-for-byte).

    ``sources[0]`` is the primary DB; later ones only supply keys the primary does not have. Keys no source holds are
    reported, not an error. ``dst`` must not exist yet (an existing DB is never replaced).
    """
    t0 = time.time()
    wanted = sorted(set(keys))
    report = CopyReport(name=name or dst.name, target=str(dst), sources=[str(s) for s in sources], requested=len(wanted))
    if (dst / "data.mdb").exists():
        msg = f"Refusing to replace existing database: {dst}"
        raise FileExistsError(msg)
    envs = [_open_read(s) for s in sources]
    txns = [e.begin() for e in envs]
    try:
        where: dict[str, int] = {}
        total = 0
        for k in wanted:
            for i, txn in enumerate(txns):
                v = txn.get(k.encode())
                if v is not None:
                    where[k] = i
                    total += len(v)
                    break
            else:
                report.missing.append(k)
        dst.mkdir(parents=True, exist_ok=True)
        with lmdb.open(str(dst), map_size=max(1 << 26, int(total * 1.5) + (1 << 24)), subdir=True) as out:
            with out.begin(write=True) as wtxn:
                for k, i in where.items():
                    wtxn.put(k.encode(), bytes(txns[i].get(k.encode())))
                    report.written += 1
                    report.from_fallback += i > 0
        report.value_bytes = total
    finally:
        for txn in txns:
            txn.abort()
        for env in envs:
            env.close()
    report.seconds = round(time.time() - t0, 2)
    return report


# Release-relative path of each source, as the db config that builds it writes it.
RELEASE_LAYOUT = {
    "protein_msa": ("msa/a3m", "lmdb/pdb/msa/a3m.lmdb"),
    "rna_msa": ("msa/msa_rna", "lmdb/pdb/msa/msa_rna.lmdb"),
    "template": ("template/chain_template", "lmdb/chain/template/template.lmdb"),
}


@dataclass
class Sources:
    """The full-release DBs to cut the FoldBench DB from."""

    protein_msa: Path
    rna_msa: Path
    template: Path
    ccd: Path | None = None
    protein_fallback: list[Path] = field(default_factory=list)  # DBs holding the from-scratch (edge) sequences

    @classmethod
    def from_release(cls, release: Path, *, ccd: Path | None = None, fallback: Sequence[Path] = ()) -> "Sources":
        """``release`` = a release tree (``BIOMOL_ROOT``, e.g. ``/data/shared/cssb_data/BioMol_clean``)."""
        paths = {name: release / rel for name, (_db, rel) in RELEASE_LAYOUT.items()}
        return cls(paths["protein_msa"], paths["rna_msa"], paths["template"], ccd, list(fallback))


def build(package: Path, out: Path, sources: Sources) -> dict:
    """Cut the four FoldBench DBs under ``out`` and write ``out/foldbench_build.json``."""
    keys = load_keys(package)
    # the from-scratch (edge) sequences are asked of the primary first, then of the fallbacks
    reports = [
        copy_keys(out / "msa" / "a3m.lmdb", keys.protein, [sources.protein_msa, *sources.protein_fallback], name="protein_msa"),
        copy_keys(out / "msa" / "a3m_rna.lmdb", keys.rna, [sources.rna_msa], name="rna_msa"),
        copy_keys(out / "template" / "template.lmdb", keys.template, [sources.template], name="template"),
    ]
    if sources.ccd is not None:
        reports.append(copy_keys(out / "CCD" / sources.ccd.name, keys.ccd, [sources.ccd], name="ccd"))
    summary = {"package": str(package), "keys": keys.counts(), "copies": [asdict(r) for r in reports]}
    (out / "foldbench_build.json").write_text(json.dumps(summary, indent=1))
    return summary


def verify(package: Path, out: Path, *, decode_sample: int = 20) -> dict:
    """Check every lookup the package's ``data.yaml`` files make resolves in the DBs under ``out``.

    Reports slots by key prefix (found / not found) and decodes ``decode_sample`` records of each DB. Not-found ``D``/``N``/``Q``
    keys are expected (DNA and short nucleic-acid oligos have no MSA) and are counted separately from real gaps.
    """
    import yaml

    msa_dbs = [out / "msa" / "a3m.lmdb", out / "msa" / "a3m_rna.lmdb"]
    tmpl_db = out / "template" / "template.lmdb"
    envs = {p: lmdb.open(str(p), readonly=True, lock=False, readahead=False) for p in [*msa_dbs, tmpl_db] if p.exists()}
    txns = {p: e.begin() for p, e in envs.items()}

    def has(db: Path, key: str) -> bool:
        return db in txns and txns[db].get(key.encode()) is not None

    msa = {"found": 0, "expected_missing": 0, "missing": 0}
    msa_missing: list[str] = []
    tmpl = {"found": 0, "missing": 0}
    ids = [ln.strip() for ln in (package / "inputs" / "index.txt").read_text().splitlines() if ln.strip()]
    for target in ids:
        spec = yaml.safe_load((package / "inputs" / target / "data.yaml").read_text())
        for k in (spec.get("a3m_key") or {}).values():
            if any(has(db, k) for db in msa_dbs):
                msa["found"] += 1
            elif str(k)[0] in _NO_MSA_PREFIXES:
                msa["expected_missing"] += 1
            else:
                msa["missing"] += 1
                msa_missing.append(f"{target}:{k}")
        for k in (spec.get("template") or {}).values():
            if has(tmpl_db, k):
                tmpl["found"] += 1
            else:
                tmpl["missing"] += 1  # a chain with no template before the cutoff has no key: not an error
    decoded = {}
    try:
        from structcooker.instructions.transforms.codecs import from_bytes

        if decode_sample <= 0:
            raise ImportError

        for p, txn in txns.items():
            n = bad = 0
            for _, v in txn.cursor():
                try:
                    from_bytes(bytes(v))
                except Exception:  # noqa: BLE001 - a diagnostic count
                    bad += 1
                n += 1
                if n >= decode_sample:
                    break
            decoded[str(p)] = {"decoded": n - bad, "failed": bad}
    except ImportError:
        decoded = {"skipped": "decode_sample=0 or codecs not importable"}
    for e in envs.values():
        e.close()
    return {"targets": len(ids), "msa_slots": msa, "msa_missing_examples": msa_missing[:10], "template_slots": tmpl, "decode": decoded}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="structcooker.foldbench", description=__doc__.split("\n", 1)[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    k = sub.add_parser("keys", help="print and optionally save the key lists of a FoldBench package")
    k.add_argument("package", type=Path)
    k.add_argument("--out", type=Path)
    b = sub.add_parser("build", help="cut the FoldBench DBs out of the full-release LMDBs")
    b.add_argument("package", type=Path)
    b.add_argument("out", type=Path)
    b.add_argument("--release", type=Path, required=True, help="release tree (BIOMOL_ROOT), e.g. .../BioMol_clean")
    b.add_argument("--ccd", type=Path, help="inference CCD LMDB (atoms/residues/chains, single-atom components in)")
    b.add_argument("--fallback", type=Path, action="append", default=[], help="DB holding the from-scratch edge MSAs")
    v = sub.add_parser("verify", help="check every data.yaml lookup resolves in OUT")
    v.add_argument("package", type=Path)
    v.add_argument("out", type=Path)
    a = p.parse_args(argv)
    if a.cmd == "keys":
        keys = load_keys(a.package)
        print(json.dumps(keys.counts()))
        if a.out:
            a.out.write_text(json.dumps(asdict(keys), indent=1))
    elif a.cmd == "build":
        out = build(a.package, a.out, Sources.from_release(a.release, ccd=a.ccd, fallback=a.fallback))
        for c in out["copies"]:
            print(f"{c['name']:<12} requested {c['requested']:>5}  written {c['written']:>5}  fallback {c['from_fallback']:>3}  "
                  f"missing {len(c['missing']):>4}  {c['value_bytes'] / 2**20:8.1f} MiB  {c['seconds']}s")
    else:
        print(json.dumps(verify(a.package, a.out), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
