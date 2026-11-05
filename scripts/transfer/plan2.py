"""Plan the 2026-10-28 refresh: re-attached (40%) cif_attached, metadata, CCD, design DBs (+ extras for H100).

    python plan2.py h100|cssb K     # writes units_<target>_rank<r>.txt (ROOT<TAB>REL<TAB>PREFIX per file),
                                    # files_<target>.tsv (dest path, bytes), prefixes_<target>.txt
Group 1 goes before group 2 in every rank (group 2 is only sent if time is left). PREFIX is the destination
directory under the remote root ('.' = the release layout itself).
"""

import os
import sys
from pathlib import Path

REL = Path("/data/shared/cssb_data/BioMol_clean")
DES = Path("/data/psk6950/design")
TR = Path("/data/psk6950/transfer")
HERE = Path(__file__).resolve().parent

PDB_CIF = [
    "lmdb/pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300.lmdb",
    "lmdb/pdb/valid/cif_attached_valid_1.lmdb",
    "lmdb/pdb/valid/cif_attached_valid_2.lmdb",
]
CCD = ["lmdb/ccd/ccd_inference.lmdb", "lmdb/ccd/biomol_CCD.lmdb"]
OTHER_CIF = [
    "lmdb/distillation_short/cif_attached/cif_short_attached.lmdb",
    "lmdb/distillation_long/cif_attached/cif_long_attached.lmdb",
    "lmdb/distillation_rna/cif_attached/cif_rna_attached.lmdb",
    "lmdb/distillation_disordered/cif_attached/cif_disordered_attached.lmdb",
]
META_H100 = [
    "metadata/seq_cluster40_PDB.tsv",
    "metadata/seq_cluster40_OFD.tsv",
    "metadata/cif_metadata.tsv",
    "metadata/train_item_20210930.tsv",
]
META_CSSB = [
    "metadata/seq_cluster40_PDB.tsv",
    "metadata/cif_metadata.tsv",
    "metadata/train_item_20210930.tsv",
]
GROUP2_H100 = [
    "lmdb/interface/pdb",
    "lmdb/teddymer/cif_attached",
    "lmdb/teddymer/template",
    "lmdb/teddymer/msa",
    "lmdb/afdb_multimer/cif_attached",
    "lmdb/afdb_multimer/template",
    "lmdb/afdb_multimer/msa",
]


def walk(p: Path):
    if p.is_file():
        yield p
        return
    for dp, _, fns in os.walk(p):
        for fn in fns:
            if fn != "lock.mdb":
                yield Path(dp) / fn


def add_release(items, rels, group):
    for rel in rels:
        p = REL / rel
        assert p.exists(), rel
        for f in walk(p):
            items.append((group, REL, str(f.relative_to(REL)), ".", f.stat().st_size))
        for suffix in (".index.tsv", ".meta.json"):
            s = Path(str(p) + suffix)
            if s.is_file():
                items.append(
                    (group, REL, str(s.relative_to(REL)), ".", s.stat().st_size)
                )


def main(target: str, k: int) -> None:
    items: list = []
    if (
        target == "h100x"
    ):  # shards.json of cif_long_attached points to ../.cif_long_attached.build/<hash>/ (outside the lmdb dir)
        b = REL / "lmdb/distillation_long/cif_attached/.cif_long_attached.build"
        for f in walk(b):
            items.append((1, REL, str(f.relative_to(REL)), ".", f.stat().st_size))
    elif target in (
        "h100c",
        "cssbc",
    ):  # only the PDB complex design DB (built after the main refresh)
        pass
    elif target == "h100":
        add_release(items, PDB_CIF + CCD + OTHER_CIF + META_H100, 1)
    else:
        add_release(items, [*PDB_CIF, "lmdb/ccd/ccd_inference.lmdb", *META_CSSB], 1)
    for f in (
        walk(TR / f"{target}_meta" / "metadata" / "resources.lmdb")
        if target in ("h100", "cssb")
        else []
    ):
        items.append(
            (
                1,
                TR / f"{target}_meta",
                str(f.relative_to(TR / f"{target}_meta")),
                ".",
                f.stat().st_size,
            )
        )
    for f in walk(REL / "docs") if target in ("h100", "cssb") else []:
        items.append((1, REL, str(f.relative_to(REL)), ".", f.stat().st_size))
    designs = (
        ["db_monomer_long", "db_monomer_short", "db_teddymer"]
        if target == "h100"
        else []
    )
    if target in ("h100c", "cssbc"):
        designs = []
    if (DES / "db_complex_pdb" / "manifest.tsv").exists():
        designs.append("db_complex_pdb")
    for d in designs:
        for f in walk(DES / d):
            items.append(
                (
                    1,
                    DES / d,
                    str(f.relative_to(DES / d)),
                    f"lmdb/design/{d}",
                    f.stat().st_size,
                )
            )
    if target == "h100":
        for rel in GROUP2_H100:
            for f in walk(REL / rel):
                items.append(
                    (
                        3 if "afdb_multimer" in rel else 2,
                        REL,
                        str(f.relative_to(REL)),
                        ".",
                        f.stat().st_size,
                    )
                )
        for rel in GROUP2_H100:
            for suffix in (".index.tsv", ".meta.json"):
                pass
    items = list(
        {(r, rel, pre): (g, r, rel, pre, b) for g, r, rel, pre, b in items}.values()
    )
    ranks = [[] for _ in range(k)]
    loads = [0] * k
    for g in (1, 2, 3):
        for it in sorted((x for x in items if x[0] == g), key=lambda x: -x[4]):
            r = loads.index(min(loads))
            ranks[r].append(it)
            loads[r] += it[4]
    for r in range(k):
        (HERE / f"units_{target}_rank{r}.txt").write_text(
            "".join(f"{root}\t{rel}\t{pre}\n" for _, root, rel, pre, _ in ranks[r])
        )
        g1 = sum(i[4] for i in ranks[r] if i[0] == 1)
        print(
            f"rank {r}: {len(ranks[r])} files, group1 {g1 / 1e9:.1f} GB, total {loads[r] / 1e9:.1f} GB"
        )
    with (HERE / f"files_{target}.tsv").open("w") as out:
        out.writelines(
            f"{(pre + '/' if pre != '.' else '')}{rel}\t{b}\n"
            for _, _, rel, pre, b in sorted(items, key=lambda x: x[2])
        )
    (HERE / f"prefixes_{target}.txt").write_text(
        "".join(sorted({f"{p}\n" for *_, p, _ in items}))
    )
    for g in (1, 2, 3):
        print(
            f"group {g}: {sum(1 for x in items if x[0] == g)} files, {sum(x[4] for x in items if x[0] == g) / 1e12:.3f} TB"
        )


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
