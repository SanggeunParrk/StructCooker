"""README.md for a server copy of the release, written from what the server really holds.

    python make_server_readme.py SERVER_NAME SERVER_LISTING RELEASE_LISTING [--metadata f1 f2 ...] > README.md

A listing is ``<db path without .lmdb> <entries> <total bytes>`` per line (from each DB's ``.meta.json``). The text never lists
a database the server lacks as present, and names what is not here (the release minus the server).
"""
from __future__ import annotations

import argparse

WHAT = {
    "pdb/cif_attached/cif_pdb_attached": "all PDB entries (unfiltered), 40% cluster ids",
    "pdb/cif_attached/cif_pdb_attached_train_20210930_res9_chain300": "PDB training structures: released by 2021-09-30, resolution <= 9 A, <= 300 chains",
    "pdb/valid/cif_attached_valid_1": "validation: released 2021-09-30 .. 2023-09-30",
    "pdb/valid/cif_attached_valid_2": "valid_1 minus clusters connected to train",
    "pdb/msa/a3m": "PDB protein MSAs, uncapped depth, keyed by seq_id",
    "pdb/msa/a3m_d16k": "PDB protein MSAs, depth <= 16,384",
    "pdb/msa/a3m_d2k": "PDB protein MSAs, depth <= 2,000",
    "pdb/msa/a3m_d512": "PDB protein MSAs, depth <= 512",
    "pdb/msa/msa_rna": "PDB RNA MSAs, uncapped (nhmmer; most are shallow)",
    "pdb/msa/msa_rna_d2k": "PDB RNA MSAs, depth <= 2,000",
    "pdb/msa/msa_rna_d512": "PDB RNA MSAs, depth <= 512",
    "chain/template/template": "PDB templates by chain (`<pdbid>_<chain>`; a chain may hold none)",
    "chain/template/seqid_template_mols": "per-sequence union of template mols (what `template` is copied out of)",
    "ccd/biomol_CCD": "ligand chemistry (CCD)",
    "ccd/ccd_inference": "CCD for inference",
    "interface/pdb/interface_pdb": "interface units of PDB assemblies (2..N chains), 40% interface clusters",
    "foldbench/msa_edge": "FoldBench from-scratch MSAs (`lmdb/foldbench/{msa,template,CCD}` is the FoldBench inference DB)",
    "distillation_short/cif_attached/cif_short_attached": "OpenFold3 short monomers: structures",
    "distillation_short/msa/msa_short_d2k": "OFD short MSAs, depth <= 2,000",
    "distillation_short/template/template_topn_short": "OFD short templates (key = record id)",
    "distillation_long/cif_attached/cif_long_attached": "OpenFold3 long monomers: structures (shards under `shards/`)",
    "distillation_long/msa/msa_long_d2k": "OFD long MSAs, depth <= 2,000",
    "distillation_long/template/template_topn_long": "OFD long templates (key = record id)",
    "distillation_rna/cif_attached/cif_rna_attached": "OpenFold3 RNA monomers: structures",
    "distillation_rna/msa/msa_rna": "OFD RNA MSAs, uncapped",
    "distillation_rna/msa/msa_rna_d2k": "OFD RNA MSAs, depth <= 2,000",
    "distillation_rna/msa/msa_rna_d512": "OFD RNA MSAs, depth <= 512",
    "distillation_rna/msa/msa_rna_d5k": "OFD RNA MSAs, depth <= 5,000",
    "distillation_disordered/cif_attached/cif_disordered_attached": "OpenFold3 disordered set: structures (with chain contacts)",
    "distillation_disordered/msa/msa_disordered": "disordered MSAs, uncapped",
    "distillation_disordered/msa/msa_disordered_d2k": "disordered MSAs, depth <= 2,000",
    "distillation_disordered/msa/msa_disordered_d512": "disordered MSAs, depth <= 512",
    "distillation_disordered/template/template_disordered": "disordered templates (key = `<pdbid>_<chain id of the structure>`)",
    "teddymer/cif_attached/cif_teddymer_attached": "Teddymer domain pairs: structures (510,454 cluster representatives; CC BY-NC 4.0)",
    "teddymer/msa/msa_teddymer": "Teddymer MSAs, uncapped",
    "teddymer/msa/msa_teddymer_d2k": "Teddymer MSAs, depth <= 2,000",
    "teddymer/template/template_teddymer": "Teddymer templates (key = seq_id)",
    "afdb_multimer/cif_attached/cif_afm_homodimer_attached": "AFDB multimer homodimers: structures",
    "afdb_multimer/cif_attached/cif_afm_heterodimer_attached": "AFDB multimer heterodimers: structures",
    "afdb_multimer/msa/msa_afm": "AFDB multimer MSAs, uncapped",
    "afdb_multimer/msa/msa_afm_d2k": "AFDB multimer MSAs, depth <= 2,000",
    "afdb_multimer/template/template_afm": "AFDB multimer templates (key = seq_id)",
    "design/db_monomer_long/design_monomer_long": "design training DB: OFD long monomers filtered, with SS fractions and foldseek clusters",
    "design/db_monomer_short/design_monomer_short": "design training DB: OFD short monomers",
    "design/db_teddymer/design_complex_teddymer": "design training DB: Teddymer complexes",
}
META = {
    "ab_seq_id.tsv": "antibody sequences (seq_id, 40% cluster, ANARCI chain type)",
    "cif_metadata.tsv": "resolution, dates, counts per assembly",
    "seq_cluster40_PDB.tsv": "40% sequence clusters of the PDB (sampling weights)",
    "seq_cluster40_OFD.tsv": "40% sequence clusters of the OpenFold3 sets",
    "seq_cluster40_TDM.tsv": "40% sequence clusters of Teddymer",
    "seq_cluster40_AFM.tsv": "40% sequence clusters of AFDB multimer",
    "interacting_seq_clusters.tsv": "interface clusters",
    "interacting_seq_ids.tsv": "interface sequence ids",
    "seq_id_map.tsv": "seq_id -> sequence, every database",
    "train_item_20210930.tsv": "MiniWorld train items (40% clusters, AF3 crystallization aids excluded)",
    "resources.lmdb": "MiniWorld ResourceLocator index (source -> DB paths)",
}


def read(path: str) -> dict[str, tuple[int, int]]:
    out = {}
    for line in open(path):
        p = line.split()
        if len(p) >= 2:
            out[p[0].removeprefix("lmdb/")] = (int(p[1]), int(p[2]) if len(p) > 2 and p[2].isdigit() else 0)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("name"), ap.add_argument("server"), ap.add_argument("release")
    ap.add_argument("--metadata", nargs="*", default=[])
    ap.add_argument("--extra", default="")
    a = ap.parse_args()
    srv, rel = read(a.server), read(a.release)
    print(f"# BioMolDB_20261014 on {a.name} — a subset of the BioMol release\n")
    print("A subset of the BioMol release (`/data/shared/cssb_data/BioMol_clean` on the BioAI cluster; release-check: 95 stages and the 7 "
          "template key contracts pass, 2026-11-02), copied with its paths unchanged, so a config written for the release works here with "
          "`BIOMOL_ROOT` pointing at this directory. `docs/` is the whole release's documentation (schemas, every database, provenance, "
          "licences); only the files below are here. Sequence clusters are 40% only.\n")
    print("| path under `lmdb/` | entries | what |\n|---|---:|---|")
    for k in sorted(srv):
        print(f"| `{k}.lmdb` | {srv[k][0]:,} | {WHAT.get(k, '')} |")
    if any(k.startswith("foldbench") for k in srv):
        print("| `foldbench/{msa,template,CCD}` | | FoldBench inference DB (`foldbench_build.json` says what was cut from the release) |")
    if a.metadata:
        print("\n| path | what |\n|---|---|")
        for m in a.metadata:
            print(f"| `metadata/{m}` | {META.get(m, '')} |")
    print("\nEach `<db>.lmdb` comes with its `.index.tsv` and `.meta.json`; sharded databases keep their shards under `<db>.lmdb/shards/` "
          "(see `shards.json`). Reading records: `docs/README.md`.")
    print("\nNot here (in the release): " + ", ".join(f"`{k}`" for k in sorted(set(rel) - set(srv))) + ".")
    if a.extra:
        print("\n" + a.extra)


if __name__ == "__main__":
    main()
