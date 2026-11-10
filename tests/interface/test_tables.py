import hashlib

import pytest

from structcooker.instructions.interface import tables as t


def _h(seq):
    return hashlib.md5(f"P:{seq}".encode()).hexdigest()[:16]  # noqa: S324


def row(uid, pdb, seqs=(), ligands=(), date="2020-01-01", res="2.0", **kw):
    """Build a manifest row of protein nodes ``seqs`` and ligand nodes ``ligands``."""
    kinds = ["protein"] * len(seqs) + ["ligand"] * len(ligands)
    return {
        "unit_id": uid,
        "pdb_id": pdb,
        "resolution": res,
        "release_date": date,
        "type": "+".join(
            x
            for x in (
                f"{len(seqs)}protein" if seqs else "",
                f"{len(ligands)}ligand" if ligands else "",
            )
            if x
        ),
        "n_nodes": len(kinds),
        "n_protein": str(len(seqs)),
        "n_na": "0",
        "n_ligand": str(len(ligands)),
        "n_ion": "0",
        "has_ab": "0",
        "n_duplex_pairs": "0",
        "truncated": "False",
        "n_atoms": "100",
        "n_residues": "10",
        "n_interface_res": "4",
        "n_atom_pairs": "20",
        "n_water": "0",
        "n_bridging_water": "0",
        "node_kinds": "|".join(kinds),
        "node_hashes": "|".join([_h(s) for s in seqs] + [""] * len(ligands)),
        "node_ccd": "|".join([""] * len(seqs) + list(ligands)),
        **kw,
    }


def test_write_manifest_sorts_by_unit_and_blanks_none(tmp_path):
    data = {
        "b": {"manifest_row": {"unit_id": "2xyz_1_1_0", "x": None}},
        "a": {"manifest_row": {"unit_id": "1abc_1_1_0", "x": 3}},
    }
    out = tmp_path / "sub" / "manifest.tsv"
    t.write_manifest(data, out)
    assert out.read_text() == "unit_id\tx\n1abc_1_1_0\t3\n2xyz_1_1_0\t\n"


def test_write_manifest_of_nothing_writes_a_header(tmp_path):
    t.write_manifest({}, tmp_path / "m.tsv")
    assert (tmp_path / "m.tsv").read_text() == "unit_id\n"


def test_hashes_sequences_ids_and_clusters(tmp_path):
    rows = [row("u1", "1abc", ["MKV", "GGG"]), row("u2", "2xyz", ["MKV"], ["ATP"])]
    assert t.protein_hashes(rows) == {_h("MKV"), _h("GGG")}
    fasta = tmp_path / "p.fasta"
    fasta.write_text(">a\nMKV\n>b\nGGG\n>c\nWWW\n")
    h2s = t.hash_to_sequence(fasta, t.protein_hashes(rows))
    assert h2s == {_h("MKV"): "MKV", _h("GGG"): "GGG"}
    idmap = tmp_path / "map.tsv"
    idmap.write_text("P1\tMKV\nR2\tMKV\nP3\tWWW\nP4\tGGG\n")
    assert t.sequence_ids(idmap, {"MKV", "GGG"}) == {"MKV": "P1", "GGG": "P4"}
    clusters = tmp_path / "c.tsv"
    clusters.write_text("cPDB_P1\tP1,P4\ncPDB_P3\tP3\n")
    assert t.cluster_of(clusters, {"P1", "P4"}) == {"P1": "cPDB_P1", "P4": "cPDB_P1"}


def test_annotate_rows_adds_seq_ids_clusters_and_counts_unmapped():
    rows = [row("u1", "1abc", ["MKV", "GGG"], ["ATP"])]
    out, unmapped = t.annotate_rows(
        rows,
        {_h("MKV"): "P1", _h("GGG"): None},
        {"P1": "cPDB_P1"},
        {"1abc": "N"},
    )
    r = out[0]
    assert r["node_seq_ids"] == "P1||"
    assert r["node_c40"] == "cPDB_P1||"
    assert r["cluster_key40"] == "cPDB_P1"
    assert (r["n_unmapped_protein"], unmapped) == (1, 1)
    assert r["including_NA_entry"] == "N"


def test_cluster_key_is_the_sorted_set_of_clusters():
    rows = [row("u1", "1abc", ["A", "B", "C"])]
    out, _ = t.annotate_rows(
        rows,
        {_h("A"): "P1", _h("B"): "P2", _h("C"): "P3"},
        {"P1": "c2", "P2": "c1", "P3": "c2"},
        {},
    )
    assert out[0]["cluster_key40"] == "c1+c2"
    assert out[0]["including_NA_entry"] == ""


def test_write_manifest_full_end_to_end(tmp_path):
    rows = [
        row("1abc_1_1_0", "1abc", ["MKV", "GGG"]),
        row("2xyz_1_1_0", "2xyz", ["MKV"], ["ATP"]),
    ]
    t.write_manifest(
        {r["unit_id"]: {"manifest_row": r} for r in rows}, tmp_path / "m.tsv"
    )
    (tmp_path / "p.fasta").write_text(">a\nMKV\n>b\nGGG\n")
    (tmp_path / "map.tsv").write_text("P1\tMKV\nP2\tGGG\n")
    (tmp_path / "c.tsv").write_text("cPDB_P1\tP1,P2\n")
    (tmp_path / "meta.tsv").write_text(
        "cif_id\tincluding_NA\n1abc_1_1_.\t0\n1abc_2_1_.\t1\n2xyz_1_1_.\t0\n3zzz_1_1_.\t0\n"
    )
    status = t.write_manifest_full(
        tmp_path / "m.tsv",
        tmp_path / "p.fasta",
        tmp_path / "map.tsv",
        tmp_path / "c.tsv",
        tmp_path / "meta.tsv",
        tmp_path / "full.tsv",
        tmp_path / "out" / "report.md",
    )
    assert status.startswith("2 units")
    full = t.read_tsv(tmp_path / "full.tsv")
    assert [r["unit_id"] for r in full] == ["1abc_1_1_0", "2xyz_1_1_0"]
    assert full[0]["node_seq_ids"] == "P1|P2"
    assert full[0]["cluster_key40"] == "cPDB_P1"
    assert full[0]["including_NA_entry"] == "0"
    assert list(full[0])[-5:] == [
        "node_seq_ids",
        "node_c40",
        "cluster_key40",
        "n_unmapped_protein",
        "including_NA_entry",
    ]
    report = (tmp_path / "out" / "report.md").read_text()
    assert "Units: **2** from **2** PDB entries (1 of the 3 entries" in report
    assert "| 2protein | 1 | 50.0% |" in report
    assert "ligands: ATP 1" in report


@pytest.mark.parametrize(
    ("date", "res", "expected"),
    [
        ("2020-05-01", "2.0", ("train", "")),
        ("2020-05-01", "9.0", ("train", "")),
        ("2020-05-01", "9.5", ("none", "resolution>9")),
        ("2021-09-30", "4.5", ("valid", "")),
        ("2023-01-12", "4.5", ("valid", "")),
        ("2023-01-13", "4.5", ("test", "")),
        ("2022-01-01", "4.6", ("none", "resolution>4.5")),
        ("2022-01-01", "None", ("none", "no_resolution")),
        ("2022-01-01", "", ("none", "no_resolution")),
        ("None", "2.0", ("none", "no_release_date")),
    ],
)
def test_classify_unit(date, res, expected):
    assert t.classify_unit(date, res) == expected


def _split_fixture():
    def r(uid, pdb, c40, date, ligands=(), res="2.0"):
        key = "+".join(sorted(set(c40)))
        return row(
            uid,
            pdb,
            ["S"] * len(c40),
            ligands,
            date,
            res,
            node_c40="|".join([*c40, *([""] * len(ligands))]),
            cluster_key40=key,
        )

    return [
        r("t1", "1aaa", ["c1", "c2"], "2020-01-01", ["ATP"]),
        r("t2", "1bbb", ["c3"], "2020-01-01"),
        r("v1", "2aaa", ["c1", "c2"], "2022-01-01", ["ATP"]),  # nothing new
        r("v2", "2bbb", ["c1", "c9"], "2022-01-01", ["NAD"]),  # new pair, new ligand
        r("v3", "2ccc", ["c8"], "2022-01-01"),  # novel protein, no ligand
        r("e1", "3aaa", ["c1"], "2024-01-01", res="9.0"),  # excluded
        r("s1", "3bbb", ["c3"], "2024-01-01"),
    ]


def test_split_rows_flag_leakage_against_the_train_units():
    out = {r["unit_id"]: r for r in t.split_rows(_split_fixture())}
    assert [out[u]["split"] for u in ("t1", "v1", "e1", "s1")] == [
        "train",
        "valid",
        "none",
        "test",
    ]
    v1, v2, v3 = out["v1"], out["v2"], out["v3"]
    assert (v1["n_prot_c40_in_train"], v1["novel_protein"], v1["novel_pair"]) == (
        2,
        0,
        0,
    )
    assert v1["ligands_seen_in_train"] == 1
    assert (v2["n_prot_c40_in_train"], v2["novel_protein"], v2["novel_pair"]) == (
        1,
        0,
        1,
    )
    assert v2["ligands_seen_in_train"] == 0
    assert (v3["n_prot_c40_in_train"], v3["novel_protein"], v3["novel_pair"]) == (
        0,
        1,
        1,
    )
    assert v3["ligands_seen_in_train"] == ""
    assert out["t1"]["novel_protein"] == out["e1"]["novel_pair"] == ""


def test_split_rows_does_not_modify_its_input():
    rows = _split_fixture()
    t.split_rows(rows)
    assert "split" not in rows[0]


def test_write_split_writes_the_table_and_the_summary(tmp_path):
    rows = [
        {**r, "node_seq_ids": "", "n_unmapped_protein": 0} for r in _split_fixture()
    ]
    t._write_rows(tmp_path / "full.tsv", list(rows[0]), rows)
    status = t.write_split(
        tmp_path / "full.tsv", tmp_path / "split.tsv", tmp_path / "s" / "summary.md"
    )
    assert status == "train 2, valid 3, test 1, none 1"
    out = t.read_tsv(tmp_path / "split.tsv")
    assert list(out[0]) == t.SPLIT_COLUMNS
    assert {r["unit_id"]: r["split"] for r in out}["v2"] == "valid"
    summary = (tmp_path / "s" / "summary.md").read_text()
    assert "| train | 2 | 2 |" in summary
    assert "excluded: resolution>4.5 1" in summary
    assert "| valid | 1 | 2 |" in summary
