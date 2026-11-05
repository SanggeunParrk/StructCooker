"""Writers that project recipe outputs to TSV/FASTA files."""

from pathlib import Path

from structcooker.instructions.writers import projections as w


def test_write_fasta_concatenates_all_entries(tmp_path):
    out = tmp_path / "sub" / "a.fa"
    w.write_fasta(
        {
            "1abc": {"fasta": {"k1": ">1abc_A\nMK\n", "k2": ">1abc_B\nGG\n"}},
            "2def": {"fasta": {"k": ">2def_A\nL\n"}},
        },
        out,
    )
    assert out.read_text() == ">1abc_A\nMK\n>1abc_B\nGG\n>2def_A\nL\n"


def test_write_seq_id_map_sorts_by_seq_id(tmp_path):
    out = tmp_path / "m.tsv"
    w.write_seq_id_map({"seq_id_map": {"P2": "GG", "P1": "MK"}}, out)
    assert out.read_text() == "P1\tMK\nP2\tGG\n"


def test_write_clusters_prefixes_db_code_and_joins_members(tmp_path):
    out = tmp_path / "c.tsv"
    data = {"db_code": "PDB", "cluster_dict": {"P1": ["P1", "P2"], "P3": ["P3"]}}
    w.write_clusters(data, out)
    assert out.read_text() == "cPDB_P1\tP1,P2\ncPDB_P3\tP3\n"
    out2 = tmp_path / "c2.tsv"
    w.write_seq_cluster_dict(data, out2)
    assert out2.read_text() == out.read_text()


def test_write_metadata_fills_missing_fields_with_na(tmp_path):
    out = tmp_path / "meta.tsv"
    w.write_metadata(
        {"1abc": {"metadata_dict": {"1_1_.": {"resolution": "2.0", "chain_num": "3"}}}},
        out,
    )
    header, row = out.read_text().splitlines()
    assert header.split("\t")[:3] == ["cif_id", "resolution", "deposition_date"]
    cols = row.split("\t")
    assert cols[0] == "1abc_1_1_."
    assert cols[1] == "2.0"
    assert cols[2] == "NA"
    assert cols[4] == "3"


def test_write_asd_pairs_variants_have_headers(tmp_path):
    a = tmp_path / "a.tsv"
    w.write_asd_pairs({"pairs": [("A1", "H1", "", "asd")]}, a)
    assert a.read_text() == (
        "antigen_seq_id\theavy_seq_id\tlight_seq_id\tdatasets\nA1\tH1\t\tasd\n"
    )
    b = tmp_path / "b.tsv"
    w.write_asd_pairs_chains(
        {"pairs": [("A1", "H1", "L1", "SEQ", "ok", "src", "asd")]}, b
    )
    lines = b.read_text().splitlines()
    assert lines[0].split("\t")[3:6] == [
        "light_sequence",
        "chain_status",
        "light_source",
    ]
    assert lines[1] == "A1\tH1\tL1\tSEQ\tok\tsrc\tasd"


def test_write_header_rows_writes_rows_verbatim(tmp_path):
    out = tmp_path / "r.tsv"
    w.write_header_rows({"rows": [("h1", "h2"), ("a", "b")]}, out)
    assert out.read_text() == "h1\th2\na\tb\n"


def test_write_filtered_seq_ids_unions_and_sorts(tmp_path):
    out = tmp_path / "f.tsv"
    w.write_filtered_seq_ids(
        {
            "x": {"filtered_seq_ids": {("P2", "P1"), ("P1", "P9")}},
            "y": {"filtered_seq_ids": {("P1", "P9"), ("P1", "P3")}},
        },
        out,
    )
    assert out.read_text() == "P1\tP3\nP1\tP9\nP2\tP1\n"


def test_write_filtered_seq_clusters_sorted(tmp_path):
    out = tmp_path / "f.tsv"
    w.write_filtered_seq_clusters(
        {"interacting_seq_clusters": [("b", "a"), ("a", "z"), ("b", "a")]}, out
    )
    assert out.read_text() == "a\tz\nb\ta\n"


def test_convert_to_type_merges_ligand_like_types_by_default():
    assert w.convert_to_type("P") == "Protein(L)"
    assert w.convert_to_type("B") == "Ligand"
    assert w.convert_to_type("L") == "Ligand"
    assert w.convert_to_type("?") == "Ligand"  # unknown code -> "Unknown" -> merged
    assert w.convert_to_type("L", merge_ligands=False) == "Non-polymer"
    assert w.convert_to_type("?", merge_ligands=False) == "Unknown"


def test_write_statistics_counts_monomer_types_and_sorted_interface_pairs(tmp_path):
    out = tmp_path / "s.txt"
    w.write_statistics(
        {
            "a": {
                "monomer_clusters": {"cPDB_P1", "cPDB_P2", "cPDB_L1"},
                "interface_clusters": {("cPDB_P1", "cPDB_L1")},
            },
            "b": {
                "monomer_clusters": {"cPDB_P2"},
                "interface_clusters": {("cPDB_L2", "cPDB_P1"), ("cPDB_P1", "cPDB_P3")},
            },
        },
        out,
    )
    text = out.read_text()
    assert "Protein(L),2\n" in text
    assert "Ligand,1\n" in text
    assert "Ligand-Protein(L),2\n" in text
    assert "Protein(L)-Protein(L),1\n" in text


def test_write_statistics_counts_each_cluster_once_across_entries(tmp_path):
    out = tmp_path / "s.txt"
    w.write_statistics(
        {
            "a": {
                "monomer_clusters": {"cPDB_P1", "cPDB_L1"},
                "interface_clusters": {("cPDB_P1", "cPDB_L1")},
            },
            "b": {
                "monomer_clusters": {"cPDB_P1"},
                "interface_clusters": {("cPDB_P1", "cPDB_L1")},
            },
        },
        out,
    )
    head, _, tail = out.read_text().partition("\nInterface Cluster Statistics:\n")
    assert sorted(head.splitlines()) == [
        "Ligand,1",
        "Monomer Cluster Statistics:",
        "Protein(L),1",
    ]
    assert tail == "Ligand-Protein(L),1\n"


def test_write_each_ccd_cif_buckets_by_first_letter(tmp_path):
    w.write_each_ccd_cif({"each_cif_lines": {"ATP": ["data_ATP\n", "x\n"]}}, tmp_path)
    assert (tmp_path / "A" / "ATP.cif").read_text() == "data_ATP\nx\n"


def _chain(pdb, chain):
    return {
        "pdb_id": pdb,
        "assembly_id": "1",
        "model_id": "1",
        "alt_id": ".",
        "chain_id": chain,
    }


def test_write_edge_node_drops_monomers_that_appear_in_an_interface(tmp_path):
    out = tmp_path / "e.tsv"
    w.write_edge_node(
        {
            "1abc": {
                "monomer_map": {
                    "c1": [_chain("1abc", "A")],
                    "c2": [_chain("1abc", "B")],
                    "c3": [_chain("1abc", "C")],
                },
                "interface_map": {
                    ("c1", "c2"): [(_chain("1abc", "A"), _chain("1abc", "B"))]
                },
            }
        },
        out,
    )
    lines = out.read_text().splitlines()
    assert lines[0].startswith("cluster1\tcluster2\tpdb_id")
    assert lines[1:] == [
        "c3\tNone\t1abc\t1\t1\t.\tC\tNone",
        "c1\tc2\t1abc\t1\t1\t.\tA\tB",
    ]


def test_unittest_writer_marks_none_as_passed(tmp_path):
    out = tmp_path / "u.tsv"
    w.unittest({"i1": {"results": None}, "i2": {"results": "boom"}}, out)
    assert out.read_text() == "i1\tpassed\ni2\tboom\n"


def test_monomer_seqid_and_afm_map_writers_sort_keys(tmp_path):
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    w.write_monomer_seqid_map({"monomer_seqid": {"E2": "P2", "E1": "P1"}}, a)
    w.write_afm_msa_seqid_map({"afm_msa_seqid": {"AF-2": "P2", "AF-1": "P1"}}, b)
    assert a.read_text() == "E1\tP1\nE2\tP2\n"
    assert b.read_text() == "AF-1\tP1\nAF-2\tP2\n"


def test_write_fasta_dict(tmp_path):
    out = tmp_path / "d.fa"
    w.write_fasta_dict({"fasta_dict": {"h1": "MK", "h2": "GG"}}, out)
    assert out.read_text() == ">h1\nMK\n>h2\nGG\n"
    assert isinstance(out, Path)


def test_cluster_mol_type_reads_the_seq_id_letter_and_the_legacy_form():
    assert w.cluster_mol_type("cPDB_P0000000000000000001") == "Protein(L)"
    assert w.cluster_mol_type("cABG_L000000000000000042") == "Ligand"
    assert (
        w.cluster_mol_type("cP00000000000000000007") == "Protein(L)"
    )  # no DB code: the type letter is at index 1
