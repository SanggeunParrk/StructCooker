"""The final FASTA (kept pairs only, recovered light chains included) and the final pair table with seq_ids and clusters."""

from structcooker.instructions.readers.abag import abag_header as asd_header
from structcooker.instructions.readers.abag_final import (
    load_final_fasta,
    load_final_pairs,
)

H1, H2, L1, L2, AG = (
    "EVQLVESGGG",
    "QVQLQESGGG",
    "DIQMTQSPSS",
    "QSALTQPPAV",
    "MKTAYIAKQR",
)
COLS = "antigen_seq_id\theavy_seq_id\tlight_seq_id\tlight_sequence\tchain_status\tlight_source\tdatasets\tpolyprotein\tpdb_structure\tkeep"


def _files(tmp_path):
    (tmp_path / "map.tsv").write_text(
        f"P1\t{AG}\nP2\t{H1}\nP3\t{H2}\nP4\t{L1}\nP5\t{L2}\n"
    )
    rows = [
        "P1\tP2\tP4\t\tpaired_original\tasd\tpatents\t0\t\t1",
        f"P1\tP3\t\t{L2}\tfilled_hiv\tantibinder_hiv\thiv\t0\t\t1",  # a recovered light chain
        "P1\tP3\tP4\t\tpaired_original\tasd\tpatents\t0\t1ABC\t0",  # not kept: it has a structure
    ]
    (tmp_path / "pairs.tsv").write_text(COLS + "\n" + "\n".join(rows) + "\n")
    return tmp_path / "pairs.tsv", tmp_path / "map.tsv"


def test_fasta_holds_the_kept_chains_only(tmp_path):
    pairs, seq_map = _files(tmp_path)
    assert list(load_final_fasta(pairs, seq_map, "antigen").values()) == [AG]
    antibody = load_final_fasta(pairs, seq_map, "antibody")
    assert set(antibody.values()) == {H1, H2, L1, L2}  # the recovered light chain is in
    assert asd_header("H", H1) in antibody
    assert asd_header("L", L2) in antibody


def test_final_pairs_give_the_recovered_light_chain_its_id_and_clusters(tmp_path):
    pairs, seq_map = _files(tmp_path)
    (tmp_path / "ag.tsv").write_text("cASD_P1\tP1\n")
    (tmp_path / "ab.tsv").write_text("cASB_A2\tP2,P3\ncASB_A4\tP4\ncASB_A5\tP5\n")
    final = load_final_pairs(
        pairs, seq_map, tmp_path / "ag.tsv", {"cdr100": tmp_path / "ab.tsv"}
    )
    assert final[0] == (
        "antigen_seq_id",
        "heavy_seq_id",
        "light_seq_id",
        "chain_status",
        "light_source",
        "datasets",
        "antigen_cluster40",
        "heavy_cluster_cdr100",
        "light_cluster_cdr100",
    )
    assert len(final) == 3  # the pair with a structure is gone
    assert final[2][:3] == ("P1", "P3", "P5")
    assert final[2][-2:] == (
        "cASB_A2",
        "cASB_A5",
    )


def test_antigen_keep_list_narrows_the_fasta_and_the_pairs(tmp_path):
    pairs, seq_map = _files(tmp_path)
    keep = tmp_path / "keep.txt"
    keep.write_text("P9\n")  # an antigen no pair has
    assert load_final_fasta(pairs, seq_map, "antigen", keep) == {}
    keep.write_text("P1\n")
    assert list(load_final_fasta(pairs, seq_map, "antigen", keep).values()) == [AG]


def test_antigen_worklist_and_seed_map(tmp_path):
    from structcooker.instructions.readers.abag_final import (
        load_antigen_worklist,
        load_seq_id_maps,
    )

    pairs, seq_map = _files(tmp_path)
    assert load_antigen_worklist(pairs, seq_map) == {"P1": AG}
    (tmp_path / "a.tsv").write_text("P1\tAAA\nP2\tCCC\n")
    (tmp_path / "b.tsv").write_text("P2\tCCC\nP3\tDDD\n")
    assert load_seq_id_maps(
        [tmp_path / "a.tsv", tmp_path / "b.tsv", tmp_path / "missing.tsv"]
    ) == {"PAAA": "P1", "PCCC": "P2", "PDDD": "P3"}
    (tmp_path / "c.tsv").write_text("P7\tCCC\n")
    import pytest

    with pytest.raises(ValueError, match="two ids"):
        load_seq_id_maps([tmp_path / "a.tsv", tmp_path / "c.tsv"])
