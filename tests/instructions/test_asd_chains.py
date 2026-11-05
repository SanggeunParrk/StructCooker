"""ASD chain FASTAs, the antibody-only split and batched CDR extraction."""
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from structcooker.instructions.readers.asd import load_asd_fasta
from structcooker.instructions.transforms.sequence import (
    _extract_H3L3_batch,
    _extract_H3L3_sequence,
    separate_antibody_sequences,
)

# trastuzumab (Chothia CDR-H3 / L3 are well defined)
VH = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
ANTIGEN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"

_MAP = pa.map_(pa.string(), pa.string())


def _table(rows):
    chains = [
        {"antibody_chains": ab, "target_chains": tg, "antibody_chains_numbered": {}} for ab, tg in rows
    ]
    struct = pa.struct(
        [("antibody_chains", _MAP), ("target_chains", _MAP), ("antibody_chains_numbered", pa.map_(pa.string(), _MAP))],
    )
    return pa.table({"chains": pa.array(chains, type=struct)})


@pytest.fixture
def asd_dir(tmp_path):
    rows = [
        ([("heavy", VH), ("light", VL)], [("full", ANTIGEN)]),
        ([("heavy", VH.lower()), ("light", None)], [("full", ANTIGEN)]),  # same sequence; a nanobody-like row
        ([("heavy", "EVQL;VESGG"), ("light", "")], [("full", ANTIGEN)]),  # not a sequence
        ([("heavy", None), ("light", VL)], [("full", "")]),
    ]
    pq.write_table(_table(rows), tmp_path / "part-00000-x-c000.snappy.parquet")
    return tmp_path


def test_antibody_fasta_has_distinct_chains_by_role(asd_dir):
    fasta = load_asd_fasta(asd_dir, "antibody")
    assert sorted(fasta.values()) == sorted([VH, VL])
    heavy = [h for h, s in fasta.items() if s == VH][0]
    assert heavy.startswith("ASDH_") and heavy.endswith("| polypeptide(L) | Auth:H")
    light = [h for h, s in fasta.items() if s == VL][0]
    assert light.startswith("ASDL_") and light.endswith("Auth:L")


def test_antigen_fasta_skips_empty_chains(asd_dir):
    fasta = load_asd_fasta(asd_dir, "antigen")
    assert list(fasta.values()) == [ANTIGEN]
    assert next(iter(fasta)).endswith("Auth:A")


def test_fasta_order_is_deterministic(asd_dir):
    assert list(load_asd_fasta(asd_dir, "antibody")) == list(load_asd_fasta(asd_dir, "antibody"))


def test_batched_cdr_extraction_matches_one_call_per_sequence():
    batch = _extract_H3L3_batch("H3", [VH, "", "GGGGGGGG", VH], ncpu=2)
    assert batch[0] == batch[3] == _extract_H3L3_sequence("H3", VH) != ""
    assert batch[1] == "" and batch[2] == ""
    assert _extract_H3L3_batch("L3", [VL], ncpu=1)[0] == _extract_H3L3_sequence("L3", VL) != ""


def test_antibody_split_reads_the_role_from_the_header(tmp_path):
    fasta = tmp_path / "ab.fasta"
    fasta.write_text(f">ASDH_1 | polypeptide(L) | Auth:H\n{VH}\n>ASDL_1 | polypeptide(L) | Auth:L\n{VL}\n"
                     f">ASDL_2 | polypeptide(L) | Auth:L\n{VH}\n")  # VH also seen as a light chain: a heavy
    seq_map = tmp_path / "map.tsv"
    seq_map.write_text(f"P00000000000000000001\t{VH}\nP00000000000000000002\t{VL}\nP00000000000000000003\tAAAA\n")
    out = separate_antibody_sequences(tmp_path / "tmp", seq_map, fasta)
    assert out["antibody_H"].read_text() == f">P00000000000000000001\n{VH}\n"
    assert out["antibody_L"].read_text() == f">P00000000000000000002\n{VL}\n"
    assert out["protein"].read_text() == ""


def test_antibody_split_fails_on_a_sequence_without_an_id(tmp_path):
    fasta = tmp_path / "ab.fasta"
    fasta.write_text(f">ASDH_1 | polypeptide(L) | Auth:H\n{VH}\n")
    seq_map = tmp_path / "map.tsv"
    seq_map.write_text(f"P00000000000000000002\t{VL}\n")
    with pytest.raises(ValueError, match="no seq_id"):
        separate_antibody_sequences(tmp_path / "tmp", seq_map, fasta)


def test_asd_cluster_ids_carry_their_db_code():
    from structcooker.utils.mapping import cluster_id

    assert cluster_id("ABG", "P00000000000019613402") == "cABG_P00000000000019613402"
    assert cluster_id("ABA", "A00000000000019613402") == "cABA_A00000000000019613402"


def test_subset_seq_id_map_by_length(tmp_path):
    from structcooker.instructions.transforms.sequence import subset_seq_id_map

    fasta = tmp_path / "ag.fasta"
    seqs = {"a": "A" * 5, "b": "C" * 20, "c": "D" * 40}
    fasta.write_text("".join(f">{k} | polypeptide(L) | Auth:A\n{v}\n" for k, v in seqs.items()))
    seq_map = tmp_path / "map.tsv"
    seq_map.write_text("".join(f"P{i:020d}\t{v}\n" for i, v in enumerate(seqs.values())))
    assert list(subset_seq_id_map(fasta, seq_map).values()) == list(seqs.values())
    assert list(subset_seq_id_map(fasta, seq_map, min_length=10, max_length=30).values()) == ["C" * 20]


def test_cdhit_at_one_clusters_identical_sequences_only(tmp_path):
    import shutil

    from structcooker.instructions.transforms.sequence import parse_cdhit_cluster, run_CDHIT

    if shutil.which("cd-hit") is None:
        pytest.skip("cd-hit not installed")
    fasta = tmp_path / "h3.fasta"
    # b is a piece of a (and 90% of its length); c is a itself; d differs by one residue
    fasta.write_text(">a\nARDYYGSSYWYFDV\n>b\nARDYYGSSYWYFD\n>c\nARDYYGSSYWYFDV\n>d\nARDYYGSSYWYFDL\n")
    out = tmp_path / "out.dat"
    run_CDHIT(fasta, out, seq_id=1.0)
    clusters = sorted(sorted(v) for v in parse_cdhit_cluster(tmp_path / "out.dat.clstr").values())
    assert clusters == [["a", "c"], ["b"], ["d"]]
