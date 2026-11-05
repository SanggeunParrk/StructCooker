"""Seq-metadata maps, SignalP loading, FASTA/TSV helpers and cif metadata extraction."""

import pickle
from types import SimpleNamespace

import numpy as np
import pytest

from structcooker.instructions.transforms import metadata as md


def test_lmdb_dict_roundtrip_pickles_to_path_only():
    d = md.seq_metadata_map_to_lmdb({"1abc_A_.": ("P1", "cPDB_P1")})
    assert "1abc_A_." in d
    assert "1abc_B_." not in d
    assert d["1abc_A_."] == ("P1", "cPDB_P1")
    with pytest.raises(KeyError):
        d["1abc_B_."]
    clone = pickle.loads(pickle.dumps(d))  # noqa: S301
    assert clone._env is None  # env handle is per-process, never pickled
    assert clone["1abc_A_."] == ("P1", "cPDB_P1")


def test_load_tsv_splits_on_comma_unless_disabled(tmp_path):
    p = tmp_path / "a.tsv"
    p.write_text("k1\ta,b\nk2\tc\n")
    assert md.load_tsv(p) == {"k1": ["a", "b"], "k2": ["c"]}
    assert md.load_tsv(p, split_by_comma=False) == {"k1": ["a,b"], "k2": ["c"]}


def test_reverse_dict_maps_members_to_their_key():
    assert md.reverse_dict({"r1": ["a", "b"], "r2": ["c"]}) == {
        "a": "r1",
        "b": "r1",
        "c": "r2",
    }


def test_load_pairs_keeps_repeated_left_keys(tmp_path):
    p = tmp_path / "p.tsv"
    p.write_text("a\t1\na\t2\nb\t3\n")
    assert md.load_pairs(p) == {"a": ["1", "2"], "b": ["3"]}


def test_build_seqid_map_is_keyed_by_mol_identifier_then_sequence():
    m = md.build_seqid_map({"P1": ["MK"], "R2": ["ACGU"]})
    assert m["P"] == {"MK": "P1"}
    assert m["R"] == {"ACGU": "R2"}
    assert m["D"] == {}


def test_build_seqid_map_rejects_unknown_identifier():
    with pytest.raises(ValueError, match="Unknown molecule identifier Z"):
        md.build_seqid_map({"Z1": ["MK"]})


def _fasta():
    return {
        "1abc_A_. | polypeptide(L) | Auth:A": "MK",
        "1abc_B_. | polyribonucleotide | Auth:B": "ACGU",
    }


def test_build_seq_metadata_map_uses_cluster_or_singleton_fallback():
    out = md.build_seq_metadata_map(
        _fasta(),
        {"P1": ["MK"], "R2": ["ACGU"]},
        {"cPDB_P1": ["P1"]},
        "PDB",
    )
    assert out == {
        "1abc_A_.": ("P1", "cPDB_P1"),
        # R2 is in no cluster -> its own singleton, named after the DB
        "1abc_B_.": ("R2", "cPDB_R2"),
    }


def test_build_seq_metadata_map_errors_on_sequence_without_seq_id():
    with pytest.raises(KeyError, match="Sequence ID not found for molecule type P"):
        md.build_seq_metadata_map(_fasta(), {"R2": ["ACGU"]}, {}, "PDB")


def test_classify_seq_clusters_returns_cluster_set():
    seqid_map = md.build_seqid_map({"P1": ["MK"], "R2": ["ACGU"]})
    got = md.classify_seq_clusters(
        _fasta(),
        {"1abc_A_. | polypeptide(L) | Auth:A": ""},
        seqid_map,
        {"cX": ["P1"]},
        "PDB",
    )
    assert got == {"cX"}
    got = md.classify_seq_clusters(_fasta(), _fasta(), seqid_map, {}, "PDB")
    assert got == {"cPDB_P1", "cPDB_R2"}
    with pytest.raises(KeyError, match="Sequence ID not found"):
        md.classify_seq_clusters(_fasta(), _fasta(), {"P": {}, "R": {}}, {}, "PDB")


def test_parse_signalp_is_zero_based_and_none_when_missing(tmp_path):
    f = tmp_path / "P1.gff3"
    f.write_text("##gff-version 3\nP1\tSignalP\tsignal_peptide\t1\t22\t0.9\t.\t.\t.\n")
    assert md.parse_signalp(f) == (0, 21)
    assert md.parse_signalp(tmp_path / "missing.gff3") is None


def test_load_signalp_reads_gff3_by_stem_and_validates_dir(tmp_path):
    (tmp_path / "P1.gff3").write_text("##\nP1\ts\tt\t3\t10\t.\t.\t.\t.\n")
    (tmp_path / "ignored.txt").write_text("x")
    assert md.load_signalp(tmp_path) == {"P1": (2, 9)}
    assert md.load_signalp(None) == {}
    with pytest.raises(FileNotFoundError, match="does not exist"):
        md.load_signalp(tmp_path / "nope")


def test_index_signalp_sequences_empty_dict_is_empty(tmp_path):
    assert md.index_signalp_sequences(tmp_path / "unused", {}) == {}


def test_index_signalp_warns_when_nothing_matches(tmp_path, caplog):
    p = tmp_path / "m.tsv"
    p.write_text("P9\tMK\n")
    with caplog.at_level("WARNING"):
        assert md.index_signalp_sequences(p, {"P1": (0, 5)}) == {}
    assert "No SignalP prediction IDs matched" in caplog.text


def test_load_fasta_joins_wrapped_sequences_and_rejects_headerless(tmp_path):
    p = tmp_path / "a.fa"
    p.write_text(">h1\nMK\nVL\n>h2\nGG\n")
    assert md.load_fasta(p) == {"h1": "MKVL", "h2": "GG"}
    bad = tmp_path / "b.fa"
    bad.write_text("MK\n>h\nA\n")
    with pytest.raises(ValueError, match="before any header"):
        md.load_fasta(bad)


def test_extract_protein_seqs_filters_nonprotein_and_all_unknown():
    out = md.extract_protein_seqs(
        {"P1": ["MK"], "P2": ["XXX"], "R3": ["ACGU"], "P4": ["MXK"]}
    )
    assert [o["seqid"] for o in out] == ["P1", "P4"]
    keep = md.extract_protein_seqs({"P2": ["XXX"]}, remove_unknown=False)
    assert keep == [{"seqid": "P2", "sequence": "XXX"}]


@pytest.mark.parametrize("seqs", [[], ["A", "B"]])
def test_extract_protein_seqs_requires_exactly_one_sequence(seqs):
    with pytest.raises(ValueError, match="sequence"):
        md.extract_protein_seqs({"P1": seqs})


def test_chunk_protein_seqs_splits_with_running_index():
    seqs = {f"P{i}": ["MK"] for i in range(5)}
    chunks = md.chunk_protein_seqs(seqs, chunk_size=2)
    assert [c["chunk_index"] for c in chunks] == [0, 1, 2]
    assert [len(c["chunk"]) for c in chunks] == [2, 2, 1]
    assert chunks[2]["chunk"] == [("P4", "MK")]
    assert md.chunk_protein_seqs({}, chunk_size=2) == []


def test_chain_template_tsv_loaders(tmp_path):
    ct = tmp_path / "ct.tsv"
    ct.write_text("1A_A\tT1,T2\n1B_A\t\nshort\n")
    assert md.load_chain_templates(ct) == {"1A_A": ["T1", "T2"]}
    sc = tmp_path / "sc.tsv"
    sc.write_text("P1\tc1,c2\nP2\t\n")
    assert md.invert_seqid_to_chains(sc) == {"c1": "P1", "c2": "P1"}
    ss = tmp_path / "ss.tsv"
    ss.write_text("P1\tMK\nP2\t\nP3\n")
    assert md.load_seq_tsv(ss) == {"P1": "MK"}


def _f(v):
    return SimpleNamespace(value=np.asarray(v))


class _Chains:
    def __init__(self, entity_types):
        self.entity_type = _f(entity_types)

    def __len__(self):
        return len(self.entity_type.value)


def _fake_cifmol(entity_types, metadata):
    return SimpleNamespace(
        metadata=metadata,
        chains=_Chains(entity_types),
        residues=[0] * 7,
        atoms=[0] * 40,
    )


def test_extract_metadata_summarises_counts_and_flags():
    mol = _fake_cifmol(
        ["polypeptide(L)", "polyribonucleotide"],
        {"resolution": 2.5, "deposition_date": "2000-01-01"},
    )
    out = md.extract_metadata({"1_1_.": {"cifmol": mol}})
    assert out["1_1_."] == {
        "resolution": "2.5",
        "deposition_date": "2000-01-01",
        "release_date": "NA",
        "chain_num": "2",
        "residue_num": "7",
        "atom_num": "40",
        "including_NA": "Yes",
        "including_Dform": "No",
    }
    d_mol = _fake_cifmol(["polypeptide(D)", "polypeptide(L)"], {})
    flags = md.extract_metadata({"k": {"cifmol": d_mol}})["k"]
    assert (flags["including_NA"], flags["including_Dform"]) == ("No", "Yes")


def _attach_mol(chain_ids):
    mol = SimpleNamespace(
        id=["1abc"],
        alt_id=".",
        chains=SimpleNamespace(chain_id=_f(chain_ids)),
    )
    mol.to_dict = lambda: {"chains": {"nodes": {}}}
    return mol


def test_attach_metadata_writes_seq_id_and_cluster_per_chain(monkeypatch):
    mol = _attach_mol(["A_1", "B_1"])
    captured = {}

    class Attached:
        @staticmethod
        def from_dict(d):
            captured["d"] = d
            return SimpleNamespace(to_dict=lambda: d)

    monkeypatch.setattr(md, "filter_water", lambda m: m)
    monkeypatch.setattr(md, "CIFMolAttached", Attached)
    out = md.attach_metadata(
        mol,
        {"1abc_A_.": ("P1", "cPDB_P1"), "1abc_B_.": ("P2", "cPDB_P2")},
    )
    nodes = out["chains"]["nodes"]
    assert nodes["seq_id"]["value"].tolist() == ["P1", "P2"]
    assert nodes["cluster_id"]["value"].tolist() == ["cPDB_P1", "cPDB_P2"]


def test_attach_metadata_all_water_is_skipped_and_missing_key_raises(monkeypatch):
    monkeypatch.setattr(md, "filter_water", lambda _m: None)
    assert md.attach_metadata(_attach_mol(["A"]), {}) is None
    monkeypatch.setattr(md, "filter_water", lambda m: m)
    with pytest.raises(KeyError, match="1abc_A_."):
        md.attach_metadata(_attach_mol(["A_1"]), {})


def test_parse_signalp_header_only_file_has_no_signal_peptide(tmp_path):
    from structcooker.instructions.transforms.metadata import parse_signalp

    path = tmp_path / "x.gff3"
    path.write_text("##gff-version 3\n")
    assert parse_signalp(path) is None
