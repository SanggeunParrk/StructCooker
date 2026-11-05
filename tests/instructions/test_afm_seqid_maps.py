"""Seq_id maps for monomer distillation entries and AFM MSA entities."""

import pytest

from structcooker.instructions.transforms import afm_msa_seqid as afm
from structcooker.instructions.transforms import monomer_seqid as mono


def test_entry_of_strips_tail_and_chain_fields_keeping_underscored_accessions():
    assert mono._entry_of("MGYP001_A_1 | protein | Auth:A") == "MGYP001"
    assert mono._entry_of("RF_00_x_A_1 | rna | Auth:A") == "RF_00_x"


def test_build_monomer_seqid_map_joins_fastas_and_skips_unknown_sequences(tmp_path):
    ids = tmp_path / "ids.tsv"
    ids.write_text("P1\tMKV\nP2\tGGG\nbroken\n")
    f1, f2 = tmp_path / "a.fa", tmp_path / "b.fa"
    f1.write_text(">E1_A_1 | p | Auth:A\nMKV\n\n>E2_A_1 | p | Auth:A\nUNKNOWN\n")
    f2.write_text(">E3_A_1 | p | Auth:A\nGGG\n")
    assert mono.build_monomer_seqid_map([str(f1), str(f2)], ids) == {
        "E1": "P1",
        "E3": "P2",
    }


def _afm_inputs(tmp_path, *, hetero_row="AF-2\tB\tacc\tAF-9\tx\tok"):
    fasta = tmp_path / "afm.fa"
    fasta.write_text(
        ">AF-1-model_v1_A_0 | x\nMKV\n"
        ">AF-1-model_v1_B_0 | x\nMKV\n"
        ">AF-2-model_v1_A_0 | x\nGGG\n"
        ">AF-2-model_v1_B_0 | x\nLLL\n"
    )
    ids = tmp_path / "ids.tsv"
    ids.write_text("P1\tMKV\nP2\tGGG\nP3\tLLL\n")
    chain_msa = tmp_path / "chain_msa.tsv"
    chain_msa.write_text(
        "het\tchain\tacc\tmonomer\thow\tverified\n" + hetero_row + "\n"
    )
    return fasta, ids, chain_msa


def test_afm_map_uses_own_entity_for_homodimers_and_monomer_for_hetero(tmp_path):
    fasta, ids, chain_msa = _afm_inputs(tmp_path)
    # AF-2 is a heterodimer: its chain B uses monomer MSA AF-9; AF-1 is a homodimer
    assert afm.build_afm_msa_seqid_map(fasta, ids, chain_msa) == {
        "AF-1": "P1",
        "AF-9": "P3",
    }


def test_afm_map_ignores_unverified_hetero_chain(tmp_path):
    fasta, ids, chain_msa = _afm_inputs(
        tmp_path, hetero_row="AF-2\tB\tacc\tAF-9\tx\tmismatch"
    )
    assert afm.build_afm_msa_seqid_map(fasta, ids, chain_msa) == {"AF-1": "P1"}


def test_afm_map_rejects_entity_used_with_two_sequences(tmp_path):
    fasta, ids, chain_msa = _afm_inputs(
        tmp_path, hetero_row="AF-2\tB\tacc\tAF-1\tx\tok"
    )
    with pytest.raises(ValueError, match="AF-1 is used by chains with different"):
        afm.build_afm_msa_seqid_map(fasta, ids, chain_msa)


def test_afm_map_raises_when_sequence_has_no_seq_id(tmp_path):
    fasta, ids, chain_msa = _afm_inputs(tmp_path)
    ids.write_text("P1\tMKV\n")
    with pytest.raises(KeyError, match="MSA sequences have no seq_id"):
        afm.build_afm_msa_seqid_map(fasta, ids, chain_msa)
