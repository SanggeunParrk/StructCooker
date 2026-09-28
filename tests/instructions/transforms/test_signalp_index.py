from structcooker.instructions.transforms.metadata import index_signalp_sequences


def test_signalp_index_keeps_molecule_identity_and_omits_unpredicted(tmp_path):
    table = tmp_path / "sequences.tsv"
    table.write_text("P1\tABC\nR2\tABC\nP3\tDEF\n")
    assert index_signalp_sequences(table, {"P1": (0, 2), "R2": (1, 1)}) == {
        ("P", "ABC"): (0, 2), ("R", "ABC"): (1, 1),
    }


def test_signalp_duplicate_sequence_uses_final_id(tmp_path):
    table = tmp_path / "sequences.tsv"
    table.write_text("P1\tABC\nP2\tABC\n")
    assert index_signalp_sequences(table, {"P1": (0, 2)}) == {}
    assert index_signalp_sequences(table, {"P2": (0, 3)}) == {("P", "ABC"): (0, 3)}


def test_signalp_counter_padding_does_not_change_identity(tmp_path):
    table = tmp_path / "sequences.tsv"
    table.write_text("P00000000000000000012\tABC\nR00000000000000000012\tDEF\n")
    assert index_signalp_sequences(table, {"P0000012": (0, 3)}) == {("P", "ABC"): (0, 3)}


def test_signalp_conflicting_counter_aliases_are_rejected(tmp_path):
    import pytest
    with pytest.raises(ValueError, match="Conflicting SignalP"):
        index_signalp_sequences(tmp_path / "unused", {"P1": (0, 1), "P0001": (0, 2)})
