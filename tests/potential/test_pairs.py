from structcooker.instructions.potential import pairs


def test_type_table_is_28_plus_unknowns():
    assert pairs.T == 31
    assert len(set(pairs.TYPES)) == 31
    assert {"X", "DN", "RN"} <= set(pairs.TYPES)
    assert pairs.NBIN == 40
    assert pairs.BIN == 0.5
