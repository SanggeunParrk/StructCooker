from structcooker.potential import pairs


def test_type_table_is_28_plus_unknowns():
    assert pairs.T == 31 and len(set(pairs.TYPES)) == 31
    assert {"X", "DN", "RN"} <= set(pairs.TYPES)
    assert pairs.NBIN == 40 and pairs.BIN == 0.5
