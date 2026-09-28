from copy import deepcopy

import pytest

from structcooker.schemas import get


def filtered_record():
    return {"1_1_.": {"cifmol_dict": {
        "atoms": {}, "residues": {}, "chains": {}, "index_table": {},
        "metadata": {"id": ["1abc"]},
    }}}


def test_filtered_cif_has_a_distinct_schema():
    value = filtered_record()
    assert get("I").validate(value) == []
    assert get("A").validate(value)
    assert get("B").validate(value)
    assert get("I").expansion == get("A").expansion


@pytest.mark.parametrize("missing", ["atoms", "residues", "chains", "index_table", "metadata"])
def test_filtered_schema_checks_every_assembly(missing):
    value = filtered_record()
    value["2_1_."] = deepcopy(value["1_1_."])
    del value["2_1_."]["cifmol_dict"][missing]
    assert get("I").validate(value)


@pytest.mark.parametrize("value", [{}, None, {"1_1_.": {"cifmol_attached_dict": {}}}])
def test_filtered_schema_rejects_other_wrappers(value):
    assert get("I").validate(value)
