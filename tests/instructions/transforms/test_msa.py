import numpy as np
import pytest
from datacooker.api import execute

from structcooker.instructions.readers.a3m import get_a3m_data
from structcooker.instructions.transforms.msa import parse_headers, parse_sequence
from structcooker.instructions.transforms.openfold import (
    adapt_msa_for_cap,
    build_msa_features,
    cap_msa_dict,
)
from structcooker.utils.mapping import ResidueMapping
from structcooker.workflows.exports.cap_msa_depth import RECIPE as CAP_RECIPE


def test_insertions_and_profile_match_independent_calculation():
    rows = ["ACD-", "aaAbCcD-dd", "AC" + "z" * 300 + "D-"]
    result = parse_sequence(rows)
    expected = np.array([[0, 0, 0, 0], [2, 1, 1, 0], [0, 0, 255, 0]], dtype=np.int32)
    np.testing.assert_array_equal(result["deletions"], expected)
    aligned = result["aligned_sequences"]
    np.testing.assert_array_equal(aligned[0], aligned[1])
    classes = ResidueMapping().MAX_INDEX + 1
    reference = np.eye(classes, dtype=np.int32)[aligned].mean(axis=0).astype(np.float32)
    np.testing.assert_array_equal(result["profile"], reference)
    deletion_mean = (2 * np.arctan(expected.astype(np.float32) / 3) / np.pi).mean(axis=0)
    np.testing.assert_array_equal(result["deletion_mean"], deletion_mean)


def test_query_insertions_do_not_change_alignment_width():
    result = parse_sequence(["aACd", "AC"])
    assert result["aligned_sequences"].shape == (2, 2)
    np.testing.assert_array_equal(result["deletions"], [[1, 0], [0, 0]])


@pytest.mark.parametrize("rows", [[], [""], ["AC", "A"], ["AC", "A."]])
def test_invalid_alignment_fails_clearly(rows):
    with pytest.raises(ValueError, match="MSA"):
        parse_sequence(rows)


def test_header_rows_do_not_inherit_previous_hit():
    headers = ["query", "UniRef100_ABC protein Tax=Homo sapiens TaxID=9606 RepID=ABC_HUMAN", "SRR4029434_2280741", "other_hit"]
    result = parse_headers(headers)
    assert result["database"].tolist() == [b"query", b"uniref100", b"bfd", b"bfd"]
    assert result["database_id"].tolist() == [b"query", b"ABC", b"SRR4029434_2280741", b"other_hit"]
    assert result["species"].tolist() == [b"query", b"Homo sapiens", b"N/A", b"N/A"]


def test_reader_handles_blank_lines_and_wrapped_sequences(tmp_path):
    source = tmp_path / "query.a3m"
    source.write_text("\n>query\nAC\nD\n\n>hit\nAaCDxx\n")
    data = get_a3m_data(source)
    assert data == {"headers": ["query", "hit"], "raw_sequences": ["ACD", "AaCDxx"]}
    assert parse_sequence(data["raw_sequences"])["aligned_sequences"].shape == (2, 3)


def test_native_cap_recipe_preserves_query_and_recomputes_statistics():
    rows = ["ACD", "AaCD", "---"]
    source = {"sequences": parse_sequence(rows), "headers": parse_headers(["q", "hit1", "hit2"])}
    result = execute(CAP_RECIPE, {**adapt_msa_for_cap({"msa_dict": source}), "max_depth": 2}, targets=["msa_dict"])["msa_dict"]
    reference = parse_sequence(rows[:2])
    for key, value in reference.items():
        np.testing.assert_array_equal(result["sequences"][key], value)
    assert result["headers"]["database_id"].tolist() == [b"query", b"hit1"]
    assert source["sequences"]["aligned_sequences"].shape == (3, 3)


@pytest.mark.parametrize("depth", [0, -1])
def test_cap_cannot_drop_query(depth):
    with pytest.raises(ValueError, match="query"):
        cap_msa_dict({}, depth)


@pytest.mark.parametrize(("kind", "rows"), [("protein", ["ACD", "AaC-D".replace("-", "")]), ("rna", ["ACGU", "AaCGU"])])
def test_openfold_features_match_a3m(kind, rows):
    parsed = parse_sequence(rows, kind)
    matrix = np.array([list("".join(c for c in row if not c.islower())) for row in rows])
    direct = build_msa_features(matrix, parsed["deletions"], kind)
    for key, value in parsed.items():
        np.testing.assert_array_equal(direct[key], value)
