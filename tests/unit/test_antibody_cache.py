import pytest

from structcooker.instructions.unit.antibody import (
    AntibodyClassifier,
    annotate_antibodies,
    unique_sequences,
)

HEAVY, LIGHT, OTHER = "H" * 80, "L" * 90, "O" * 100


def fake_classify(chunk):
    return [(s, {"H": "H", "L": "K"}.get(s[0])) for s in chunk]


@pytest.fixture
def fasta(tmp_path):
    p = tmp_path / "p.fasta"
    p.write_text(
        f">1\n{HEAVY}\n>2\n{LIGHT}\n>3\n{HEAVY}\n>4\n{OTHER}\n>5\n{'S' * 20}\n>6\n{'Z' * 700}\n"
    )
    return p


def test_unique_sequences_are_sorted_and_length_filtered(fasta):
    assert unique_sequences(fasta, 70, 600) == sorted([HEAVY, LIGHT, OTHER])


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_annotate_writes_one_sorted_line_per_sequence(fasta, tmp_path, n_jobs):
    out = tmp_path / "x" / "cache.tsv"
    status = annotate_antibodies(
        fasta, out, n_jobs=n_jobs, chunk_size=1, classify=fake_classify
    )
    assert status == "2 of 3 sequences classified as antibody/TCR"
    assert out.read_text().splitlines() == [
        f"{HEAVY}\tH",
        f"{LIGHT}\tK",
        f"{OTHER}\t",
    ]
    assert not list(out.parent.glob("*.part"))


def test_the_cache_is_what_the_classifier_reads(fasta, tmp_path):
    out = tmp_path / "cache.tsv"
    annotate_antibodies(fasta, out, classify=fake_classify)
    ab = AntibodyClassifier(out, lazy=False)
    assert (ab(HEAVY), ab(LIGHT), ab(OTHER)) == ("H", "K", None)
    assert ab("A" * 80) is None  # a frozen cache: unknown means not an antibody
