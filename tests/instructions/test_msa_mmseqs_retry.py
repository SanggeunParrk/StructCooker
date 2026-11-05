"""A failed MMseqs2 batch is split and retried; only a query that fails alone raises."""

import pytest

from structcooker.instructions.transforms import msa_mmseqs
from structcooker.instructions.transforms.msa_mmseqs import run_mmseqs_msa_search


def _fake_search(max_ok, bad=()):
    calls = []

    def fake(pending, work_dir, output_dir, db_uniref, mmseqs_bin, threads):
        calls.append(len(pending))
        if len(pending) > max_ok or any(seqid in bad for seqid, _ in pending):
            msg = "expandaln failed (-11)"
            raise RuntimeError(msg)
        for seqid, _ in pending:
            path = msa_mmseqs._a3m_path(output_dir, seqid)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(">q\nAAA\n")
        return len(pending)

    return fake, calls


def _run(tmp_path, n):
    chunk = [(f"P{i:020d}", "ACDEFGHIK") for i in range(n)]
    return run_mmseqs_msa_search(
        chunk,
        0,
        tmp_path / "work",
        tmp_path / "out",
        db_uniref=tmp_path,
        mmseqs_bin=tmp_path,
    )


def test_a_failed_batch_is_split_until_it_passes(tmp_path, monkeypatch):
    fake, calls = _fake_search(max_ok=5)
    monkeypatch.setattr(msa_mmseqs, "_search_chunk", fake)
    result = _run(tmp_path, 25)
    assert result["written"] == 25
    assert result["searched"] == 25
    assert calls[0] == 25
    assert max(calls[1:]) <= 7
    assert len(list((tmp_path / "out").rglob("*.a3m"))) == 25


def test_a_query_that_fails_alone_raises(tmp_path, monkeypatch):
    fake, _ = _fake_search(max_ok=5, bad={f"P{3:020d}"})
    monkeypatch.setattr(msa_mmseqs, "_search_chunk", fake)
    with pytest.raises(RuntimeError, match="expandaln"):
        _run(tmp_path, 10)


def test_a_finished_batch_is_not_searched_again(tmp_path, monkeypatch):
    fake, calls = _fake_search(max_ok=100)
    monkeypatch.setattr(msa_mmseqs, "_search_chunk", fake)
    _run(tmp_path, 10)
    assert _run(tmp_path, 10)["searched"] == 0
    assert calls == [10]
