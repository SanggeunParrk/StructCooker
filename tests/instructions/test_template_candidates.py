"""Template Phase 1/2 precompute: hit parsing, date filters, top-k, reduced hmm."""

import datetime as dt

from structcooker.instructions.transforms import template_candidates as tc


def _row(target: str, score: float = 100.0) -> str:
    return f"  1e-30  {score}  0.1  1e-30  {score}  0.1  1.0  1  {target}  | desc\n"


def _hmm(targets: list[str]) -> str:
    text = "# header\n  --- full sequence ---\n"
    text += "".join(_row(t) for t in targets)
    text += "\nDomain annotation\n"
    for t in targets:
        text += f">> {t} | desc\n   dom line for {t}\n"
    text += "Internal pipeline statistics\n//\n"
    return text


def test_reduce_id_keeps_pdbid_and_chain():
    assert tc._reduce_id("1B4U_A_.") == "1B4U_A"
    assert tc._reduce_id("1B4U") == "1B4U"


def test_parse_date_rejects_garbage_and_none():
    assert tc._parse_date(" 2020-01-02 ") == dt.date(2020, 1, 2)
    assert tc._parse_date("nope") is None
    assert tc._parse_date(None) is None  # type: ignore[arg-type]


def test_ordered_hits_dedups_and_keeps_order():
    text = _hmm(["1AAA_A_.", "2BBB_B_.", "1AAA_A_x"])
    assert tc._ordered_hits(text) == ["1AAA_A", "2BBB_B"]


def test_load_pdb_dates_prefers_release_then_deposition_then_positional(tmp_path):
    a = tmp_path / "a.tsv"
    a.write_text(
        "cif_id\tdeposition_date\trelease_date\n"
        "1ABC_1_0_0\t2000-01-01\t2001-01-01\n"
        "2ABC_1_0_0\t2000-01-01\tbad\n"
        "short\n"
    )
    assert tc._load_pdb_dates(a) == {"1abc": dt.date(2001, 1, 1)}

    b = tmp_path / "b.tsv"
    b.write_text("pdbid\tdeposition_date\n3XYZ\t1999-05-05\n")
    assert tc._load_pdb_dates(b) == {"3xyz": dt.date(1999, 5, 5)}

    c = tmp_path / "c.tsv"
    c.write_text("cif_id\tresolution\tfoo\n4ABC_x\t2.0\t2010-10-10\n")
    assert tc._load_pdb_dates(c) == {"4abc": dt.date(2010, 10, 10)}


def test_read_chain_list_skips_blank_lines(tmp_path):
    p = tmp_path / "c.txt"
    p.write_text("1ABC_A\n\n  \n2DEF_B\n")
    assert tc._read_chain_list(p) == ["1ABC_A", "2DEF_B"]


FASTA = (
    ">1ABC_A_1 | polypeptide(L)\nMKV\nLL\n"
    ">1ABC_B_1 | polypeptide(L)\nGGG\n"
    ">2DEF_A_1 | RNA\nACGU\n"
)


def test_parse_cif_fasta_joins_wrapped_lines_and_filters_wanted(tmp_path):
    p = tmp_path / "f.fa"
    p.write_text(FASTA)
    assert tc._parse_cif_fasta(p, {"1ABC_A", "2DEF_A"}) == {
        "1ABC_A": "MKVLL",
        "2DEF_A": "ACGU",
    }


def test_stream_seq_to_id_only_needed(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("s1\tAAA\ns2\tBBB\nbroken\n")
    assert tc._stream_seq_to_id(p, {"BBB"}) == {"BBB": "s2"}


def test_hmm_path_layout(tmp_path):
    assert tc._hmm_path(tmp_path, "s123456") == tmp_path / "s" / "456" / "s123456.out"


def test_reduce_hmm_text_drops_nonkept_rows_and_domain_blocks():
    out = tc._reduce_hmm_text(_hmm(["1AAA_A_.", "2BBB_B_."]), {"2BBB_B"})
    assert "2BBB_B_." in out
    assert "1AAA_A_." not in out
    assert "dom line for 2BBB_B_." in out
    assert "dom line for 1AAA_A_." not in out
    assert "# header" in out
    assert out.rstrip().endswith("//")


def _setup_hmm(tmp_path, seq_id, targets):
    path = tc._hmm_path(tmp_path / "hmm", seq_id)
    path.parent.mkdir(parents=True)
    path.write_text(_hmm(targets))
    return tmp_path / "hmm"


def test_process_seqid_applies_cutoff_daydiff_and_topk(tmp_path):
    hmm = _setup_hmm(
        tmp_path,
        "s000001",
        ["9NEW_A_.", "1OLD_A_.", "2OLD_A_.", "3OLD_A_.", "4UNK_A_."],
    )
    dates = {
        "q": dt.date(2020, 1, 1),
        "9new": dt.date(2022, 1, 1),  # after global cutoff
        "1old": dt.date(2019, 12, 15),  # only 17 days older than query
        "2old": dt.date(2010, 1, 1),
        "3old": dt.date(2011, 1, 1),
    }
    red = tmp_path / "red"
    res = tc._process_seqid(
        "s000001",
        ["Q_A", "ZZZ_A"],
        hmm,
        red,
        dates,
        dt.date(2021, 9, 30),
        60,
        1,
    )
    # unknown-date query chain gets nothing; topk=1 keeps the first passing hit
    assert res == [("Q_A", ["2OLD_A"]), ("ZZZ_A", [])]
    reduced = (red / "s" / "001" / "s000001.reduced.out").read_text()
    assert "2OLD_A_." in reduced
    assert "3OLD_A_." not in reduced


def test_process_seqid_missing_hmm_returns_empty_per_chain(tmp_path):
    res = tc._process_seqid(
        "s000009", ["A_A", "B_B"], tmp_path, tmp_path, {}, dt.date(2021, 1, 1), 60, 5
    )
    assert res == [("A_A", []), ("B_B", [])]


def _inputs(tmp_path):
    fasta = tmp_path / "cif.fa"
    fasta.write_text(
        ">1QRY_A_1 | polypeptide(L)\nMKV\n"
        ">1QRY_B_1 | polypeptide(L)\nMKV\n"
        ">2OLD_A_1 | polypeptide(L)\nGGG\n"
        ">3NA_A_1 | RNA\nACGU\n"
    )
    chains = tmp_path / "chains.txt"
    chains.write_text("1QRY_A\n1QRY_B\n")
    smap = tmp_path / "seq_id_map.tsv"
    smap.write_text("s000001\tMKV\ns000002\tGGG\n")
    dates = tmp_path / "dates.tsv"
    dates.write_text(
        "cif_id\tdeposition_date\trelease_date\n"
        "1QRY_1_0_0\t2020-01-01\t2020-06-01\n"
        "2OLD_1_0_0\t2005-01-01\t2005-06-01\n"
    )
    hmm = _setup_hmm(tmp_path, "s000001", ["2OLD_A_.", "1QRY_A_."])
    return fasta, chains, smap, dates, hmm


def test_precompute_candidates_end_to_end_and_phase1_reuse(tmp_path):
    fasta, chains, smap, dates, hmm = _inputs(tmp_path)
    sc = tmp_path / "out" / "seqid_chains.tsv"
    ct = tmp_path / "out" / "chain_templates.tsv"
    msg = tc.precompute_candidates(
        fasta, chains, smap, dates, hmm, sc, ct, tmp_path / "red", n_jobs=1
    )
    assert sc.read_text() == "s000001\t1QRY_A,1QRY_B\n"
    # the query itself (day diff 0) is excluded; the 2005 template stays
    assert ct.read_text() == "1QRY_A\t2OLD_A\n1QRY_B\t2OLD_A\n"
    assert msg.startswith("candidates: 2 chains, 2 with >=1 template")

    # a second run reuses the written seqid_chains file even if the inputs vanish
    chains.unlink()
    again = tc._phase1(fasta, chains, smap, sc)
    assert again == {"s000001": ["1QRY_A", "1QRY_B"]}


def test_precompute_seqs_writes_only_needed_and_protein_chains(tmp_path):
    fasta, _chains, smap, _dates, _hmm_dir = _inputs(tmp_path)
    sc = tmp_path / "sc.tsv"
    sc.write_text("s000002\tX_A\n")
    out_s, out_c = tmp_path / "ss.tsv", tmp_path / "cs.tsv"
    msg = tc.precompute_seqs(sc, smap, fasta, out_s, out_c)
    assert out_s.read_text() == "s000002\tGGG\n"
    assert out_c.read_text() == "1QRY_A\tMKV\n1QRY_B\tMKV\n2OLD_A\tGGG\n"
    assert msg == "seqs: seqid_to_seq=1, chain_to_seq=3"


def test_select_seqid_ignores_query_age_and_caps_candidates(tmp_path):
    hmm = _setup_hmm(tmp_path, "s000001", ["1AAA_A_.", "2BBB_A_.", "3CCC_A_."])
    dates = {
        "1aaa": dt.date(2000, 1, 1),
        "2bbb": dt.date(2030, 1, 1),
        "3ccc": dt.date(2001, 1, 1),
    }
    sid, kept = tc._select_seqid(
        "s000001", hmm, tmp_path / "r", dates, dt.date(2021, 9, 30), 5
    )
    assert (sid, kept) == ("s000001", ["1AAA_A", "3CCC_A"])
    _, capped = tc._select_seqid(
        "s000001", hmm, tmp_path / "r", dates, dt.date(2021, 9, 30), 1
    )
    assert capped == ["1AAA_A"]
    assert tc._select_seqid(
        "s999999", hmm, tmp_path, dates, dt.date(2021, 1, 1), 5
    ) == (
        "s999999",
        [],
    )


def test_precompute_seqid_candidates_writes_list_of_seqids_with_hits(tmp_path):
    fasta, _c, _s, dates, hmm = _inputs(tmp_path)
    ids = tmp_path / "ids.tsv"
    ids.write_text("s000001\tMKV\ns000002\tGGG\n\n")
    out_t, out_l = tmp_path / "t.tsv", tmp_path / "l.txt"
    out_cs = tmp_path / "cs.tsv"
    msg = tc.precompute_seqid_candidates(
        ids, dates, hmm, fasta, out_t, out_l, out_cs, tmp_path / "red", n_jobs=1
    )
    assert out_t.read_text() == "s000001\t2OLD_A,1QRY_A\ns000002\t\n"
    assert out_l.read_text() == "s000001\n"
    assert "2 seq_ids, 1 with >=1 template; chain_to_seq=3" in msg
    assert out_cs.read_text().splitlines()[0] == "1QRY_A\tMKV"
