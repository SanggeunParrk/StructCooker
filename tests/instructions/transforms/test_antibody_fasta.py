"""ab_seq / cdr_region: the SAbDab antibody rule on a fasta, and Chothia CDR positions from ANARCI."""
import pytest

from structcooker.instructions.transforms.antibody_fasta import (
    cdr_spans_of,
    extract_cdr_table,
    filter_fasta_antibody,
    read_fasta,
)

# Trastuzumab Fab variable domains.
VH = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC"
      "SRWGGDGFYAMDYWGQGTLVTVSS")
VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
OTHER = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"   # ubiquitin
CDRS = {"H1": "GFNIKDT", "H2": "RIYPTNGYTRYADSVKG", "H3": "WGGDGFYAMDY",
        "L1": "RASQDVNTAVA", "L2": "SASFLYS", "L3": "QQHYTTPPT"}


def _fasta(path, records):
    path.write_text("".join(f">{h}\n{s}\n" for h, s in records))


def _summary(path, rows):
    header = "pdb\tHchain\tLchain\tmodel\n"
    path.write_text(header + "".join(f"{p}\t{h}\t{light}\t0\n" for p, h, light in rows))


def test_filter_keeps_every_record_of_a_sabdab_sequence(tmp_path):
    records = [
        ("1ABC_H_. | polypeptide(L) | Auth:H", VH),
        ("1ABC_L_. | polypeptide(L) | Auth:L", VL),
        ("1ABC_A_. | polypeptide(L) | Auth:A", OTHER),          # the antigen of the entry
        ("1ABC_Z_. | non-polymer | Auth:H", "(SO4)"),
        ("2XYZ_B_. | polypeptide(L) | Auth:B", VH),               # not in SAbDab, but an antibody sequence
        ("2XYZ_C_. | polypeptide(L) | Auth:C", OTHER),
        ("3DEF_A_. | polyribonucleotide | Auth:A", "ACGU"),
    ]
    _fasta(tmp_path / "all.fasta", records)
    _summary(tmp_path / "sabdab.tsv", [("1abc", "H", "L"), ("4zzz", "NA", "Q")])
    message = filter_fasta_antibody(tmp_path / "all.fasta", tmp_path / "sabdab.tsv", tmp_path / "ab.fasta")
    assert list(read_fasta(tmp_path / "ab.fasta")) == [records[0], records[1], records[4]]   # order and headers kept
    assert "3 records" in message
    assert (tmp_path / "ab.fasta").read_text().startswith(">1ABC_H_. | polypeptide(L) | Auth:H\n" + VH)


def test_filter_needs_the_polypeptide_l_record_of_the_listed_chain(tmp_path):
    _fasta(tmp_path / "all.fasta", [("1ABC_H_. | polypeptide(D) | Auth:H", VH)])
    _summary(tmp_path / "sabdab.tsv", [("1abc", "H", "NA")])
    filter_fasta_antibody(tmp_path / "all.fasta", tmp_path / "sabdab.tsv", tmp_path / "ab.fasta")
    assert (tmp_path / "ab.fasta").read_text() == ""


def test_cdr_spans_are_first_and_last_position_and_none_for_an_absent_cdr():
    residues = [(i, 100 + i) for i in range(20, 40)]                  # Chothia number i sits at position 100 + i
    ranges = {"H1": [26, 27], "H2": [38, 39], "H3": [99]}
    assert cdr_spans_of(("H", residues), ranges) == [(126, 128), (138, 140), None]


def _write_map(path, sequences):
    path.write_text("".join(f"{i}\t{s}\n" for i, s in sequences))


def _rows(path):
    return [line.split(",") for line in path.read_text().splitlines()]


def test_cdr_table_positions_slice_the_cdrs_out_of_the_sequence(tmp_path):
    tagged = "XXMGSSHHHH" + VH                                    # a tag with X residues in front: positions must map back
    _fasta(tmp_path / "ab.fasta", [("1N8Z_B_. | polypeptide(L) | Auth:B", VH), ("1N8Z_B_A | polypeptide(L) | Auth:B", VH),
                                   ("1N8Z_A_. | polypeptide(L) | Auth:A", VL), ("2TAG_A_. | polypeptide(L) | Auth:A", tagged),
                                   ("9XXX_A_. | polypeptide(L) | Auth:A", OTHER)])
    _write_map(tmp_path / "map.tsv", [("P00000000000000000007", VH), ("P00000000000000000003", VL),
                                      ("P00000000000000000011", tagged), ("P00000000000000000013", OTHER),
                                      ("R00000000000000000001", "ACGU"), ("Q00000000000000000002", VH)])
    message = extract_cdr_table(tmp_path / "ab.fasta", tmp_path / "map.tsv", tmp_path / "cdr.csv", batch_size=2, ncpu=1,
                                failed_path=tmp_path / "failed.txt", extra_path=tmp_path / "extra.csv")
    rows = _rows(tmp_path / "cdr.csv")
    assert rows[0] == ["seqid", "cdr1_start", "cdr1_end", "cdr2_start", "cdr2_end", "cdr3_start", "cdr3_end"]
    assert [r[0] for r in rows[1:]] == ["P00000000000000000003", "P00000000000000000007", "P00000000000000000011"]   # sorted, once each
    by_id = {r[0]: [int(x) for x in r[1:]] for r in rows[1:]}
    sequences = {"P00000000000000000003": (VL, "L"), "P00000000000000000007": (VH, "H"), "P00000000000000000011": (tagged, "H")}
    for seq_id, (sequence, region) in sequences.items():
        spans = by_id[seq_id]
        got = [sequence[spans[2 * k] : spans[2 * k + 1]] for k in range(3)]
        assert got == [CDRS[f"{region}{n}"] for n in (1, 2, 3)], seq_id
    assert (tmp_path / "failed.txt").read_text() == "P00000000000000000013\n"
    assert _rows(tmp_path / "extra.csv") == [["seqid", "domain", "chain_type", *rows[0][1:]]]
    assert "3 sequences of 4" in message


def test_a_second_domain_goes_to_the_extra_file(tmp_path):
    scfv = VH + "GGGGSGGGGSGGGGS" + VL
    _fasta(tmp_path / "ab.fasta", [("3SCF_A_. | polypeptide(L) | Auth:A", scfv)])
    _write_map(tmp_path / "map.tsv", [("P00000000000000000005", scfv)])
    extract_cdr_table(tmp_path / "ab.fasta", tmp_path / "map.tsv", tmp_path / "cdr.csv", ncpu=1, extra_path=tmp_path / "extra.csv")
    first = [int(x) for x in _rows(tmp_path / "cdr.csv")[1][1:]]
    assert [scfv[first[2 * k] : first[2 * k + 1]] for k in range(3)] == [CDRS["H1"], CDRS["H2"], CDRS["H3"]]
    extra = _rows(tmp_path / "extra.csv")[1]
    assert extra[:3] == ["P00000000000000000005", "2", "L"]
    second = [int(x) for x in extra[3:]]
    assert [scfv[second[2 * k] : second[2 * k + 1]] for k in range(3)] == [CDRS["L1"], CDRS["L2"], CDRS["L3"]]


def test_a_sequence_without_a_seq_id_is_refused(tmp_path):
    _fasta(tmp_path / "ab.fasta", [("1N8Z_B_. | polypeptide(L) | Auth:B", VH)])
    _write_map(tmp_path / "map.tsv", [("P00000000000000000003", VL)])
    with pytest.raises(ValueError, match="no seq_id"):
        extract_cdr_table(tmp_path / "ab.fasta", tmp_path / "map.tsv", tmp_path / "cdr.csv", ncpu=1)
