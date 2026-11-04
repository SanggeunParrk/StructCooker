"""ab_seq / cdr_region: the SAbDab antibody rule on a fasta, and Chothia CDRs from ANARCI."""
import pytest

from structcooker.instructions.transforms.antibody_fasta import (
    cdr_regions_of,
    extract_cdr_regions,
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
    path.write_text(header + "".join(f"{p}\t{h}\t{l}\t0\n" for p, h, l in rows))


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


def test_cdr_regions_select_by_chain_type_and_name_a_second_domain():
    numbered = [(i, aa) for i, aa in zip(range(20, 40), "ABCDEFGHIJKLMNOPQRST", strict=True)]
    # scFv: a heavy domain then a light one, then a second heavy domain
    domains = [("H", numbered), ("L", numbered), ("H", numbered)]
    ranges = {"H1": [26, 27], "H2": [38, 39], "H3": [99], "L1": [30]}
    regions = dict(cdr_regions_of(domains, cdrs=("H1", "H2", "H3", "L1"), ranges=ranges))
    assert regions == {"H1": "GH", "H2": "ST", "L1": "K", "H1.2": "GH", "H2.2": "ST"}
    assert "H3" not in regions                       # no residue 99: the CDR is left out


def test_anarci_chothia_cdrs_of_trastuzumab(tmp_path):
    _fasta(tmp_path / "ab.fasta", [("1N8Z_B_. | polypeptide(L) | Auth:B", VH), ("1N8Z_A_. | polypeptide(L) | Auth:A", VL)])
    extract_cdr_regions(tmp_path / "ab.fasta", tmp_path / "cdr.fasta", batch_size=1, ncpu=1)
    got = {(h.split("|")[0].strip(), h.split("CDR:")[1]): s for h, s in read_fasta(tmp_path / "cdr.fasta")}
    expected = {("1N8Z_B_.", c): CDRS[c] for c in ("H1", "H2", "H3")} | {("1N8Z_A_.", c): CDRS[c] for c in ("L1", "L2", "L3")}
    assert got == expected


def test_extract_writes_one_record_per_chain_and_cdr_and_lists_failures(tmp_path):
    records = [
        ("1N8Z_B_. | polypeptide(L) | Auth:B", VH),
        ("1N8Z_B_A | polypeptide(L) | Auth:B", VH),            # the same sequence under another altloc
        ("1N8Z_A_. | polypeptide(L) | Auth:A", VL),
        ("9XXX_A_. | polypeptide(L) | Auth:A", OTHER),         # ANARCI finds no antibody domain
    ]
    _fasta(tmp_path / "ab.fasta", records)
    message = extract_cdr_regions(tmp_path / "ab.fasta", tmp_path / "cdr.fasta", cdrs=["H3", "L3"], batch_size=2, ncpu=1,
                                  failed_path=tmp_path / "failed.txt")
    out = list(read_fasta(tmp_path / "cdr.fasta"))
    assert out == [
        ("1N8Z_B_. | polypeptide(L) | Auth:B | CDR:H3", CDRS["H3"]),
        ("1N8Z_B_A | polypeptide(L) | Auth:B | CDR:H3", CDRS["H3"]),
        ("1N8Z_A_. | polypeptide(L) | Auth:A | CDR:L3", CDRS["L3"]),
    ]
    assert (tmp_path / "failed.txt").read_text() == "9XXX_A_. | polypeptide(L) | Auth:A\n"
    assert "3 of 4 chains" in message


def test_unknown_cdr_name_is_refused(tmp_path):
    _fasta(tmp_path / "ab.fasta", [("1N8Z_B_. | polypeptide(L) | Auth:B", VH)])
    with pytest.raises(ValueError, match="unknown CDR names"):
        extract_cdr_regions(tmp_path / "ab.fasta", tmp_path / "cdr.fasta", cdrs=["H4"])
