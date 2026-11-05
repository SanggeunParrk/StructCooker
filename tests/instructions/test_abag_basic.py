"""The experimental-structure rule (a point mutant of a crystallised antibody is not removed) and the polyprotein name rule."""

import pytest

from structcooker.instructions.readers.abag_basic import POLYPROTEIN, PdbComplexIndex

H = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
L = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
AG = "TQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQGYVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPVTAEDGTQRCEKCSKPCARVCYGLGMEHLREVRAV"
CH1 = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAAVQSSG"


@pytest.fixture
def index(tmp_path):
    summary = tmp_path / "sabdab.tsv"
    summary.write_text(
        "pdb\tHchain\tLchain\tmodel\tantigen_chain\tantigen_type\n1abc\tH\tL\t0\tA\tprotein\n2def\tH\tNA\t0\tA\tprotein\n3ghi\tH\tL\t0\tA\tprotein\n"
    )
    fasta = tmp_path / "pdb.fasta"
    fasta.write_text(
        f">1ABC_A_. | polypeptide(L) | Auth:H\n{H + CH1}\n>1ABC_B_. | polypeptide(L) | Auth:L\n{L}RTVAAPSVFIFPPS\n"
        f">1ABC_C_. | polypeptide(L) | Auth:A\n{AG[20:110]}\n"
        f">2DEF_A_. | polypeptide(L) | Auth:H\n{H}\n>2DEF_B_. | polypeptide(L) | Auth:A\n{AG[20:110]}\n"
        f">3GHI_A_. | polypeptide(L) | Auth:H\n{H}\n>3GHI_B_. | polypeptide(L) | Auth:L\n{L}\n>3GHI_C_. | polypeptide(L) | Auth:A\n{AG}\n",
    )
    return PdbComplexIndex(summary, fasta)


def test_exact_needs_every_chain_identical(index):
    assert (
        index.find(H, L, AG) == "3GHI"
    )  # heavy, light and antigen are the chains of 3GHI
    assert (
        index.find(H, L, AG[20:110]) is None
    )  # 1ABC has a constant domain on its chains, the antigen here is a fragment of 3GHI's: not identical
    assert (
        index.find(H, "", AG[20:110]) == "2DEF"
    )  # a single-chain pair against an entry without a light chain


def test_a_point_mutant_does_not(index):
    mutant = H.replace("DYWGQG", "DYAGQG")
    assert (
        index.find(mutant, L, AG) is None
    )  # the parent has a structure, the variant has none
    assert index.find(H, L[:50] + "G" + L[51:], AG) is None
    assert index.find(mutant, L, AG, "contained") is None


def test_another_antigen_does_not(index):
    assert index.find(H, L, "MKTAYIAKQR" * 20) is None
    assert index.find(H, L, "MKTAYIAKQR" * 20, "contained") is None


def test_contained_also_takes_constant_domains_and_fragments(index):
    assert index.find(H, L, AG, "contained") in ("1ABC", "3GHI")
    assert (
        index.find(H, L, AG[20:110], "contained") == "1ABC"
    )  # variable domain against a Fab chain, an antigen fragment
    assert index.find(H, "", AG[20:110], "contained") == "2DEF"


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("Genome polyprotein", True),
        ("ORF1ab polyprotein", True),
        ("Gag-Pol polyprotein", True),
        ("pp1ab", True),
        ("envelope glycoprotein gp120", False),
        ("Spike glycoprotein", False),
        ("poly(A) polymerase", False),
    ],
)
def test_polyprotein_names(name, expected):
    assert bool(POLYPROTEIN.search(name)) is expected
