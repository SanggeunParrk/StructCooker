"""Chain status of an ASD pair and the light chains recovered from the AntiBinder tables."""

import pytest

from structcooker.instructions.readers.abag_chains import (
    AntiBinderIndex,
    classify_chains,
)

H, L1, L2, ANT = (
    "EVQLVESGGG" * 12,
    "DIQMTQSPSS" * 11,
    "QSALTQPPAV" * 11,
    "MKTAYIAKQR" * 5,
)


@pytest.fixture
def antibinder(tmp_path):
    for d in ("COVID-19", "HIV", "BioMap"):
        (tmp_path / d).mkdir()
    (tmp_path / "COVID-19/CoV-AbDab_only_sars2_filter.csv").write_text(
        f"Heavy,Light,Target,Protein_Epitope,Label,antigen\n{H},{L1},x,y,1,{ANT}\n"
    )
    (tmp_path / "HIV/dataset_hiv_cls.csv").write_text(
        f"antibody,virus_seq,ANT_Binding,split\n{H + L2},{ANT},1,train\n"
    )
    (tmp_path / "BioMap/biomap.csv").write_text(
        f"pdb,antibody_seq_a,antibody_seq_b,antigen_seq,delta_g\n1abc,{L1},{H},{ANT},-9\n"
    )
    return AntiBinderIndex(tmp_path)


def _status(ab, heavy=H, light="", antigen=ANT, **flags):
    return classify_chains(
        heavy,
        light,
        nanobody=flags.get("nanobody", False),
        scfv=flags.get("scfv", False),
        antibinder=ab,
        antigen=antigen,
    )


def test_the_flags_and_a_given_light_come_first(antibinder):
    assert _status(antibinder, light=L1) == ("paired_original", "source", "")
    assert _status(antibinder, nanobody=True)[0] == "vhh"
    assert _status(antibinder, scfv=True)[0] == "scfv"
    assert (
        _status(antibinder, heavy=H + "GGGGSGGGGS" + L1)[0] == "scfv"
    )  # a linker, no flag


def test_a_single_candidate_across_sources_is_filled(antibinder):
    # covid and biomap both give L1; hiv gives L2: three sources, two lights -> ambiguous
    assert _status(antibinder)[0] == "ambiguous"


def test_covid_alone(tmp_path):
    (tmp_path / "COVID-19").mkdir()
    (tmp_path / "HIV").mkdir()
    (tmp_path / "BioMap").mkdir()
    (tmp_path / "COVID-19/CoV-AbDab_only_sars2_filter.csv").write_text(
        f"Heavy,Light,Target,Protein_Epitope,Label,antigen\n{H},{L1},x,y,1,{ANT}\n"
    )
    (tmp_path / "HIV/dataset_hiv_cls.csv").write_text(
        "antibody,virus_seq,ANT_Binding,split\n"
    )
    (tmp_path / "BioMap/biomap.csv").write_text(
        "pdb,antibody_seq_a,antibody_seq_b,antigen_seq,delta_g\n"
    )
    ab = AntiBinderIndex(tmp_path)
    assert _status(ab) == ("filled_covid", "antibinder_covid", L1)
    assert (
        _status(ab, antigen="OTHER")[0] == "unresolved"
    )  # the antigen has to match too


def test_hiv_needs_a_light_sized_rest(tmp_path):
    for d in ("COVID-19", "HIV", "BioMap"):
        (tmp_path / d).mkdir()
    (tmp_path / "COVID-19/CoV-AbDab_only_sars2_filter.csv").write_text(
        "Heavy,Light,Target,Protein_Epitope,Label,antigen\n"
    )
    (tmp_path / "BioMap/biomap.csv").write_text(
        "pdb,antibody_seq_a,antibody_seq_b,antigen_seq,delta_g\n"
    )
    (tmp_path / "HIV/dataset_hiv_cls.csv").write_text(
        f"antibody,virus_seq,ANT_Binding,split\n{H + L2},{ANT},1,t\n{H + 'GGGG'},{ANT},1,t\n"
    )
    ab = AntiBinderIndex(tmp_path)
    assert _status(ab) == (
        "filled_hiv",
        "antibinder_hiv",
        L2,
    )  # the 4-aa rest is not a light chain
