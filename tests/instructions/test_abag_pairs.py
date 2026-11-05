"""The pair table of ASD, AbRank and POISK: S2, antigen length, the binding verdicts of the three sources and the pair rule."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from structcooker.instructions.readers.abag import abag_header, load_abag_fasta
from structcooker.instructions.readers.abag_pairs import load_abag_pairs

H1, H2, H3 = "EVQLVESGGG" * 3, "QVQLQESGGG" * 3, "QVQLVQSGAE" * 3
L1 = "DIQMTQSPSS" * 3
AG1, AG2, AG3 = "MKTAYIAKQR" * 2, "GSHMLEDPVD" * 2, "AAAAKKKKLL" * 2
SHORT = "MKTAYIAKQ"  # 9 aa: under the minimum
SW_LIGHTISH, SW_HEAVYISH = (
    "EVQL" + "XXXXXXXXXFGQGTKVEIK",
    "QVQL" + "XXXXXXXXXWGQGTLVTVSS",
)  # AbRank gave them in the wrong columns

_MAP = pa.map_(pa.string(), pa.string())
_CHAINS = pa.struct(
    [
        ("antibody_chains", _MAP),
        ("target_chains", _MAP),
        ("antibody_chains_numbered", pa.map_(pa.string(), _MAP)),
    ]
)
_META = pa.struct([("target_name", pa.string()), ("target_pdb", pa.string())])


def _asd_row(dataset, kind, value, heavy, light, antigen, name="", nanobody=False):
    ab = [("heavy", heavy)] + ([("light", light)] if light else [])
    chains = {
        "antibody_chains": ab,
        "target_chains": [("full", antigen)],
        "antibody_chains_numbered": {},
    }
    return {
        "dataset": dataset,
        "affinity_type": kind,
        "affinity": value,
        "chains": chains,
        "metadata": {"target_name": name, "target_pdb": ""},
        "nanobody": nanobody,
        "scfv": False,
    }


@pytest.fixture
def sources(tmp_path):
    rows = [
        _asd_row("patents", "bool", "1", H1, L1, AG1, "Spike"),
        _asd_row("hiv", "bool", "0", H2, "", AG1),  # says no
        _asd_row("alphaseq", "alphaseq", "2.5", H3, "", AG2, nanobody=True),
        _asd_row("patents", "bool", "1", H3, "", SHORT),  # antigen too short
        _asd_row("patents", "bool", "1", "", L1, AG1),  # no heavy chain
    ]
    table = pa.table(
        {
            "dataset": [r["dataset"] for r in rows],
            "affinity_type": [r["affinity_type"] for r in rows],
            "affinity": [r["affinity"] for r in rows],
            "chains": pa.array([r["chains"] for r in rows], type=_CHAINS),
            "metadata": pa.array([r["metadata"] for r in rows], type=_META),
            "nanobody": [r["nanobody"] for r in rows],
            "scfv": [r["scfv"] for r in rows],
        }
    )
    asd = tmp_path / "asd"
    asd.mkdir()
    pq.write_table(table, asd / "part-00000-x-c000.snappy.parquet")
    header = "Ab_name\tAg_name\tIC50 [ug/mL]\tAffinity_Kd [nM]\tAb_heavy_chain_seq\tAb_light_chain_seq\tAg_seq\tSource\tAff_op\n"
    lines = [
        f"a\tSARS-CoV-2\t\t3.4\t{H1}\t{L1}\t{AG1}\tAlphaSeq\t=",  # the same triple as the ASD patents row: one pair
        f"b\tSARS-CoV-2\t\t5000\t{H2}\t{L1}\t{AG2}\tAlphaSeq\t=",  # too weak
        f"c\tX\t\t\t{H2}\t{L1}\t{AG3}\tRBD-escape\t=",  # escape only: no verdict
        f"d\tX\t\t2.0\t{H2}\t{L1}\t{AG1}:{AG2}\tAlphaSeq\t=",  # two antigen chains: left out
        f"e\tX\t\t2.0\t{SW_LIGHTISH}\t{SW_HEAVYISH}\t{AG3}\tCATNAP\t=",  # swapped columns
    ]
    (tmp_path / "abrank.tsv").write_text(header + "\n".join(lines) + "\n")
    poisk_header = (
        "antibody_type,vh.sequence_from_patent,vl.sequence_from_patent,target.sequence_from_patent,target.protein_name,"
        "target.binding_affinity.affinity_value,audit.verdict\n"
    )
    (tmp_path / "poisk.csv").write_text(
        poisk_header
        + f"Conventional,{H1},{L1},{AG3},PD-1,,CONFIRMED\nVHH,{H3},,{AG3},PD-1,Kd 3 nM,PLAUSIBLE\n"
        f"Conventional,{H2},{L1},{AG3},PD-1,,CONTRADICTED\n"
    )
    seq_ids = {
        "P1": AG1,
        "P2": AG2,
        "P3": AG3,
        "A1": H1,
        "A2": H2,
        "A3": H3,
        "A4": L1,
        "A5": SW_HEAVYISH,
        "A6": SW_LIGHTISH,
    }
    (tmp_path / "map.tsv").write_text(
        "".join(f"{k}\t{v}\n" for k, v in seq_ids.items())
    )
    return asd, tmp_path / "abrank.tsv", tmp_path / "poisk.csv", tmp_path / "map.tsv"


def test_pairs(sources):
    rows = load_abag_pairs(*sources)
    assert rows[0] == (
        "antigen_seq_id",
        "heavy_seq_id",
        "light_seq_id",
        "nanobody",
        "scfv",
        "datasets",
        "names",
    )
    pairs = {r[:3]: r for r in rows[1:]}
    # patents + AbRank AlphaSeq: one pair with both datasets; the ASD hiv row (no) is another triple and does not matter here
    assert pairs[("P1", "A1", "A4")][5] == "abrank:AlphaSeq,asd:patents"
    assert pairs[("P1", "A1", "A4")][6] == "SARS-CoV-2;Spike"
    assert ("P1", "A2", "A4") not in pairs  # the hiv row says no
    assert pairs[("P2", "A3", "")][3] == "1"  # the alphaseq row is a nanobody
    assert ("P2", "A2", "A4") not in pairs  # Kd 5000 nM
    assert (
        "P3",
        "A2",
        "A4",
    ) not in pairs  # escape only (and the POISK row is CONTRADICTED)
    assert (
        pairs[("P3", "A1", "A4")][5] == "poisk:Conventional"
    )  # a patent-asserted pair
    assert pairs[("P3", "A3", "")][5] == "poisk:VHH"  # Kd 3 nM
    assert (
        pairs[("P3", "A5", "A6")][5] == "abrank:CATNAP"
    )  # the swapped columns are put right: heavy and light
    assert all(len(r) == 7 for r in rows)


def test_the_antigen_has_to_be_long_enough_and_one_chain(sources):
    antigens = {r[0] for r in load_abag_pairs(*sources)[1:]}
    assert "P1" in antigens
    assert antigens <= {"P1", "P2", "P3"}


def test_fasta_holds_every_usable_chain_once(sources):
    asd, abrank, poisk, _ = sources
    antigens = load_abag_fasta(asd, abrank, poisk, "antigen")
    assert set(antigens.values()) >= {AG1, AG2, AG3}
    assert abag_header("T", AG1) in antigens
    antibody = load_abag_fasta(asd, abrank, poisk, "antibody")
    assert H1 in antibody.values()
    assert abag_header("L", L1) in antibody
