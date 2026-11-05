"""Key contracts: the template key a consumer derives from a structure record, and the disordered key translation."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from structcooker.instructions.transforms import filelists
from structcooker.release import key_contracts


def _record(chain_ids, entity_types, seq_ids):
    nodes = {
        "chain_id": {"value": np.array(chain_ids)},
        "entity_type": {"value": np.array(entity_types)},
        "seq_id": {"value": np.array(seq_ids)},
    }
    return {"1_1_.": {"cifmol_attached_dict": {"chains": {"nodes": nodes}}}}


def test_rules_use_polypeptide_chains_only():
    rec = _record(
        ["A_1", "B_1", "C_1"],
        ["polypeptide(L)", "polypeptide(L)", "polydeoxyribonucleotide"],
        ["P1", "P2", "D3"],
    )
    assert key_contracts.record_keys("101m", rec, "label") == ["101m_A", "101m_B"]
    assert key_contracts.record_keys("MGYP1", rec, "record") == ["MGYP1", "MGYP1"]
    assert key_contracts.record_keys("x", rec, "seq_id") == ["P1", "P2"]


def test_label_rule_on_numbered_chains():
    rec = _record(["1", "2", "10"], ["polypeptide(L)"] * 3, ["P1", "P1", "P2"])
    assert key_contracts.record_keys("1a0c", rec, "label") == [
        "1a0c_1",
        "1a0c_2",
        "1a0c_10",
    ]


def test_every_contract_has_a_rule():
    assert {c.rule for c in key_contracts.CONTRACTS} <= set(key_contracts.RULES)


def _disordered_entry(root: Path, entry: str, structure, templates):
    """structure: [(chain_id, label, entity)], templates: {label: entity}."""
    (root / "structure_files" / entry).mkdir(parents=True)
    cid, lab, ent = (np.array(x) for x in zip(*structure, strict=False))
    np.savez(
        root / "structure_files" / entry / f"{entry}.npz",
        chain_id=cid,
        label_asym_id=lab,
        entity_id=ent,
    )
    (root / "templates" / entry).mkdir(parents=True)
    for label, entity in templates.items():
        np.savez(
            root / "templates" / entry / f"{entry}_{label}.npz",
            entity_id=np.array([entity]),
        )


def _links(root: Path):
    return filelists.disordered_template_links(
        root / "structure_files",
        root / "templates",
        root / "links",
        root / "list.txt",
        workers=2,
    )


def _listed(root: Path):
    return sorted((root / "list.txt").read_text().split())


def test_disordered_links_by_label_with_copies(tmp_path):
    _disordered_entry(
        tmp_path,
        "1a0e",
        [("1", "A", "1"), ("11", "A", "1"), ("2", "B", "1"), ("3", "C", "2")],
        {"A": "1", "B": "1"},
    )
    _links(tmp_path)
    assert _listed(tmp_path) == [
        "1a0e/1a0e_1.npz",
        "1a0e/1a0e_11.npz",
        "1a0e/1a0e_2.npz",
    ]
    assert (tmp_path / "links/1a0e/1a0e_11.npz").resolve() == (
        tmp_path / "templates/1a0e/1a0e_A.npz"
    ).resolve()
    assert not (tmp_path / "links/1a0e/1a0e_11.npz").readlink().is_absolute()


def test_disordered_links_fall_back_to_the_first_file_of_the_entity(tmp_path):
    # 1a0o: the structure keeps one copy per entity (labels A, B); the templates are the asymmetric unit A..D
    _disordered_entry(
        tmp_path,
        "1a0o",
        [("1", "A", "1"), ("2", "B", "2"), ("3", "I", "3")],
        {"C": "1", "D": "2", "E": "1", "F": "2"},
    )
    msg = _links(tmp_path)
    assert _listed(tmp_path) == ["1a0o/1a0o_1.npz", "1a0o/1a0o_2.npz"]
    assert (tmp_path / "links/1a0o/1a0o_1.npz").resolve() == (
        tmp_path / "templates/1a0o/1a0o_C.npz"
    ).resolve()
    assert (tmp_path / "links/1a0o/1a0o_2.npz").resolve() == (
        tmp_path / "templates/1a0o/1a0o_D.npz"
    ).resolve()
    assert "2 by entity" in msg
    assert "2 template files unused" in msg


def test_disordered_links_chain_without_a_template_gets_no_key(tmp_path):
    _disordered_entry(tmp_path, "1abc", [("1", "A", "1"), ("2", "B", "2")], {"A": "1"})
    _links(tmp_path)
    assert _listed(tmp_path) == ["1abc/1abc_1.npz"]


def test_rna_msa_t_is_u():
    from structcooker.instructions.transforms.msa import parse_sequence

    out = parse_sequence(["GAUT", "GA-T", "GAUU"], a3m_type="rna")
    assert "".join(out["query_sequence"]) == "GAUU"
    assert out["aligned_sequences"][0].tolist() == [
        23,
        21,
        22,
        22,
    ]  # G A U U, no "other" token
    assert out["aligned_sequences"][1].tolist() == [23, 21, 31, 22]
    protein = parse_sequence(["MKT"], a3m_type="protein")
    assert (
        "".join(protein["query_sequence"]) == "MKT"
    )  # proteins keep their T (threonine)


def test_rna_openfold_msa_features_t_is_u():
    import numpy as np

    from structcooker.instructions.transforms.openfold import build_msa_features

    msa = np.array([list("GAUT"), list("GAUU")])
    out = build_msa_features(msa, np.zeros(msa.shape, dtype=np.int32), a3m_type="rna")
    assert "".join(out["query_sequence"]) == "GAUU"
    assert out["aligned_sequences"][0].tolist() == [23, 21, 22, 22]
    assert out["aligned_sequences"][1].tolist() == [23, 21, 22, 22]


def test_msa_alphabet_follows_the_search_sources():
    from structcooker.instructions.transforms.openfold import resolve_a3m_type

    assert (
        resolve_a3m_type({"rnacentral_hits": {}, "rfam_hits": {}}, "protein") == "rna"
    )
    assert resolve_a3m_type({"nucleotide_collection_hits": {}}, "protein") == "rna"
    assert resolve_a3m_type({"mmseqs_colabfold": {}}, "protein") == "protein"
    assert resolve_a3m_type({"mmseqs_colabfold": {}}, None) == "protein"
