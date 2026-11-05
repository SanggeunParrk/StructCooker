"""Key contracts: can a consumer find a chain's template from the chain's structure record?

A template DB only helps if the key a consumer derives from a structure record is the key the DB was written
under. Every DB below is keyed differently (by chain, by sequence, by entry), and a mismatch is silent: the lookup
finds nothing and the model trains without templates. Disordered templates were keyed by label letter (``1a0c_A``)
while the structures number their chains (``1a0c_1``), and for a month no chain matched. ``release-check`` now
samples structure records, derives each protein chain's key by the rule below and requires the key to be in the
template DB.

Rules (``chain_id`` and ``seq_id`` are the structure's ``chains.nodes`` columns)::

    label   {record}_{chain_id minus its copy suffix}    PDB, disordered   (``A_1`` -> ``101m_A``, ``1`` -> ``1a0c_1``)
    record  {record}                                      short / long OpenFold monomers (one chain per record)
    seq_id  {seq_id}                                      Teddymer, AFDB multimer (templates are per sequence)

A new template DB or a new key rule is added here, with the structure DB it is looked up from.
"""
from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np

SAMPLE_RECORDS = 300


@dataclass(frozen=True)
class Contract:
    name: str
    structure: str  # attached cif DB, relative to ``<BIOMOL_ROOT>/lmdb``, without ``.lmdb``
    template: str
    rule: str
    min_hit: float = 0.9  # a chain may lack a template key (e.g. not a polypeptide the template build covers)


CONTRACTS = (
    Contract("pdb", "pdb/cif_attached/cif_pdb_attached", "chain/template/template", "label"),
    Contract("disordered", "distillation_disordered/cif_attached/cif_disordered_attached",
             "distillation_disordered/template/template_disordered", "label"),
    Contract("short", "distillation_short/cif_attached/cif_short_attached",
             "distillation_short/template/template_topn_short", "record"),
    Contract("long", "distillation_long/cif_attached/cif_long_attached",
             "distillation_long/template/template_topn_long", "record"),
    Contract("teddymer", "teddymer/cif_attached/cif_teddymer_attached",
             "teddymer/template/template_teddymer", "seq_id"),
    Contract("afm_homodimer", "afdb_multimer/cif_attached/cif_afm_homodimer_attached",
             "afdb_multimer/template/template_afm", "seq_id"),
    Contract("afm_heterodimer", "afdb_multimer/cif_attached/cif_afm_heterodimer_attached",
             "afdb_multimer/template/template_afm", "seq_id"),
)


def _label(record: str, chain_id: str, seq_id: str) -> str:
    return f"{record}_{chain_id.rsplit('_', 1)[0]}"


def _record(record: str, chain_id: str, seq_id: str) -> str:
    return record


def _seq_id(record: str, chain_id: str, seq_id: str) -> str:
    return seq_id


RULES: dict[str, Callable[[str, str, str], str]] = {"label": _label, "record": _record, "seq_id": _seq_id}


def _column(nodes: dict[str, Any], name: str) -> list[str]:
    col = nodes[name]
    return [str(x) for x in np.asarray(col["value"] if isinstance(col, dict) else col).ravel()]


def record_keys(record: str, value: dict, rule: str) -> list[str]:
    """Template keys of the polypeptide chains of one decoded structure record."""
    derive = RULES[rule]
    keys = []
    for assembly in value.values():
        nodes = assembly["cifmol_attached_dict"]["chains"]["nodes"]
        for cid, etype, sid in zip(_column(nodes, "chain_id"), _column(nodes, "entity_type"), _column(nodes, "seq_id")):
            if "polypeptide" in etype:
                keys.append(derive(record, cid, sid))
    return keys


def check(contract: Contract, lmdb_root: Path, records: int = SAMPLE_RECORDS) -> dict[str, Any]:
    """Sample structure records of ``contract`` and count the chains whose derived key is in the template DB."""
    from datacooker.lmdb.sharded import open_env

    from structcooker.instructions.transforms.codecs import from_bytes

    report: dict[str, Any] = {"contract": contract.name, "rule": contract.rule, "problems": []}
    struct, templ = lmdb_root / f"{contract.structure}.lmdb", lmdb_root / f"{contract.template}.lmdb"
    for path in (struct, templ):
        if not path.exists():
            report["problems"].append(f"missing {path}")
            return report
    index = Path(f"{struct}.index.tsv")
    with index.open() as handle:
        names = [ln.split("\t", 1)[0] for ln in handle if not ln.startswith(("#", "key\t"))]
    chosen = random.Random(contract.name).sample(names, min(records, len(names)))
    chains = hit = 0
    missing: list[str] = []
    with open_env(str(struct), readonly=True, lock=False, readahead=False) as senv, \
            open_env(str(templ), readonly=True, lock=False, readahead=False) as tenv, \
            senv.begin() as stxn, tenv.begin() as ttxn:
        for name in chosen:
            raw = stxn.get(name.encode())
            if raw is None:
                report["problems"].append(f"{name} is in the index but not in the DB")
                continue
            for key in record_keys(name, from_bytes(bytes(raw)), contract.rule):
                chains += 1
                if ttxn.get(key.encode()) is not None:
                    hit += 1
                elif len(missing) < 5:
                    missing.append(key)
    report.update(records=len(chosen), chains=chains, hit=hit, missing_examples=missing)
    rate = hit / chains if chains else 0.0
    report["hit_rate"] = round(rate, 4)
    if chains == 0:
        report["problems"].append("no polypeptide chain sampled")
    elif rate < contract.min_hit:
        report["problems"].append(
            f"only {rate:.1%} of {chains} protein chains find their template key (rule {contract.rule}; "
            f"need {contract.min_hit:.0%}); e.g. {missing[:3]}")
    return report


def check_all(lmdb_root: Path, records: int = SAMPLE_RECORDS) -> list[dict[str, Any]]:
    return [check(c, lmdb_root, records) for c in CONTRACTS]
