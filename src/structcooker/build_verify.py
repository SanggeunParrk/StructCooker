"""Full output checks required before writing a reusable build receipt."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from datacooker.config.runtime import activate_environment
from datacooker.lmdb import default_lmdb_key
from datacooker.lmdb.index import read_index
from datacooker.lmdb.sharded import open_env
from datacooker.utils.importing import resolve_object

from structcooker import preflight, schemas
from structcooker.build_state import Fingerprints
from structcooker.instructions.transforms.cifmol import (
    convert_to_cifmol_attached_transformed,
    convert_to_cifmol_dict,
    convert_to_cifmol_transformed,
)
from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.instructions.transforms.msa import msa_statistics
from structcooker.mols import CIFMol


def verify_structure(value: Any, schema: str) -> None:
    """Check every assembly's coordinate, parent-index and bond-index invariants."""
    adapters = {"A": convert_to_cifmol_dict, "B": convert_to_cifmol_attached_transformed,
                "I": convert_to_cifmol_transformed}
    if schema in adapters:
        assemblies = adapters[schema](value)
        mols = [item["cifmol"] for item in assemblies.values()]
    else:
        mols = [CIFMol.from_dict(value)]
    if not mols:
        msg = "Empty assembly map"
        raise ValueError(msg)
    for mol in mols:
        xyz = mol.atoms.xyz.value
        if xyz.ndim != 2 or xyz.shape[1] != 3 or not len(xyz) or not np.isfinite(xyz).any():
            msg = "Invalid coordinate array"
            raise ValueError(msg)
        chain = mol.index_table.atoms_to_chains(np.arange(len(xyz)))
        if chain.min() < 0 or chain.max() >= len(mol.chains.chain_id.value):
            msg = "Atom-to-chain index out of range"
            raise ValueError(msg)
        bonds = mol.atoms.bond_type
        for endpoints in (bonds.src_indices, bonds.dst_indices):
            if len(endpoints) and (endpoints.min() < 0 or endpoints.max() >= len(xyz)):
                msg = "Bond endpoint out of range"
                raise ValueError(msg)
        if schema == "B":
            for field in (mol.chains.seq_id.value, mol.chains.cluster_id.value):
                if len(field) != len(mol.chains.chain_id.value) or any(not str(v) for v in field):
                    msg = "Missing chain sequence/cluster identity"
                    raise ValueError(msg)


def verify_output(cfg: dict[str, Any], workdir: Path, allowed_failures: dict[str, dict[str, str]]) -> dict[str, Any]:
    """Check all records, complete index coverage, and exact build-input key coverage."""
    activate_environment(cfg)
    target = cfg.get("new_env_path") or cfg.get("env_path")
    if target is None:
        outputs = preflight.output_paths(cfg)
        if not outputs or any(not p.is_file() or not p.stat().st_size for p in outputs):
            msg = "Projection outputs are missing or empty"
            raise ValueError(msg)
        return {"outputs": list(map(str, outputs)), "verification": "nonempty files"}
    path = Path(target)
    indices = read_index(path)
    sizes = {item.key: item.value_bytes for item in indices}
    if len(sizes) != len(indices) or not sizes:
        msg = "Missing, empty, or duplicate-key index"
        raise ValueError(msg)
    metadata = json.loads(Path(str(path) + ".meta.json").read_text())
    schema = schemas.get(str(cfg["schema"]))
    if (metadata["schema"] != schema.name or metadata["entries"] != len(sizes)
            or metadata["total_bytes"] != sum(sizes.values())
            or metadata["max_bytes"] != max(sizes.values())):
        msg = "Index/metadata mismatch"
        raise ValueError(msg)
    checked = 0
    with open_env(str(path), readonly=True, lock=False, readahead=False) as env, env.begin(buffers=True) as txn:
        for raw_key, raw_value in txn.cursor():
            key = bytes(raw_key).decode()
            if sizes.get(key) != len(raw_value):
                msg = f"DB/index mismatch for {key}"
                raise ValueError(msg)
            value = bytes(raw_value) if schema.codec == schemas.Codec.raw else from_bytes(bytes(raw_value))
            issues = schema.validate(value)
            if issues:
                msg = f"Invalid {schema.name} record {key}: {issues}"
                raise ValueError(msg)
            if schema.name in {"A", "B", "C", "H", "I"}:
                if not isinstance(value, dict):
                    msg = f"Invalid structure payload for {key}"
                    raise ValueError(msg)
                verify_structure(value, schema.name)
            if schema.name == "E":
                if not isinstance(value, dict):
                    msg = f"Invalid MSA payload type for {key}"
                    raise ValueError(msg)
                msa = value["msa_dict"]
                sequence = msa["sequences"]
                aligned = sequence["aligned_sequences"]
                max_depth = cfg.get("parameters", {}).get("max_depth")
                if max_depth is not None and len(aligned) > int(max_depth):
                    msg = f"MSA depth exceeds configured cap for {key}"
                    raise ValueError(msg)
                profile = sequence["profile"]
                if (aligned.ndim != 2 or min(aligned.shape) == 0
                        or sequence["deletions"].shape != aligned.shape
                        or profile.shape[0] != aligned.shape[1]
                        or not np.isfinite(profile).all()
                        or not np.allclose(profile.sum(axis=1), 1)
                        or aligned.min() < 0 or aligned.max() >= profile.shape[1]
                        or sequence["deletions"].min() < 0 or sequence["deletions"].max() > 255
                        or any(len(v) != len(aligned) for v in msa["headers"].values())):
                    msg = f"Invalid MSA arrays for {key}"
                    raise ValueError(msg)
                deletion_mean, expected_profile = msa_statistics(aligned, sequence["deletions"], profile.shape[1])
                if (not np.array_equal(profile, expected_profile)
                        or not np.array_equal(sequence["deletion_mean"], deletion_mean)):
                    msg = f"Incorrect MSA statistics for {key}"
                    raise ValueError(msg)
            checked += 1
    if checked != len(sizes):
        msg = "Index contains keys absent from the DB"
        raise ValueError(msg)

    manifests = list(workdir.rglob("merge_shards.txt"))
    if len(manifests) != 1:
        msg = "Require exactly one native shard manifest for this attempt"
        raise ValueError(msg)
    manifest = manifests[0]
    failures: set[str] = set()
    totals = {"attempted": 0, "written": 0, "skipped_existing": 0, "skipped_empty": 0, "failed": 0}
    for shard in manifest.read_text().splitlines():
        report = json.loads(Path(shard + ".build-report.json").read_text())
        if report["failed"] != len(report["failed_keys"]):
            msg = f"Incomplete failure ledger for {shard}"
            raise ValueError(msg)
        failures.update(report["failed_keys"])
        for field in totals:
            totals[field] += report.get(field, 0)
    unexpected = failures - allowed_failures.keys()
    if unexpected:
        msg = f"Unapproved input failures: {sorted(unexpected)}"
        raise ValueError(msg)
    if (totals["written"] + totals["skipped_existing"] != checked or totals["attempted"] !=
            totals["written"] + totals["skipped_existing"] + totals["skipped_empty"] + totals["failed"]):
        msg = "Shard accounting does not match output records"
        raise ValueError(msg)
    if not cfg.get("old_env_path"):
        builder = cfg.get("key_builder")
        key_fn = resolve_object(builder) if isinstance(builder, str) else (builder if callable(builder) else default_lmdb_key)
        expected: set[str] = set()
        paths: dict[str, Path] = {}
        count = 0
        for filelist in manifest.parent.glob("items_*.txt"):
            for line in filelist.read_text().splitlines():
                key = str(key_fn(Path(line)))
                expected.add(key)
                paths[key] = Path(line)
                count += 1
        if count != len(expected):
            msg = "Duplicate source keys"
            raise ValueError(msg)
        if (not expected or not failures <= expected or count != totals["attempted"]
                or expected - failures != sizes.keys()):
            msg = "Input/output key coverage mismatch"
            raise ValueError(msg)
        fingerprints = Fingerprints()
        for key in failures:
            if fingerprints.file(paths[key]) != allowed_failures[key]["sha256"]:
                msg = f"Accepted failure input changed: {key}"
                raise ValueError(msg)
    elif failures:
        msg = "Failed rebuild records cannot be waived by a raw-input exception policy"
        raise ValueError(msg)
    if cfg.get("old_env_path") and schema.name == "E":
        source_keys = {item.key for item in read_index(Path(cfg["old_env_path"]))}
        if source_keys != sizes.keys():
            msg = "MSA rebuild must preserve every source key"
            raise ValueError(msg)
    return {"records": checked, "schema": schema.name, "build_counts": totals,
            "accepted_failures": {k: allowed_failures[k] for k in sorted(failures)}}
