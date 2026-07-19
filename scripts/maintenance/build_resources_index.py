"""Build the resources.lmdb location index for the MiniWorld train_item loader.

The index is the single location-authority the training dataloader
(``ResourceLocator``) consults to turn a record's ``source`` (+ record_id / seq_id
/ template key) into the exact LMDB path for cif / msa / template. Storage is
minimal and base-relative so the index is portable across mounts:

  * cif / template / single-shard msa are ONE lmdb per source -> recorded in a
    small ``__meta__`` json (source -> base-relative paths). No per-key rows.
  * MULTI-shard msa (distillation_long's 80 seqid shards) is the only thing that
    needs a per-key index: the shards were rekeyed to seq_id from record-id
    ranges, so a seq_id is not computable to a shard. Those get
    ``{source}{sep}{seq_id} -> shard_idx`` rows (2-byte value).

Paths are read from a data-config yaml's ``train_db`` block (OmegaConf only, no
MiniWorld import), so the index always matches what training consumes.

    pixi run python -m scripts.maintenance.build_resources_index \
        --data-config /home/psk6950/MiniWorld/configs/miniworld/data/\
full_train_af3_mix_no_disordered_bioai.yaml \
        --output /data/psk6950/filtering/edge_node/resources.lmdb
"""

from __future__ import annotations

import argparse
import json
import struct
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import lmdb
from omegaconf import OmegaConf

from structcooker.utils.seq_id import seq_id_from_name

KEY_SEP = "|"
META_KEY = b"__meta__"
DEFAULT_BASE = "/data/shared/cssb_data/BioMol/"
_COMMIT_EVERY = 1_000_000


def _rel(path: str | None, base: str) -> str | None:
    """Store paths relative to ``base`` (a resolver-time default the config can
    override). Paths outside base stay absolute; the resolver keeps those as-is."""
    if path is None:
        return None
    text = str(path)
    if text.startswith(base):
        return text[len(base):].lstrip("/")
    return text


def _shard_keys(path: str) -> list[str]:
    """Return all keys of one MSA shard LMDB (keys only, no values)."""
    env = lmdb.open(
        path, readonly=True, lock=False, readahead=False, max_readers=2048,
    )
    with env.begin() as txn:
        keys = [k.decode() for k in txn.cursor().iternext(keys=True, values=False)]
    env.close()
    return keys


def _sources_from_config(db: object, base: str) -> dict[str, dict]:
    """Build the per-source base-relative path table from a train_db block."""
    sources: dict[str, dict] = {}
    if db.get("pdb") is not None:
        pdb = db.pdb
        sources["pdb"] = {
            "cif": _rel(pdb.cif_db_path, base),
            "template": _rel(pdb.get("template_db_path"), base),
            "msa": [_rel(pdb.a3m_db_path, base)] if pdb.get("a3m_db_path") else [],
            "msa_rna": _rel(pdb.get("a3m_rna_db_path"), base),
            "sharded_msa": False,
        }
    for src in db.get("distillation_sources", []) or []:
        msa = [_rel(p, base) for p in (src.get("a3m_db_paths") or [])]
        sources[src.name] = {
            "cif": _rel(src.cif_db_path, base),
            "template": _rel(src.get("template_db_path"), base),
            "msa": msa,
            "msa_rna": None,
            "sharded_msa": len(msa) > 1,
        }
    return sources


def build(data_config: Path, output: Path, workers: int, base: str) -> int:
    """Write meta + the seq_id -> shard index for every multi-shard msa source."""
    if not base.endswith("/"):
        base += "/"
    db = OmegaConf.load(data_config).train_db
    sources = _sources_from_config(db, base)
    meta = {"version": 1, "key_sep": KEY_SEP, "base": base, "sources": sources}

    output.parent.mkdir(parents=True, exist_ok=True)
    env = lmdb.open(str(output), map_size=8 * 1024**3, subdir=True)
    with env.begin(write=True) as txn:
        txn.put(META_KEY, json.dumps(meta).encode())

    total = 0
    for name, info in sources.items():
        if not info["sharded_msa"]:
            continue
        shards = [p if p.startswith("/") else base + p for p in info["msa"]]
        print(f"[{name}] {len(shards)} msa shards -> seq_id index", flush=True)
        n_src = skipped = 0
        pending = 0
        txn = env.begin(write=True)
        # Read shards in parallel (IO-bound); write seq_id->shard_idx serially.
        # First-shard-wins dedup (overwrite=False): identical seq_ids share an MSA.
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for shard_idx, keys in enumerate(pool.map(_shard_keys, shards)):
                val = struct.pack("<H", shard_idx)
                for seq_id in keys:
                    if seq_id_from_name(seq_id) is None:
                        skipped += 1
                        continue
                    txn.put(f"{name}{KEY_SEP}{seq_id}".encode(), val, overwrite=False)
                    pending += 1
                    if pending >= _COMMIT_EVERY:
                        txn.commit()
                        txn = env.begin(write=True)
                        pending = 0
                n_src += len(keys)
                print(
                    f"  [{name}] shard {shard_idx}: {len(keys):,} keys "
                    f"(src total {n_src:,})",
                    flush=True,
                )
        txn.commit()
        total += n_src
        note = f" ({skipped:,} malformed keys skipped)" if skipped else ""
        print(f"[{name}] done: {n_src:,} seq_id rows{note}", flush=True)

    env.sync()
    env.close()
    print(f"DONE: {total:,} seq_id rows + meta -> {output}", flush=True)
    return 0


def main() -> int:
    """CLI entry: build resources.lmdb from a data-config yaml."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-config", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument(
        "--base", default=DEFAULT_BASE,
        help="Default resolver root; paths stored relative to it.",
    )
    args = ap.parse_args()
    return build(args.data_config, args.output, args.workers, args.base)


if __name__ == "__main__":
    raise SystemExit(main())
