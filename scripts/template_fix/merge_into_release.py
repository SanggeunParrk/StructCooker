"""Append the template_fix delta (chains / seq_ids that the 2026-07 template build never processed) to the release.

    python merge_into_release.py DELTA_ROOT RELEASE_ROOT [--apply]

A chain key already in the release is left alone unless its record holds no template and the new one does (``--replace-empty``). A seq_id already in the release (its older chains were built, the new ones
were not) gets the union of the old and the new template mols, so Phase 4 can be re-run from the release's seq_id DB. Chain keys are written the way the release has them
(pdb id lower-case, chain id as is), values are copied byte for byte. The ``.index.tsv`` and ``.meta.json`` next to each DB
are rewritten from the merged sizes; the old pair is kept as ``*.pre_template_fix``.
"""

from __future__ import annotations

import argparse
import json
import shutil
import statistics
import sys
from pathlib import Path

import lmdb

DBS = ("template", "seqid_template_mols")


def human(n: int) -> str:
    x = float(n)
    for u in ("B", "K", "M", "G", "T"):
        if x < 1024 or u == "T":
            return f"{x:.1f}{u}" if u != "B" else f"{int(x)}B"
        x /= 1024
    return str(n)


def chain_key(k: bytes) -> bytes:
    pdb, _, chain = k.decode().partition("_")
    return f"{pdb.lower()}_{chain}".encode()


def read_index(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open() as f:
        for line in f:
            if line.startswith(("#", "key\t")):
                continue
            k, n, *_ = line.rstrip("\n").split("\t")
            sizes[k] = int(n)
    return sizes


def write_index(path: Path, name: str, sizes: dict[str, int]) -> None:
    order = sorted(sizes.items(), key=lambda kv: -kv[1])
    vals = [n for _, n in order]
    total = sum(vals)
    with path.open("w") as f:
        f.write(
            f"# db\t{name}\n# entries\t{len(vals)}\n# total_bytes\t{total}\t{human(total)}\n"
        )
        f.write(f"# max_bytes\t{vals[0]}\t{human(vals[0])}\t{order[0][0]}\n")
        med = int(statistics.median(vals))
        f.write(
            f"# median_bytes\t{med}\t{human(med)}\n# n_value_gt_100MB\t{sum(v > 100 * 2**20 for v in vals)}\n"
        )
        f.write(
            f"# n_value_gt_1GB\t{sum(v > 2**30 for v in vals)}\nkey\tvalue_bytes\tvalue_human\n"
        )
        for k, n in order:
            f.write(f"{k}\t{n}\t{human(n)}\n")


def union_value(old: bytes, new: bytes) -> bytes:
    from structcooker.instructions.transforms.codecs import from_bytes, to_bytes

    a, b = from_bytes(old), from_bytes(new)
    mols = {**b.get("template_mols", {}), **a.get("template_mols", {})}
    return to_bytes({"template_mols": mols, "template_ids": list(mols)})


def is_empty(value: bytes) -> bool:
    from structcooker.instructions.transforms.codecs import from_bytes

    return not from_bytes(value).get("template_mols")


def merge(
    delta: Path,
    rel: Path,
    name: str,
    keyfn,
    apply: bool,
    union: bool = False,
    replace_empty: bool = False,
) -> dict:
    src = lmdb.open(
        str(delta / f"{name}.lmdb"), readonly=True, lock=False, readahead=False
    )
    dst = lmdb.open(
        str(rel / f"{name}.lmdb"),
        map_size=int(300e9),
        readahead=False,
        metasync=True,
        sync=True,
    )
    sizes = read_index(rel / f"{name}.lmdb.index.tsv")
    added = skipped = merged = 0
    batch: list[tuple[bytes, bytes]] = []

    def flush() -> None:
        with dst.begin(write=True) as t:
            for k, v in batch:
                t.put(k, v)
        batch.clear()

    with src.begin() as t, dst.begin() as r:
        for k, v in t.cursor():
            nk = keyfn(k)
            old = r.get(nk)
            if old is not None and union:
                nv = union_value(bytes(old), bytes(v))
                if nv != bytes(old):
                    merged += 1
                    sizes[nk.decode()] = len(nv)
                    if apply:
                        batch.append((nk, nv))
                        if len(batch) >= 200:
                            flush()
                else:
                    skipped += 1
                continue
            if (
                old is not None
                and replace_empty
                and is_empty(bytes(old))
                and not is_empty(bytes(v))
            ):
                merged += 1
                sizes[nk.decode()] = len(v)
                if apply:
                    batch.append((nk, bytes(v)))
                    if len(batch) >= 2000:
                        flush()
                continue
            if old is not None or nk.decode() in sizes:
                skipped += 1
                continue
            sizes[nk.decode()] = len(v)
            added += 1
            if apply:
                batch.append((nk, v))
                if len(batch) >= 2000:
                    flush()
    if apply and batch:
        flush()
    if apply:
        dst.sync(True)
        for suffix in (".index.tsv", ".meta.json"):
            p = rel / f"{name}.lmdb{suffix}"
            shutil.copy2(p, p.with_name(p.name + ".pre_template_fix"))
        write_index(rel / f"{name}.lmdb.index.tsv", f"{name}.lmdb", sizes)
        mp = rel / f"{name}.lmdb.meta.json"
        meta = json.loads(mp.read_text())
        vals = list(sizes.values())
        meta.update(entries=len(vals), total_bytes=sum(vals), max_bytes=max(vals))
        mp.write_text(json.dumps(meta, indent=2))
    return {
        "db": name,
        "added": added,
        "unioned": merged,
        "already_present": skipped,
        "entries_after": len(sizes),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    (
        ap.add_argument("delta"),
        ap.add_argument("release"),
        ap.add_argument("--apply", action="store_true"),
    )
    ap.add_argument(
        "--replace-empty",
        action="store_true",
        help="replace release chain records that hold no template",
    )
    a = ap.parse_args()
    d, r = (
        Path(a.delta) / "lmdb/chain/template",
        Path(a.release) / "lmdb/chain/template",
    )
    out = [
        merge(d, r, "seqid_template_mols", lambda k: k, a.apply, union=True),
        merge(d, r, "template", chain_key, a.apply, replace_empty=a.replace_empty),
    ]
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
