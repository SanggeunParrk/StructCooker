"""Stage ASD v2.0 under materials/raw/asd/.

download/                       the four Drive zips and their extraction (extract/asd_2.0/), as downloaded
v2.0/                           the live Delta table: the parquet files the log's current version adds, plus _delta_log/
                                (hard links into download/extract, so nothing is copied)
SOURCE.tsv                      path, bytes, sha256 of every file under download/*.zip and v2.0/
README.md

Run through sbatch (hashing 5 GB): see stage_asd.sbatch.
"""

import argparse
import hashlib
import json
import os
import shutil

import pyarrow.parquet as pq


def live_files(delta_dir):
    log = f"{delta_dir}/_delta_log"
    ck = json.load(open(f"{log}/_last_checkpoint"))["version"]
    live = {}
    for r in pq.read_table(
        f"{log}/{ck:020d}.checkpoint.parquet", columns=["add"]
    ).to_pylist():
        if r["add"]:
            live[r["add"]["path"]] = r["add"]["size"]
    v = ck
    while os.path.exists(f"{log}/{v + 1:020d}.json"):
        v += 1
        for line in open(f"{log}/{v:020d}.json"):
            a = json.loads(line)
            if "add" in a:
                live[a["add"]["path"]] = a["add"]["size"]
            if "remove" in a:
                live.pop(a["remove"]["path"], None)
    return v, live


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="materials/raw/asd")
    ap.add_argument("--expect-rows", type=int, default=1177819)
    a = ap.parse_args()
    src = f"{a.root}/download/extract/asd_2.0"
    dst = f"{a.root}/v2.0"
    version, live = live_files(src)
    rows = 0
    for p, size in live.items():
        assert os.path.getsize(f"{src}/{p}") == size, p
        rows += pq.ParquetFile(f"{src}/{p}").metadata.num_rows
    assert rows == a.expect_rows, (rows, a.expect_rows)
    if os.path.exists(dst):
        shutil.rmtree(dst)
    os.makedirs(dst)
    for p in live:
        os.link(f"{src}/{p}", f"{dst}/{p}")
    shutil.copytree(f"{src}/_delta_log", f"{dst}/_delta_log", copy_function=os.link)
    with open(f"{a.root}/SOURCE.tsv", "w") as out:
        out.write("path\tbytes\tsha256\n")
        files = sorted(
            f"download/{z}"
            for z in os.listdir(f"{a.root}/download")
            if z.endswith(".zip")
        )
        for top, _, names in os.walk(dst):
            files += [os.path.relpath(f"{top}/{n}", a.root) for n in names]
        for rel in sorted(files):
            out.write(
                f"{rel}\t{os.path.getsize(f'{a.root}/{rel}')}\t{sha256(f'{a.root}/{rel}')}\n"
            )
    print(f"delta version {version}, {len(live)} live parquet, {rows} rows")


main()
