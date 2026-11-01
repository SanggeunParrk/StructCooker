"""Anti-leakage check of stored PDB templates: every template of a chain must be released on or before the cutoff (2021-09-30),
at least 60 days before the query chain's entry, and not be the query's own entry.

    python template_date_audit.py DB [--sample N]      DB: a chain-keyed template.lmdb (keys <pdbid>_<chain>)
"""
from __future__ import annotations

import argparse
import datetime as dt
import random

from datacooker.lmdb.sharded import open_env

from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.instructions.transforms.template_candidates import _load_pdb_dates

META = "/data/shared/cssb_data/BioMol_clean/metadata/cif_metadata.tsv"
CUTOFF, GAP = dt.date(2021, 9, 30), 60


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("db")
    ap.add_argument("--sample", type=int, default=0)
    a = ap.parse_args()
    dates = _load_pdb_dates(META)
    with open_env(a.db, readonly=True, lock=False) as env, env.begin() as txn:
        keys = [bytes(k) for k, _ in txn.cursor()]
        if a.sample:
            keys = random.Random(0).sample(keys, min(a.sample, len(keys)))
        chains = mols = late = close = self_ = nodate = 0
        ex = []
        for k in keys:
            chain = k.decode()
            q = dates.get(chain.split("_")[0].lower())
            ids = list(from_bytes(bytes(txn.get(k)))["template_mols"])
            chains += 1
            for t in ids:
                mols += 1
                td = dates.get(t.split("_")[0].lower())
                if q is None or td is None:
                    nodate += 1
                elif t.split("_")[0].lower() == chain.split("_")[0].lower():
                    self_ += 1
                    ex.append((chain, t, "same entry"))
                elif td > CUTOFF:
                    late += 1
                    ex.append((chain, t, f"released {td} > cutoff"))
                elif (q - td).days < GAP:
                    close += 1
                    ex.append((chain, t, f"only {(q - td).days} days before the query"))
    print(f"{a.db}: chains={chains} template_mols={mols} after_cutoff={late} within_{GAP}d={close} same_entry={self_} "
          f"no_date={nodate} examples={ex[:5]}")


if __name__ == "__main__":
    main()
