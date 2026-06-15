"""Post-build QC for the ingest LMDBs (msa / cif / template).

Scans a built LMDB, decodes every record with the StructCooker codec and runs
the matching per-record validator, then reports how many records are valid and
a breakdown of issue codes. Use after a build to confirm the DB conforms to the
reference BioMolDB schema.

    pixi run python -m scripts.maintenance.validate_lmdb --type cif DB_PATH --limit 100
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import lmdb

from structcooker.instructions.transforms.codecs import from_bytes
from structcooker.instructions.transforms.validate import (
    validate_cif_record,
    validate_msa_record,
    validate_template_record,
)

VALIDATORS = {
    "msa": validate_msa_record,
    "cif": validate_cif_record,
    "template": validate_template_record,
}


def validate_lmdb(env_path: Path, db_type: str, limit: int | None) -> int:
    """Validate records in an LMDB and print a summary; return the invalid count."""
    validator = VALIDATORS[db_type]
    env = lmdb.open(str(env_path), readonly=True, lock=False, subdir=True, max_dbs=0)
    total = valid = invalid = decode_errors = 0
    issue_codes: Counter[str] = Counter()
    sample_failures: list[tuple[str, list[dict]]] = []
    with env.begin() as txn:
        for key, value in txn.cursor():
            if limit is not None and total >= limit:
                break
            total += 1
            try:
                record = from_bytes(value)
            except Exception as exc:  # noqa: BLE001 - report, do not abort the scan
                decode_errors += 1
                issue_codes["decode_error"] += 1
                if len(sample_failures) < 10:  # noqa: PLR2004 - keep a short sample
                    sample_failures.append((key.decode(errors="replace"), [{"code": "decode_error", "message": str(exc)}]))
                continue
            issues = validator(record)
            if issues:
                invalid += 1
                for issue in issues:
                    issue_codes[issue["code"]] += 1
                if len(sample_failures) < 10:  # noqa: PLR2004 - keep a short sample
                    sample_failures.append((key.decode(errors="replace"), issues))
            else:
                valid += 1
    env.close()

    out = sys.stdout.write
    out(f"DB:       {env_path}\n")
    out(f"type:     {db_type}\n")
    out(f"scanned:  {total}\n")
    out(f"valid:    {valid}\n")
    out(f"invalid:  {invalid} (decode errors: {decode_errors})\n")
    if issue_codes:
        out("issue codes:\n")
        for code, count in issue_codes.most_common():
            out(f"  {code}: {count}\n")
    if sample_failures:
        out("sample failures:\n")
        for key, issues in sample_failures:
            shown = "; ".join(f"{i['code']}:{i['message']}" for i in issues[:5])
            out(f"  [{key}] {shown}\n")
    out("RESULT: " + ("PASS" if invalid == 0 and decode_errors == 0 else "FAIL") + "\n")
    return invalid + decode_errors


def main() -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Validate a built ingest LMDB.")
    parser.add_argument("env_path", type=Path, help="Path to the LMDB directory.")
    parser.add_argument("--type", required=True, choices=sorted(VALIDATORS), dest="db_type")
    parser.add_argument("--limit", type=int, default=None, help="Validate at most N records.")
    args = parser.parse_args()
    return 1 if validate_lmdb(args.env_path, args.db_type, args.limit) else 0


if __name__ == "__main__":
    raise SystemExit(main())
