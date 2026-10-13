"""Recover the 14 OFD long structure.npz files that end in stale bytes (2026-09-14).

Each raw `preprocessed/<id>/structure.npz` holds a complete zip archive followed by
32-479 trailing bytes that contain a later end-of-directory marker, so zipfile picks the
wrong one ("Bad magic number for central directory"). The valid archive is the first
`retained_bytes` bytes. This writes those prefixes under intermediate (raw stays as
downloaded), checks both hashes against the 2026-09-14 recovery record, and writes
`overrides.json` for OPENFOLD_STRUCTURE_OVERRIDES (docs/history/long-cif-recovery-2026-09-14.md):

    {"<path under raw/openfold_distillation>": "<path relative to overrides.json>"}

Usage: python scripts/staging/recover_ofd_long_npz.py [--root /data/shared/cssb_data/BioMol_clean]
"""
import argparse
import os
import hashlib
import json
from pathlib import Path

# id, sha256 of the raw file, bytes of the valid archive, sha256 of that prefix
RECOVERED = (
    ("MGYP000485802788", "a67cc536c07347b0880fa8adf28c00cf6bfac769b452ee8fdc6905fce41c2a58", 32686, "b3cfcd9c5465133a36e581854d1667d81adc683bea0ae42a3ae117900797d597"),
    ("MGYP006508416081", "f0ed58e4af554dc0935c5ae01e2b08da6261b877694faaf6c99fbbe7069d4be6", 31046, "8c91cfa8b521d6895b210fcf1e26786f81fb76b10dd6e45f096ddc745ab4566a"),
    ("MGYP006787128943", "ab0685743ff5536bc471c6868c817d4e048ac95e306590ff30928180935944bf", 27091, "10cff37630dbc0b566a313ba3fc34eed98237cae8e6c179b4ab2844cfffa1272"),
    ("MGYP006867818277", "d8974fc20be03a36fc87081f8bd3ec22942dd6b5fc8bdc3445aebbf10b3944e8", 30104, "180fa34ab2291c55f320f394ea674d0497c33b3efe48093d13b0c7e130ff382f"),
    ("MGYP006893305931", "07896a160a2a7a8f45e310cb030fdbc11740b813543dda39563a50df580ae85c", 49605, "7aced8b5f393e66533c07fe3ffb029bc5720f7e1b6eea133c7cdb16439c29b13"),
    ("MGYP006909030808", "5ac0e1280de982483786b1448288ce98be39ba9f4a0de00f1c736b1d5ffd81d2", 63486, "b76c2b9ae216f7223516a0db4f84e048a651829b1b17470566ccfe1660f21544"),
    ("MGYP006915949634", "8d1eeb697c0266722d2e0de733218409859728eb2554abb69edf84b60988dc2b", 43407, "51fc5f512beb1e6098d158dabfd3629e1f41989b95c34339569cef4caf8b87e9"),
    ("MGYP006919632377", "3374cfdcebbcd6843f8898576f9f73cc158bd255e590b5951763c19e490e2b78", 40706, "7f72b890df55888781897d1d7da7aacf147db9b85d5278727ea235f00e735e54"),
    ("MGYP006923605014", "7a7d54a25ceda6214b1ddadad6558e9bf36459bb66a9fc75fad786e333fb0308", 46650, "d3a7de4c9fc2bc1110986701b8a4751a8fe3a90bfa4ea3d63ef0bf3c0382397f"),
    ("MGYP006928561759", "f662e2e19917e79bd9d18f043bbb6cc0cb3abb892a5a64f6e38f97b68bb79ae8", 53043, "d2e7373c62d7c9cd4d06fc6f362d2b597651b9a4894e156cfc6633eca836c999"),
    ("MGYP006935494712", "65f5c19f9a6a17dc903942b619dfcbe0028655ad8c3e859942c28e7b31081801", 54271, "1a6a74eb3d0240b70718de46d09238e855bd7a71aef110a9a3ed3b56dbd27f96"),
    ("MGYP006938809930", "0d8a8575a64ea68bb7738c993cf068793e8c2d5b7ac316a814c65a7153134ff7", 32996, "14de4faa260502d793fb2682aeda3c93c3892db3ef85e965d0a6b51feef23ad1"),
    ("MGYP006949394884", "ead426b9c3405c769da1d4ae42855b054ac6d11940e5321cb40fccf2bde23eea", 65329, "a3f78b869dbf4d646a08fd5bf2689462c97352ca646035570aac7e538a935bb9"),
    ("MGYP006957396418", "ab179de3c07a02abe8fe7dfee39f4ba5fc33329dd20531593e371373ea3d0692", 34302, "8494e775e01993e3f599307fd91360bd13d27f9da32a795e3c76583cb3a7ad6e"),)
SET = Path("monomer_distillation_sets_v2/long_monomers/preprocessed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(os.environ.get("BIOMOL_ROOT", "/data/shared/cssb_data/BioMol_clean")))
    args = parser.parse_args()
    raw = args.root / "materials/raw/openfold_distillation"
    out = args.root / "materials/intermediate/openfold_distillation/long_recovered"
    overrides = {}
    for key, original_sha, retained, recovered_sha in RECOVERED:
        data = (raw / SET / key / "structure.npz").read_bytes()
        if hashlib.sha256(data).hexdigest() != original_sha:
            raise SystemExit(f"{key}: raw structure.npz is not the file the recovery was made from")
        prefix = data[:retained]
        if hashlib.sha256(prefix).hexdigest() != recovered_sha:
            raise SystemExit(f"{key}: recovered prefix does not match the recorded hash")
        target = out / key / "structure.npz"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(prefix)
        overrides[str(SET / key / "structure.npz")] = str(target.relative_to(out))
    (out / "overrides.json").write_text(json.dumps(overrides, indent=1, sort_keys=True) + "\n")
    print(f"{len(overrides)} recovered structure.npz -> {out}")


if __name__ == "__main__":
    main()
