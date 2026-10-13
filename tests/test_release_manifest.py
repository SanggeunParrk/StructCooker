import yaml

from structcooker.manifests import RELEASE, merged


def test_release_manifest_is_the_union_of_the_set_manifests():
    # Regenerate with: python -c 'from structcooker.manifests import write; write()'
    assert yaml.safe_load(RELEASE.read_text()) == merged()


def test_release_exceptions_name_release_stages_with_reason_and_hash():
    import json

    from structcooker.manifests import DB_ROOT

    release = json.loads((DB_ROOT / "release-exceptions.json").read_text())
    assert len(release["pdb/cif"]) == 5
    stages = merged()
    for node, entries in release.items():
        assert node in stages
        assert all(e["reason"] and len(e["sha256"]) == 64 for e in entries.values())
