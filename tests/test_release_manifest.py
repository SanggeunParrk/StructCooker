import yaml

from structcooker.manifests import RELEASE, merged


def test_release_manifest_is_the_union_of_the_set_manifests():
    # Regenerate with: python -c 'from structcooker.manifests import write; write()'
    assert yaml.safe_load(RELEASE.read_text()) == merged()


def test_release_exceptions_cover_the_pdb_policy_and_name_release_stages():
    import json

    from structcooker.manifests import DB_ROOT

    release = json.loads((DB_ROOT / "release-exceptions.json").read_text())
    pdb = json.loads((DB_ROOT / "pdb-exceptions.json").read_text())
    for node, entries in pdb.items():
        assert entries.items() <= release[node].items()
    stages = merged()
    for node, entries in release.items():
        assert node in stages
        assert all(e["reason"] and len(e["sha256"]) == 64 for e in entries.values())
