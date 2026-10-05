import yaml

from structcooker.manifests import RELEASE, merged


def test_release_manifest_is_the_union_of_the_set_manifests():
    # Regenerate with: python -c 'from structcooker.manifests import write; write()'
    assert yaml.safe_load(RELEASE.read_text()) == merged()
