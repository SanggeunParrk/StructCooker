import numpy as np

from structcooker.instructions.transforms.geometry import chain_contacts_grid


def test_unique_contact_pairs_match_dense_reference_and_keep_legacy_mode():
    rng = np.random.default_rng(17)
    xyz = rng.uniform(-12, 12, size=(90, 3))
    xyz[0] = np.nan
    xyz[1] = [0, 0, 0]
    xyz[2] = [6, 0, 0]
    chain = np.arange(len(xyz)) % 3
    expected = {}
    for i in range(len(xyz)):
        for j in range(i + 1, len(xyz)):
            if chain[i] != chain[j] and np.sum((xyz[i] - xyz[j]) ** 2) <= 36:
                pair = tuple(sorted((int(chain[i]), int(chain[j]))))
                expected[pair] = expected.get(pair, 0) + 1
    def edges(once):
        src, dst, counts = chain_contacts_grid(xyz, chain, 6, count_atom_pairs_once=once)
        return {(int(a), int(b)): int(v) for a, b, v in zip(src, dst, counts, strict=True)}
    assert edges(True) == expected
    assert edges(False) == {key: 2 * value for key, value in expected.items()}
