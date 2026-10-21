import numpy as np

from structcooker.interface.detect import (
    ContactCfg,
    UnitCfg,
    dedup,
    find_units,
    inter_node_pairs,
    select_waters,
)
from structcooker.interface.ligands import is_binding_ligand, load_excluded


def _chain(x0, n=10, y=0.0, z=0.0):
    """N residues, one atom each, 3.8 A apart along x."""
    return np.array([[x0 + 3.8 * k, y, z] for k in range(n)])


def _build(chains, kinds):
    xyz = np.vstack(chains)
    node = np.concatenate([np.full(len(c), i) for i, c in enumerate(chains)])
    return xyz, node, np.arange(len(xyz)), kinds


def test_pairs_match_bruteforce():
    rng = np.random.default_rng(0)
    xyz = rng.uniform(0, 30, (300, 3))
    node = rng.integers(0, 5, 300)
    i, j, _ = inter_node_pairs(xyz, node, 5.0)
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    ii, jj = np.nonzero((d < 5.0) & (node[:, None] != node[None]) & (np.arange(300)[:, None] < np.arange(300)[None]))
    assert set(zip(i.tolist(), j.tolist(), strict=False)) == set(zip(ii.tolist(), jj.tolist(), strict=False))


def test_antibody_antigen_is_one_unit_of_three():
    # H-Ag and H-L share H residues 0..9 -> L joins the H-Ag unit
    H, L, A = _chain(0), _chain(0, y=4.0), _chain(0, y=-4.0)
    xyz, node, res, kinds = _build([H, L, A], ["protein"] * 3)
    assert find_units(xyz, node, res, kinds)[0].nodes == [0, 1, 2]


def test_far_chain_does_not_join():
    # B touches A at A's residues 0-9; C touches B only at B's far end: C must not join A-B via R
    A, B = _chain(0), _chain(0, y=4.0)
    C = _chain(40, y=4.0 + 4.0)  # contacts B only if x overlaps: here 40..76 vs B 0..34 -> no contact
    xyz, node, res, kinds = _build([A, B, C], ["protein"] * 3)
    units = find_units(xyz, node, res, kinds)
    assert [u.nodes for u in units] == [[0, 1]]


def test_one_hop_no_percolation():
    chains = [_chain(0, y=4.0 * k) for k in range(5)]  # a stack, each touches neighbours only
    xyz, node, res, kinds = _build(chains, ["protein"] * 5)
    for u in find_units(xyz, node, res, kinds):
        assert len(u.nodes) < 5


def test_cap_flags_truncated():
    chains = [_chain(0, y=4.0 * k) for k in range(2)] + [_chain(0, y=2.0, z=3.5 * k) for k in range(1, 6)]
    xyz, node, res, kinds = _build(chains, ["protein"] * 7)
    units = find_units(xyz, node, res, kinds, ucfg=UnitCfg(max_polymer_chains=3))
    assert all(len(u.nodes) <= 3 for u in units) and any(u.truncated for u in units)


def test_small_contact_rejected():
    a, b = _chain(0, n=10), _chain(0, n=3, y=4.0)
    xyz, node, res, kinds = _build([a, b], ["protein"] * 2)
    assert find_units(xyz, node, res, kinds, ContactCfg(min_residues=5)) == []


def test_ligand_needs_three_atoms():
    p = _chain(0)
    lig3 = np.array([[0, 4, 0], [3.8, 4, 0], [7.6, 4, 0.0]])
    xyz, node, res, kinds = _build([p, lig3], ["protein", "ligand"])
    assert len(find_units(xyz, node, res, kinds)) == 1
    xyz, node, res, kinds = _build([p, lig3[:2]], ["protein", "ligand"])
    assert find_units(xyz, node, res, kinds) == []


def test_ion_bridges_but_does_not_seed_alone():
    A = _chain(0, n=3)
    B = _chain(0, n=3, y=5.6)
    zn = np.array([[3.8, 2.8, 0.0]])                    # 2.8 A from both chains
    xyz, node, res, kinds = _build([A, B, zn], ["protein", "protein", "ion"])
    assert find_units(xyz, node, res, kinds)[0].nodes == [0, 1, 2]
    xyz, node, res, kinds = _build([A, np.array([[50.0, 50, 50]]), zn], ["protein", "protein", "ion"])
    assert find_units(xyz, node, res, kinds) == []      # ion touching one chain only: no unit
    far = np.array([[3.8, 3.3, 0.0]])                    # 3.3 A: outside the 3.0 A coordination shell
    xyz, node, res, kinds = _build([A, _chain(0, n=3, y=6.6), far], ["protein", "protein", "ion"])
    assert find_units(xyz, node, res, kinds) == []


def test_ion_only_counts_donor_atoms():
    A = _chain(0, n=3)
    B = _chain(0, n=3, y=5.6)
    zn = np.array([[3.8, 2.8, 0.0]])
    xyz, node, res, kinds = _build([A, B, zn], ["protein", "protein", "ion"])
    carbon_only = np.zeros(len(xyz), dtype=bool)
    assert find_units(xyz, node, res, kinds, donor=carbon_only) == []
    assert find_units(xyz, node, res, kinds, donor=~carbon_only)


def test_antibody_pair_is_forced_even_without_touching_the_epitope():
    # antigen A touches H only; L touches H but not A. L must still come along.
    H = _chain(0, y=0.0)
    L = _chain(0, y=-4.0)
    A = _chain(0, y=4.0)
    xyz, node, res, kinds = _build([H, L, A], ["protein"] * 3)
    # shorten L so it shares no residues with the H-A interface: move L away along x
    xyz[node == 1, 0] += 60.0
    # ... and give it a contact with H elsewhere (a tail of H next to L)
    tail = np.array([[60.0 + 3.8 * k, 0.0, 0.0] for k in range(10)])
    xyz = np.vstack([xyz, tail]); node = np.concatenate([node, np.zeros(10, int)]); res = np.arange(len(xyz))
    plain = find_units(xyz, node, res, kinds, ucfg=UnitCfg(join_min_residues=3))
    assert all(1 not in u.nodes or 2 not in u.nodes for u in plain)          # without the rule L is not bundled with A
    forced = find_units(xyz, node, res, kinds, ab_class=["H", "K", None])
    assert any(set(u.nodes) == {0, 1, 2} for u in forced)
    assert not any(set(u.nodes) == {0, 1} for u in forced)                    # the bare H-L unit is redundant now


def test_duplex_partner_is_forced():
    p = _chain(0, y=0.0)
    d1 = _chain(0, y=4.0)           # touches the protein
    d2 = _chain(0, y=40.0)          # far from the protein
    xyz, node, res, kinds = _build([p, d1, d2], ["protein", "na", "na"])
    assert all(2 not in u.nodes for u in find_units(xyz, node, res, kinds))
    units = find_units(xyz, node, res, kinds, duplex={(1, 2): 12})
    assert units[0].nodes == [0, 1, 2]


def test_waters_ignore_unit_unit_contacts():
    # two unit atoms 3 A apart (a unit-unit pair inside the cutoff) must not leak into the water result
    xyz = np.array([[0.0, 0, 0], [3.0, 0, 0]]); node = np.array([0, 1])
    keep, bridge = select_waters(xyz, node, np.array([[1.5, 3.0, 0.0], [50.0, 0, 0]]))
    assert keep.tolist() == [0] and bridge.tolist() == [True]


def test_waters_context_and_bridging():
    A, B = np.array([[0.0, 0, 0]]), np.array([[6.0, 0, 0]])
    xyz = np.vstack([A, B]); node = np.array([0, 1])
    w = np.array([[3.0, 0, 0], [0.0, 3.8, 0], [30.0, 0, 0]])
    keep, bridge = select_waters(xyz, node, w)
    assert keep.tolist() == [0, 1] and bridge.tolist() == [True, False]


def test_homodimer_symmetry_dedups():
    chains = [_chain(0), _chain(0, y=4.0), _chain(0, y=2.0, z=3.4)]
    xyz, node, res, kinds = _build(chains, ["protein"] * 3)
    units = find_units(xyz, node, res, kinds, ucfg=UnitCfg(max_polymer_chains=2))
    assert len(dedup(units, ["s1"] * 3)) == 1


def test_exclusion_list():
    ex = load_excluded(["af3_crystallization_aids", "solvent"])
    assert {"SO4", "GOL", "PEG", "HOH"} <= ex
    assert not is_binding_ligand("SO4", 20, ex, 6)
    assert not is_binding_ligand("ATP", 3, ex, 6)
    assert is_binding_ligand("ATP", 31, ex, 6)


def _helix_strands(n_bp):
    """Synthetic duplex: strand A N1/N6 (atoms 2.9 A from strand B's N3/O4), C1' 10.4 A apart."""
    names, xyz, res, strand = [], [], [], []
    for k in range(n_bp):
        x = 3.4 * k
        for sid, y_c1, y_n, code in ((0, 0.0, 3.5, "A"), (1, 10.4, 6.6, "T")):
            n_a, n_b = ("N1", "N6") if sid == 0 else ("N3", "O4")
            names += ["C1'", n_a, n_b]
            xyz += [[x, y_c1, 0.0], [x, y_n, 0.0], [x + 1.2, y_n, 0.0]]
            res += [2 * k + sid] * 3
            strand += [sid] * 3
    return np.array(names), np.array(xyz), np.array(res), np.array(strand)


def test_find_duplex_counts_base_pairs():
    from structcooker.interface.duplex import find_duplex

    names, xyz, res, strand = _helix_strands(8)
    code = np.where(strand == 0, "A", "T")
    got = find_duplex(xyz, strand, res, names, code, np.ones(len(xyz), bool), min_bp=4)
    assert got == {(0, 1): 8}
    names, xyz, res, strand = _helix_strands(3)
    code = np.where(strand == 0, "A", "T")
    assert find_duplex(xyz, strand, res, names, code, np.ones(len(xyz), bool), min_bp=4) == {}


def test_find_duplex_ignores_stacked_chains():
    from structcooker.interface.duplex import find_duplex

    names, xyz, res, strand = _helix_strands(8)
    xyz[strand == 1, 1] -= 6.5   # strand B 3.9 A from strand A in y: stacked, C1' too close
    code = np.where(strand == 0, "A", "T")
    assert find_duplex(xyz, strand, res, names, code, np.ones(len(xyz), bool), min_bp=4) == {}
