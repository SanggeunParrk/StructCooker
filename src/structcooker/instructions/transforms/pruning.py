"""Project a CIFMol onto the features one downstream model actually reads.

A model-input DB is a *view* of the general CIFMol DB, not a new parse of it: an
MPNN reads backbone coordinates, a residue index, a sequence token and the chain
contact graph, and pays disk + per-sample decompression for every other feature
it never touches. :func:`select_features` keeps a declared feature list and drops
the rest, including the derived half of the index table -- the CSR maps are a
pure function of the two parent maps, so storing them is redundant.

The two halves belong together: :func:`select_features` writes the pruned record
and :func:`restore_index_table` (or :func:`load_pruned_cifmol`) reads it back,
rebuilding what was dropped. A consumer that loads such a DB by hand must go
through them, because :meth:`biomol.core.index.IndexTable.from_dict` requires all
six index arrays.
"""

from __future__ import annotations

from dataclasses import fields
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from biomol.core.feature import EdgeFeature, NodeFeature
from biomol.enums import StructureLevel
from biomol.exceptions import FeatureKeyError

from structcooker.mols import CIFMolAttached

if TYPE_CHECKING:
    from collections.abc import Sequence

    from biomol.core.container import FeatureContainer
    from biomol.core.feature import Feature

    from structcooker.mols import CIFMol

_LEVELS: dict[str, StructureLevel] = {
    "atoms": StructureLevel.ATOM,
    "residues": StructureLevel.RESIDUE,
    "chains": StructureLevel.CHAIN,
}

_PARENT_FIELDS = ("atom_to_res", "res_to_chain")
_CSR_FIELDS = (
    "res_atom_indptr",
    "res_atom_indices",
    "chain_res_indptr",
    "chain_res_indices",
)


def select_features(
    cifmol: CIFMol | CIFMolAttached | None,
    atom_features: Sequence[str],
    residue_features: Sequence[str],
    chain_features: Sequence[str],
    atom_edge_features: Sequence[str] = (),
    residue_edge_features: Sequence[str] = (),
    chain_edge_features: Sequence[str] = (),
    *,
    compact_index_table: bool = True,
) -> dict | None:
    """Keep only the named features, as a BioMol dict ready to serialize.

    Parameters
    ----------
    cifmol
        The molecule to project. ``None`` propagates (the entry is dropped).
    atom_features, residue_features, chain_features
        Node feature names to keep at each level. At least one per level is
        required -- a :class:`~biomol.core.container.FeatureContainer` with no
        node features cannot be reconstructed.
    atom_edge_features, residue_edge_features, chain_edge_features
        Edge feature names to keep at each level.
    compact_index_table
        Store only the ``atom_to_res`` / ``res_to_chain`` parent maps and let the
        reader rebuild the CSR maps (see :func:`restore_index_table`). The CSR
        maps are the larger half of the index table, so keeping them roughly
        doubles a pruned record.

    Raises
    ------
    FeatureKeyError
        If any requested feature is absent, so a missing feature fails the build
        loudly instead of silently producing a DB the consumer cannot read.
    """
    if cifmol is None:
        return None

    requested = {
        "atoms": (atom_features, atom_edge_features),
        "residues": (residue_features, residue_edge_features),
        "chains": (chain_features, chain_edge_features),
    }
    data: dict[str, Any] = {}
    missing: list[str] = []
    for level, (node_names, edge_names) in requested.items():
        container = cifmol.get_container(_LEVELS[level])
        nodes = _collect(container, level, "nodes", node_names, NodeFeature, missing)
        edges = _collect(container, level, "edges", edge_names, EdgeFeature, missing)
        data[level] = {"nodes": nodes, "edges": edges}
    if missing:
        msg = f"cannot keep absent features: {', '.join(sorted(missing))}"
        raise FeatureKeyError(msg)

    index = cifmol.index_table
    fields = _PARENT_FIELDS if compact_index_table else _PARENT_FIELDS + _CSR_FIELDS
    # int32 throughout: these index the record's own atoms/residues/chains, and
    # no assembly comes close to 2**31 of any of them.
    data["index_table"] = {
        name: np.asarray(getattr(index, name), dtype=np.int32) for name in fields
    }
    data["metadata"] = dict(cifmol.metadata)
    return data


def restore_index_table(value: dict) -> dict:
    """Rebuild the CSR index maps a compact record dropped.

    Returns ``value`` unchanged when the record already carries the full index
    table, so this is safe to apply to pruned and unpruned records alike.
    """
    index = value["index_table"]
    if all(name in index for name in _CSR_FIELDS):
        return value

    atom_to_res = np.asarray(index["atom_to_res"], dtype=int)
    res_to_chain = np.asarray(index["res_to_chain"], dtype=int)
    # The chain count comes from the chain container, not from res_to_chain: a
    # chain holding no residues is never named there, and inferring the count
    # from the maps alone would silently drop it.
    n_chain = _node_count(value["chains"])
    res_atom_indptr, res_atom_indices = _build_csr(atom_to_res, len(res_to_chain))
    chain_res_indptr, chain_res_indices = _build_csr(res_to_chain, n_chain)

    restored = dict(value)
    restored["index_table"] = {
        "atom_to_res": atom_to_res,
        "res_to_chain": res_to_chain,
        "res_atom_indptr": res_atom_indptr,
        "res_atom_indices": res_atom_indices,
        "chain_res_indptr": chain_res_indptr,
        "chain_res_indices": chain_res_indices,
    }
    return restored


def load_pruned_cifmol(value: dict) -> CIFMolAttached:
    """Deserialize one pruned record into a CIFMolAttached."""
    return CIFMolAttached.from_dict(cast("Any", restore_index_table(value)))


def convert_to_pruned_cifmol(value: dict) -> dict[str, CIFMolAttached]:
    """Reader adapter for a pruned DB: one record -> ``{"cifmol": CIFMolAttached}``."""
    return {"cifmol": load_pruned_cifmol(value)}


def _collect(
    container: FeatureContainer,
    level: str,
    kind: str,
    names: Sequence[str],
    expected: type[Feature],
    missing: list[str],
) -> dict[str, dict[str, Any]]:
    """Pull the named features of one kind out of a container.

    Names that are absent -- or present at the other kind, e.g. an edge feature
    asked for as a node feature -- are recorded in ``missing`` rather than
    raised on, so one error names every problem in the config at once.
    """
    collected: dict[str, dict[str, Any]] = {}
    for name in names:
        if name not in container or not isinstance(container[name], expected):
            missing.append(f"{level}.{kind}.{name}")
            continue
        collected[name] = _feature_dict(container[name])
    return collected


def _feature_dict(feature: Feature) -> dict[str, Any]:
    # Read the dataclass fields rather than calling dataclasses.asdict, which
    # deep-copies every array -- pointless here, since they go straight to the
    # serializer -- and rather than naming them, which would pin this to one
    # biomol version's Feature layout.
    return {field.name: getattr(feature, field.name) for field in fields(feature)}


def _node_count(level: dict) -> int:
    nodes = level["nodes"]
    if not nodes:
        msg = "pruned record has a level with no node features"
        raise FeatureKeyError(msg)
    return len(next(iter(nodes.values()))["value"])


def _build_csr(
    parent_of_child: np.ndarray,
    n_parent: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized twin of ``biomol.core.index._build_csr``.

    biomol's version loops in Python over every child, which is a real cost when
    it runs once per training sample; this produces the same arrays -- children
    grouped by parent, ascending within a group -- with numpy.
    """
    counts = np.bincount(parent_of_child, minlength=n_parent)
    indptr = np.empty(n_parent + 1, dtype=int)
    indptr[0] = 0
    np.cumsum(counts, dtype=int, out=indptr[1:])
    indices = np.argsort(parent_of_child, kind="stable")
    return indptr, indices
