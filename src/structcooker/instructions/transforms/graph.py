import itertools
from collections import Counter, defaultdict
from collections.abc import Callable
from typing import TypeVar, cast

import networkx as nx
from joblib import Parallel, delayed

from structcooker.instructions.transforms.sequence import (
    extract_sequence_from_cifmol,
    filter_water,
)
from structcooker.mols import CIFMol
from structcooker.utils.mapping import cluster_id, mol_type_map

InputType = TypeVar("InputType", str, int, float)
FeatureType = TypeVar("FeatureType")
NumericType = TypeVar("NumericType", int, float)


def build_graph_hash(
    n_jobs: int = -1,
) -> Callable[..., tuple[dict[int, nx.Graph], dict[str | int, int]]]:
    """Deduplicate graphs up to isomorphism using node 'label' attribute."""

    def _worker(
        graph_map: dict[str, nx.Graph],
    ) -> tuple[dict[int, nx.Graph], dict[str | int, int]]:
        """
        Deduplicate graphs up to isomorphism using node 'label' attribute.

        - Uses strong invariants (WL hash + degree multiset + label multiset) to bucket graphs.
        - Within each bucket, runs exact isomorphism checks with categorical node match.
        - Assigns contiguous global cluster IDs starting from `start_idx`.

        Args:
            graph_map: Mapping from user graph key -> NetworkX Graph.
                    Each node must have a 'label' attribute (hashable).
            start_idx: Starting integer for cluster IDs.
            n_jobs:    Parallel workers for invariant computation; use -1 for all cores.

        Returns
        -------
            unique_map: {cluster_id -> representative Graph}
            gid_to_id:  {input_key -> cluster_id}

        Notes
        -----
            - Bucketing by invariants prunes most non-isomorphic pairs before exact checks.
            - Exact checks use GraphMatcher with categorical_node_match('label', None).
            - Determinism: per-bucket keys are processed in sorted order of input keys.
        """
        if not graph_map:
            return {}, {}

        # ---------- 1) Compute invariants in parallel ----------
        # Invariants are strong but cheap: WL hash (label-aware), degree multiset, label multiset.
        def _invariants(item: tuple[str | int, nx.Graph]) -> tuple[str | int, tuple]:
            gid, G = item
            # WL hash aware of node label
            wl = nx.weisfeiler_lehman_graph_hash(G, node_attr="label")
            # degree multiset
            degs = tuple(sorted(d for _, d in G.degree))
            # node-label multiset
            lbls = tuple(sorted(Counter(nx.get_node_attributes(G, "label").values()).items()))
            return gid, (wl, degs, lbls)

        items = list(graph_map.items())
        inv_list = Parallel(n_jobs=n_jobs, verbose=5)(
            delayed(_invariants)(it) for it in items
        )
        inv_map: dict[str | int, tuple] = dict(
            cast("list[tuple[str | int, tuple]]", inv_list),
        )

        # ---------- 2) Bucket by invariants ----------
        buckets: dict[tuple, list[str | int]] = defaultdict(list)
        for gid, key in inv_map.items():
            buckets[key].append(gid)

        # Sort keys inside each bucket for determinism
        for values in buckets.values():
            values.sort(key=lambda x: (str(type(x)), str(x)))

        # ---------- 3) Per-bucket exact clustering (isomorphism with labels) ----------
        # Use threading backend to avoid heavy pickling of NetworkX graphs.
        node_match = nx.algorithms.isomorphism.categorical_node_match("label", None)

        def _cluster_bucket(
            gids: list[str | int],
        ) -> list[tuple[list[str | int], str | int]]:
            """
            Cluster one bucket into isomorphism classes.

            Returns a list of clusters as (member_gids, representative_gid).
            """
            reps: list[tuple[str | int, nx.Graph]] = []  # [(rep_gid, rep_graph)]
            clusters: list[list[str | int]] = []

            for gid in gids:
                G = graph_map[cast("str", gid)]
                placed = False
                # Try match against existing representatives
                for c_idx, (_, repG) in enumerate(reps):
                    GM = nx.isomorphism.GraphMatcher(G, repG, node_match=node_match)
                    if GM.is_isomorphic():
                        clusters[c_idx].append(gid)
                        placed = True
                        break
                if not placed:
                    reps.append((gid, G))
                    clusters.append([gid])

            # Pair each cluster with its representative gid (first in deterministic order)
            return [(cls, reps[i][0]) for i, cls in enumerate(clusters)]

        bucket_results = Parallel(n_jobs=n_jobs, backend="threading", verbose=5)(
            delayed(_cluster_bucket)(gids) for _, gids in buckets.items()
        )

        # Flatten list of clusters
        clusters_all: list[tuple[list[str | int], str | int]] = list(
            itertools.chain.from_iterable(
                cast("list[list[tuple[list[str | int], str | int]]]", bucket_results),
            ),
        )

        # ---------- 4) Assign global cluster IDs and build outputs ----------
        unique_map: dict[int, nx.Graph] = {}
        gid_to_id: dict[str | int, int] = {}

        for cur, (members, rep_gid) in enumerate(clusters_all):
            unique_map[cur] = graph_map[cast("str", rep_gid)]
            for gid in members:
                gid_to_id[gid] = cur
        return unique_map, gid_to_id

    return _worker


def interacting_seq_ids(
    cifmol_dict: dict[str, dict[str, CIFMol]],
    seqid_map: dict,
) -> set[tuple[str, str]]:
    """Extract interacting sequence IDs from cifmol dict."""
    interacting_seq_ids: set = set()
    for cifmol_wrapper in cifmol_dict.values():
        cifmol = cifmol_wrapper["cifmol"]
        cifmol_wo_water = filter_water(cifmol)
        if cifmol_wo_water is None:
            continue
        seq_dict, entity_type_dict = extract_sequence_from_cifmol(cifmol_wo_water)
        seq_id_dict = {}

        for chain_id, seq in seq_dict.items():
            entity_type = entity_type_dict[chain_id]
            mol_identifier = mol_type_map.get(entity_type, "X")
            seq_id = seqid_map[mol_identifier][seq]
            seq_id_dict[chain_id] = seq_id
        contact_graph = cifmol_wo_water.chains.contact
        for src_idx, dst_idx in zip(
            contact_graph.src_indices,
            contact_graph.dst_indices,
            strict=False,
        ):
            src_chain = cifmol_wo_water.chains.chain_id.value[src_idx]
            dst_chain = cifmol_wo_water.chains.chain_id.value[dst_idx]
            src_seq_id = seq_id_dict[src_chain]
            dst_seq_id = seq_id_dict[dst_chain]
            pair_key = sorted((src_seq_id, dst_seq_id))
            interacting_seq_ids.add((pair_key[0], pair_key[1]))

    return interacting_seq_ids


def interacting_seq_ids_from_attached(record: dict) -> set[tuple[str, str]]:
    """Read interacting sequence pairs from attached chain IDs and contact edges.

    The attach stage has already removed water and assigned sequence IDs on the
    full sequences. No coordinate reconstruction or global sequence lookup is needed.
    """
    pairs: set[tuple[str, str]] = set()
    for wrapped in record.values():
        chains = wrapped["cifmol_attached_dict"]["chains"]
        seqids = chains["nodes"]["seq_id"]["value"]
        contact = chains["edges"]["contact"]
        for src, dst in zip(contact["src_indices"], contact["dst_indices"], strict=True):
            left, right = sorted((str(seqids[src]), str(seqids[dst])))
            pairs.add((left, right))
    return pairs


def filter_seq_ids(
    interacting_seq_ids: set[tuple[str, str]],
    *,
    valid_entity_types: set[str],
) -> set[tuple[str, str]]:
    """Filter interacting sequence IDs by valid entity types."""
    filtered_seq_ids: set = set()
    for seq_id1, seq_id2 in interacting_seq_ids:
        entity_type1 = seq_id1[0]
        entity_type2 = seq_id2[0]
        if entity_type1 in valid_entity_types or entity_type2 in valid_entity_types:
            filtered_seq_ids.add((seq_id1, seq_id2))
    return filtered_seq_ids


def build_interacting_seq_clusters(
    interacting_seq_ids: dict[str, list[str]],
    seqclusters2seqids: dict,
    db_code: str,
) -> set[tuple[str, str]]:
    """Build interacting sequence clusters from interacting sequence IDs."""
    interacting_seq_clusters: set = set()
    interacting_seq_ids_pair = {
        (key, partner) for key, partners in interacting_seq_ids.items()
        for partner in partners
    }
    seqid2cluster = {}
    for cluster, seqids in seqclusters2seqids.items():
        for seqid in seqids:
            seqid2cluster[seqid] = cluster
    for seq_id1, seq_id2 in interacting_seq_ids_pair:
        # Match attachment and validation classification: a sequence absent from
        # the supplied clustering keeps its own singleton cluster. Keep its edges
        # so a contact with training data still excludes it from validation.
        cluster1 = seqid2cluster.get(seq_id1, cluster_id(db_code, seq_id1))
        cluster2 = seqid2cluster.get(seq_id2, cluster_id(db_code, seq_id2))
        pair_key = sorted((cluster1, cluster2))
        interacting_seq_clusters.add((pair_key[0], pair_key[1]))
    return interacting_seq_clusters
