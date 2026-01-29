import itertools
from collections import Counter, OrderedDict, defaultdict
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import TypeVar

import networkx as nx
import numpy as np
from biomol.core.container import FeatureContainer
from biomol.core.feature import EdgeFeature, NodeFeature
from joblib import Parallel, delayed

from pipelines.cifmol import CIFMol, CIFMolAttached
from pipelines.instructions.seq_instructions import graph_to_canonical_sequence
from pipelines.utils.convert import to_bytes

InputType = TypeVar("InputType", str, int, float)
FeatureType = TypeVar("FeatureType")
NumericType = TypeVar("NumericType", int, float)


def extract_graph_per_cifmol() -> Callable[..., bytes]:
    """Read CIFMol and extract chain-level contact graphs with cluster labels."""

    def _worker(
        cifmol: CIFMol,
        seq_to_seq_hash_list: dict[str, list[str]],
        seq_hash_to_cluster: dict[str, str],
    ) -> type[InputType]:
        chain_id_to_cluster = {}
        chain_ids = cifmol.chains.chain_id.value

        for full_chain_id in chain_ids:
            entity_type = cifmol.chains[
                cifmol.chains.chain_id == full_chain_id
            ].entity_type.value[0]
            if entity_type == "non-polymer":
                seq = cifmol.chains[
                    cifmol.chains.chain_id == full_chain_id
                ].residues.chem_comp_id.value
                seq = f"({seq[0]})"
            elif entity_type == "branched":
                seq_list = cifmol.chains[
                    cifmol.chains.chain_id == full_chain_id
                ].residues.chem_comp_id.value
                bonds = cifmol.chains[
                    cifmol.chains.chain_id == full_chain_id
                ].residues.bond
                seq = graph_to_canonical_sequence(
                    seq_list,
                    bonds.src_indices,
                    bonds.dst_indices,
                )
            else:  # polymer
                seq = cifmol.chains[
                    cifmol.chains.chain_id == full_chain_id
                ].residues.one_letter_code_can.value
                seq = "".join(seq)

            match entity_type:
                case "polypeptide(L)":
                    mol_identifier = "P"
                case "polypeptide(D)":
                    mol_identifier = "Q"
                case "polydeoxyribonucleotide":
                    mol_identifier = "D"
                case "polyribonucleotide":
                    mol_identifier = "R"
                case "polydeoxyribonucleotide/polyribonucleotide hybrid":
                    mol_identifier = "N"
                case "branched":
                    mol_identifier = "B"
                case "non-polymer":
                    mol_identifier = "L"
                case _:
                    mol_identifier = "X"

            seq_hash_list = seq_to_seq_hash_list[seq]
            seq_hash = next(h for h in seq_hash_list if h.startswith(mol_identifier))
            chain_id_to_cluster[full_chain_id] = seq_hash_to_cluster[seq_hash]

        # chain level contact graph
        contact_graph = cifmol.chains.contact
        chain_id_list = cifmol.chains.chain_id.value
        cluster_list = [chain_id_to_cluster[chain_id] for chain_id in chain_id_list]
        src, dst = contact_graph.src_indices, contact_graph.dst_indices

        cluster_list = np.array(cluster_list)
        features = {
            "chain_ids": NodeFeature(np.array(chain_id_list)),
            "seq_clusters": NodeFeature(cluster_list),
            "contact_edges": EdgeFeature(
                value=np.array([1] * len(src)),
                src_indices=np.array(src),
                dst_indices=np.array(dst),
            ),
        }
        container = FeatureContainer(features)
        return to_bytes({"cluster_graph": container})

    return _worker



def extract_graph_per_cifmol_attached() -> Callable[..., bytes]:
    """Read CIFMolAttached and extract chain-level contact graphs with cluster labels."""

    def _worker(
        cifmol: CIFMolAttached,
    ) -> type[InputType]:
        # chain level contact graph
        contact_graph = cifmol.chains.contact
        chain_id_list = cifmol.chains.chain_id.value
        cluster_list = cifmol.chains.cluster_id.value
        src, dst = contact_graph.src_indices, contact_graph.dst_indices

        cluster_list = np.array(cluster_list)
        features = {
            "chain_ids": NodeFeature(np.array(chain_id_list)),
            "seq_clusters": NodeFeature(cluster_list),
            "contact_edges": EdgeFeature(
                value=np.array([1] * len(src)),
                src_indices=np.array(src),
                dst_indices=np.array(dst),
            ),
        }
        container = FeatureContainer(features)
        return to_bytes({"cluster_graph": container})

    return _worker


def extract_graphs(n_jobs: int = -1) -> Callable[..., type[InputType]]:
    """Read CIFMol and extract chain-level contact graphs with cluster labels."""

    def _single_function(
        cifmol: CIFMol,
        seq_to_cluster: dict[str, str],
    ) -> type[InputType]:
        chain_id_to_cluster = {}
        chain_ids = cifmol.chains.chain_id.value
        for full_chain_id in chain_ids:
            entity_type = cifmol.chains[
                cifmol.chains.chain_id == full_chain_id
            ].entity_type.value[0]
            if entity_type == "non-polymer":
                seq = cifmol.chains[
                    cifmol.chains.chain_id == full_chain_id
                ].residues.chem_comp_id.value
                seq = f"({seq[0]})"
            elif entity_type == "branched":
                seq_list = cifmol.chains[
                    cifmol.chains.chain_id == full_chain_id
                ].residues.chem_comp_id.value
                bonds = cifmol.chains[
                    cifmol.chains.chain_id == full_chain_id
                ].residues.bond
                seq = graph_to_canonical_sequence(
                    seq_list,
                    bonds.src_indices,
                    bonds.dst_indices,
                )
            else:  # polymer
                seq = cifmol.chains[
                    cifmol.chains.chain_id == full_chain_id
                ].residues.one_letter_code_can.value
                seq = "".join(seq)
            chain_id_to_cluster[full_chain_id] = seq_to_cluster[seq]

        # chain level contact graph
        contact_graph = cifmol.chains.contact
        chain_id_list = cifmol.chains.chain_id.value
        src, dst = contact_graph.src_indices, contact_graph.dst_indices
        graph = nx.Graph()
        for chain_id in chain_id_list:
            graph.add_node(chain_id, label=chain_id_to_cluster[chain_id])
        for s, d in zip(src, dst, strict=True):
            graph.add_edge(chain_id_list[s], chain_id_list[d])

        return graph

    def _worker(
        cifmol_dict: dict[str, CIFMol],
        seq_hash_map: Path,  # seq to seq hash
        seq_cluster_map: Path,  # seq hash to cluster
    ) -> dict[str, nx.Graph]:
        seq_to_hash: dict[str, int] = {}
        with seq_hash_map.open("r") as f:
            for line in f:
                seq_hash, seq = line.strip().split("\t")
                seq_to_hash[seq] = seq_hash
        seq_hash_to_cluster: dict[int, str] = {}
        with seq_cluster_map.open("r") as f:
            for line in f:
                rep, members = line.strip().split("\t")
                members = members.split(",")
                for m in members:
                    seq_hash_to_cluster[m] = rep

        seq_to_cluster: dict[str, str] = {}
        for seq, seq_hash in seq_to_hash.items():
            cluster_id = seq_hash_to_cluster[seq_hash]
            seq_to_cluster[seq] = cluster_id
        results = Parallel(n_jobs=n_jobs, verbose=10)(
            delayed(_single_function)(cifmol, seq_to_cluster)
            for cifmol in cifmol_dict.values()
        )
        return dict(zip(cifmol_dict.keys(), results, strict=True))

    return _worker


def has_common_node(g1: nx.Graph, g2: nx.Graph) -> bool:
    """Check if two graphs have common node labels."""
    labels_g1 = {data.get("label") for node, data in g1.nodes(data=True)}
    labels_g2 = {data.get("label") for node, data in g2.nodes(data=True)}
    return bool(labels_g1.intersection(labels_g2))


def get_edge_labels(graph: nx.Graph) -> set[tuple[str, str]]:
    """Get edge labels from a graph."""
    edge_labels = set()
    for u, v in graph.edges():
        label_u = graph.nodes[u].get("label")
        label_v = graph.nodes[v].get("label")
        edge_labels.add(tuple(sorted((label_u, label_v))))
    return edge_labels


def has_common_edge(g1: nx.Graph, g2: nx.Graph) -> bool:
    """Check if two graphs have common edge labels."""
    edge_labels_g1 = get_edge_labels(g1)
    edge_labels_g2 = get_edge_labels(g2)
    return bool(edge_labels_g1.intersection(edge_labels_g2))


def graph_isomorphism(graph1: nx.Graph, graph2: nx.Graph) -> bool:
    """Check if two graphs are isomorphic."""

    def node_match(attr1: dict, attr2: dict) -> bool:
        return attr1.get("label") == attr2.get("label")

    gm = nx.isomorphism.GraphMatcher(graph1, graph2, node_match=node_match)
    return gm.is_isomorphic()


def build_graph_hash(n_jobs: int = -1) -> Callable[..., type[InputType]]:
    """Deduplicate graphs up to isomorphism using node 'label' attribute."""

    def _worker(
        graph_map: dict[str, nx.Graph],
    ) -> dict[str, dict[str, str]]:
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
            degs = tuple(sorted(d for _, d in G.degree()))
            # node-label multiset
            lbls = tuple(sorted(Counter(nx.get_node_attributes(G, "label")).items()))
            return gid, (wl, degs, lbls)

        items = list(graph_map.items())
        inv_list = Parallel(n_jobs=n_jobs, verbose=5)(
            delayed(_invariants)(it) for it in items
        )
        inv_map: dict[str | int, tuple] = dict(inv_list)

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
                G = graph_map[gid]
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
            itertools.chain.from_iterable(bucket_results),
        )

        # ---------- 4) Assign global cluster IDs and build outputs ----------
        unique_map: dict[int, nx.Graph] = {}
        gid_to_id: dict[str | int, int] = {}

        for cur, (members, rep_gid) in enumerate(clusters_all):
            unique_map[cur] = graph_map[rep_gid]
            for gid in members:
                gid_to_id[gid] = cur
        return unique_map, gid_to_id

    return _worker


class _UnionFind:
    """Disjoint-set (Union-Find) with path compression and union by rank."""

    def __init__(self, elems: Iterable[int]) -> None:
        self._index = {e: i for i, e in enumerate(sorted(elems))}
        n = len(self._index)
        self.parent = list(range(n))
        self.rank = [0] * n
        self.reverse = {i: e for e, i in self._index.items()}

    def _id(self, e: int) -> int:
        return self._index[e]

    def find(self, e: int) -> int:
        i = self._id(e)
        # path compression
        if self.parent[i] != i:
            self.parent[i] = self.find(self.reverse[self.parent[i]])
        return self.parent[i]

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        # union by rank
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[rb] < self.rank[ra]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1

    def components(self) -> list[set[int]]:
        roots_to_members: dict[int, set[int]] = defaultdict(set)
        for e in self.reverse.values():
            r = self.find(e)
            roots_to_members[r].add(e)
        # stable order by smallest member for determinism
        comps = list(roots_to_members.values())
        comps.sort(key=lambda s: min(s))
        return comps


def graph_edge_cluster(n_jobs: int = -1) -> Callable[..., type[InputType]]:
    """Cluster graphs that share at least one edge defined by node label pairs."""

    def _worker(
        unique_graphs: dict[int, nx.Graph],
    ) -> tuple[list[set[int]], dict[int, int]]:
        """
        Cluster graphs if they share at least one edge defined by node label pairs.

        An "edge key" is an unordered pair (min(label_u, label_v), max(label_u, label_v)).
        Two graphs belong to the same cluster if they share >= 1 such edge key.

        Args:
            unique_graphs: {graph_id -> NetworkX Graph} where each node has attribute 'label'.
            n_jobs:        Number of parallel workers for edge extraction (-1 = all cores).

        Returns
        -------
            clusters: List of sets of graph_ids, one set per connected component.
            comp_map: Mapping {graph_id -> component_index} where components are 0..C-1
                    assigned in ascending order of each component's smallest graph_id.

        Notes
        -----
            - Graphs with no edges (or single-node) become singleton clusters.
            - Uses threading backend implicitly (joblib default) to avoid heavy pickling.
            - Deterministic component indices via sorting by smallest member id.
        """
        if not unique_graphs:
            return [], {}

        graph_ids = sorted(unique_graphs.keys())

        # --- 1) Extract per-graph labeled edge keys in parallel -------------------
        def _extract_edge_keys(gid: int) -> tuple[int, list[tuple[object, object]]]:
            G = unique_graphs[gid]
            labels = nx.get_node_attributes(G, "label")
            # Build a set to avoid duplicates within the same graph
            keys: set[tuple[object, object]] = set()
            # If directed graphs appear, treat as undirected by sorting endpoints
            for u, v in G.edges():
                lu = labels.get(u, None)
                lv = labels.get(v, None)
                # Require both endpoints to have 'label'
                if lu is None or lv is None:
                    continue
                key = (lu, lv) if lu <= lv else (lv, lu)
                keys.add(key)
            return gid, sorted(keys)

        edge_key_lists: list[tuple[int, list[tuple[object, object]]]] = Parallel(
            n_jobs=n_jobs, verbose=5,
        )(delayed(_extract_edge_keys)(gid) for gid in graph_ids)

        # --- 2) Build edge_key -> [graph_id, ...] inverted index ------------------
        edge_to_graphs: dict[tuple[object, object], list[int]] = defaultdict(list)
        for gid, keys in edge_key_lists:
            for k in keys:
                edge_to_graphs[k].append(gid)

        # --- 3) Union-Find: union all graphs that share any edge key --------------
        uf = _UnionFind(graph_ids)
        for gids in edge_to_graphs.values():
            if len(gids) <= 1:
                continue
            base = gids[0]
            for other in gids[1:]:
                uf.union(base, other)

        # --- 4) Collect components and build comp_map -----------------------------
        clusters = uf.components()  # sorted by smallest member id
        comp_map: dict[int, int] = {}
        for c_idx, members in enumerate(clusters):
            for gid in members:
                comp_map[gid] = c_idx

        return clusters, comp_map

    return _worker


def convert_graph_to_bytes(n_jobs: int = -1) -> Callable[..., bytes]:
    """Convert a NetworkX graph to a FeatureContainer."""

    def _function(
        graph: nx.Graph,
    ) -> bytes:
        label_list = []
        src_list, dst_list = [], []
        for _, data in graph.nodes(data=True):
            label = data.get("label")
            label_list.append(label)
        for u, v in graph.edges():
            src_list.append(u)
            dst_list.append(v)

        label_list = np.array(label_list)
        features = {
            "seq_clusters": NodeFeature(label_list),
            "contact_edges": EdgeFeature(
                value=np.array([1] * len(src_list)),
                src_indices=np.array(src_list),
                dst_indices=np.array(dst_list),
            ),
        }
        container = FeatureContainer(features)
        return to_bytes({"cluster_graph": container})

    def _worker(
        graph_map: dict[str, nx.Graph],
    ) -> dict[str, bytes]:
        results = Parallel(n_jobs=n_jobs, verbose=10)(
            delayed(_function)(graph) for graph in graph_map.values()
        )
        return dict(zip(graph_map.keys(), results, strict=True))

    return _worker


def build_whole_graph(edge_tsv_path: Path, ignore_nodes:list|None=None) -> tuple[nx.Graph, nx.Graph]:
    """Build a whole graph from edge TSV file."""
    edges = []
    polymer_edges = []
    ligand_tags = {"L", "B", "X"}
    with edge_tsv_path.open("r", encoding="utf-8") as f:
        for line in f:
            src, dst = line.strip().split("\t")[:2]
            if ignore_nodes and (src in ignore_nodes or dst in ignore_nodes):
                continue
            # test
            if src[1] not in ligand_tags and dst[1] not in ligand_tags:
                polymer_edges.append((src, dst))
            edges.append((src, dst))
    whole, polymer = nx.Graph(), nx.Graph()
    whole.add_edges_from(edges)
    polymer.add_edges_from(polymer_edges)
    return whole, polymer


def build_category()-> tuple[dict[str, int], dict[str, str], dict[str, int]]:
    """Build category pairs and priority map."""
    special_map = {"B" : "L", "Q": "P"}
    category_priority = ["P", "A", "D", "R", "N", "L"]
    priority_map = {c: i for i, c in enumerate(category_priority)}

    # 1. Create all possible ordered pairs by priority
    pair_keys = []
    for c1 in category_priority:
        for c2 in category_priority:
            # normalize: the higher priority (smaller index) goes to the BACK
            # so we force order by comparing priority index
            key = c2 + c1 if priority_map[c1] > priority_map[c2] else c1 + c2
            if key not in pair_keys:
                pair_keys.append(key)

    counts = dict.fromkeys(pair_keys, 0)
    return counts, special_map, priority_map


def count_category_count(
    edge_list: list[tuple[str, str]],
) -> list[str]:
    """Summarize split results into edge list strings."""
    counts, special_map, priority_map = build_category()

    # 2. Count edges
    for src, dst in edge_list:
        c1, c2 = src[1], dst[1]
        if c1 in special_map:
            c1 = special_map[c1]
        if c2 in special_map:
            c2 = special_map[c2]
        if c1 not in priority_map or c2 not in priority_map:
            continue
        # normalize AP, PA → PA
        key = c2 + c1 if priority_map[c1] > priority_map[c2] else c1 + c2
        counts[key] += 1

    # 3. OrderedDict
    ordered = OrderedDict()
    for k in counts:
        ordered[k] = counts[k]

    return ordered


def split_graph_by_components(
    whole_graph: nx.Graph,
) -> tuple[list[nx.Graph], dict[str, int]]:
    """Split a whole graph into connected components."""
    edge_wo_ligand_counts, _, _ = build_category()

    components = list(nx.connected_components(whole_graph))
    subgraphs = []
    for comp in components:
        subgraph = whole_graph.subgraph(comp).copy()

        subgraph_edge_list = list(subgraph.edges())
        if not subgraph_edge_list:
            continue  # skip empty subgraphs
        _counts = count_category_count(subgraph_edge_list)
        for k, v in _counts.items():
            edge_wo_ligand_counts[k] += v

        subgraphs.append(subgraph)
    # sort by size descending
    subgraphs.sort(key=lambda g: g.number_of_nodes(), reverse=True)

    return subgraphs, edge_wo_ligand_counts


def split_train_valid(
    whole_graph: nx.Graph,
    edge_wo_ligand_counts: int,
    subgraphs: list[nx.Graph],
    train_ratio: float = 0.9,  # fraction of train set
    min_valid_size: int = 100,  # min #edges in valid set
) -> tuple[list[nx.Graph], list[nx.Graph]]:
    """
    Split subgraphs into train and valid sets.

    Train:
        - Always includes the largest component
        - Then adds smallest components (by #edges) one by one
          until total #edges in train <= total_edges * train_ratio for each types

    Valid:
        - All remaining components
    """
    if not subgraphs:
        return [], []

    min_edge_counts = {k: min(int(v * (1 - train_ratio)), min_valid_size) for k, v in edge_wo_ligand_counts.items()} # min edges per category in valid set

    # test (DD except)
    remaining_edges = edge_wo_ligand_counts.copy()

    # find must-have train components and add to train set
    train_indices = set()
    train_edge_count = dict.fromkeys(remaining_edges, 0)
    for ii, subg in enumerate(subgraphs):
        comp_edges = count_category_count(list(subg.edges()))
        must_have = {k : comp_edges[k] > min_edge_counts[k] for k in edge_wo_ligand_counts}
        if any(must_have.values()):
            train_indices.add(ii)
            remaining_edges = {k: remaining_edges[k] - comp_edges[k] for k in remaining_edges}
            train_edge_count = {k: train_edge_count[k] + comp_edges[k] for k in comp_edges}

    remaining_indices = list(range(1, len(subgraphs)))
    remaining_indices.sort(key=lambda i: subgraphs[i].number_of_edges())

    for idx in remaining_indices:
        comp_edges = count_category_count(list(subgraphs[idx].edges()))
        _remaining_edges = {k: remaining_edges[k] - comp_edges[k] for k in remaining_edges}
        ok = {k : _remaining_edges[k] >= min_edge_counts[k] for k in remaining_edges}
        train_edge_count = {k: train_edge_count[k] + comp_edges[k] for k in comp_edges}
        if not all(ok.values()):
            continue
        train_indices.add(idx)
        remaining_edges = _remaining_edges

    train_subgraphs = [subgraphs[i] for i in sorted(train_indices)]
    valid_subgraphs = [subgraphs[i] for i in range(len(subgraphs)) if i not in train_indices]

    train_nodes: set = set().union(*(g.nodes for g in train_subgraphs))
    valid_nodes: set = set().union(*(g.nodes for g in valid_subgraphs))

    # get train & valid edges from whole graph (at least one node in train/valid)
    train_edges = []
    valid_edges = []
    for u, v in whole_graph.edges():
        if u in train_nodes or v in train_nodes:
            train_edges.append((u, v))
        if u in valid_nodes or v in valid_nodes:
            valid_edges.append((u, v))
        if u in train_nodes and v in valid_nodes:
            msg = f"Edge ({u}, {v}) connects train and valid nodes!"
            raise ValueError(msg)
        if v in train_nodes and u in valid_nodes:
            msg = f"Edge ({u}, {v}) connects train and valid nodes!"
            raise ValueError(msg)
    train_edges = list(set(train_edges))
    valid_edges = list(set(valid_edges))

    return train_edges, valid_edges

def extract_edges(edge_tsv_path: Path, to_be_extracted: list[tuple[str, str]]) -> list[str]:
    """Extract edges from edge TSV file."""
    extracted_edges = []
    to_be_extracted = set(to_be_extracted)
    with edge_tsv_path.open("r", encoding="utf-8") as f:
        for line in f:
            src, dst = line.strip().split("\t")[:2]
            if (src, dst) in to_be_extracted or (dst, src) in to_be_extracted:
                extracted_edges.append(line)
    return extracted_edges
