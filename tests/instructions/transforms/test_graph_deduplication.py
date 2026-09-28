import networkx as nx

from structcooker.instructions.transforms.graph import build_graph_hash


def test_isomorphic_graphs_deduplicate_after_node_renaming():
    a = nx.Graph()
    a.add_edge(0, 1)
    nx.set_node_attributes(a, {0: "protein", 1: "ligand"}, "label")
    b = nx.relabel_nodes(a, {0: 12, 1: 34})
    c = a.copy()
    c.nodes[1]["label"] = "protein"
    unique, assignments = build_graph_hash(n_jobs=1)({"a": a, "b": b, "c": c})
    assert len(unique) == 2
    assert assignments["a"] == assignments["b"]
    assert assignments["a"] != assignments["c"]
