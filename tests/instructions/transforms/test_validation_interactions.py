from structcooker.instructions.transforms.filtering import filter_valid_2_clusters
from structcooker.instructions.transforms.graph import build_interacting_seq_clusters
from structcooker.instructions.transforms.metadata import load_pairs


def test_all_partners_survive_loading_and_cluster_mapping(tmp_path):
    path=tmp_path/"pairs.tsv"
    path.write_text("P1\tP2\nP1\tP3\n")
    pairs=load_pairs(path)
    assert pairs == {"P1":["P2","P3"]}
    assert build_interacting_seq_clusters(pairs,{"cPDB_P1":["P1"],"cPDB_P2":["P2"],"cPDB_P3":["P3"]},"PDB") \
        == {("cPDB_P1","cPDB_P2"),("cPDB_P1","cPDB_P3")}


def test_transitive_train_contacts_exclude_validation_clusters():
    remaining=filter_valid_2_clusters({"cP1"},{"cP2","cP3","cP4","cL1"},
                                     {"cP1":["cP2"],"cP2":["cP3"],"cP4":["cP5"]})
    assert remaining == {"cP4","cL1"}


def test_unclustered_sequences_keep_contacts_to_training():
    edges = build_interacting_seq_clusters(
        {"Pnew": ["Pknown", "Pother"]}, {"cPDB_Ptrain": ["Pknown"]}, "PDB",
    )
    # Unclustered ids fall back to their own singleton, in this DB's cluster namespace.
    assert edges == {("cPDB_Pnew", "cPDB_Ptrain"), ("cPDB_Pnew", "cPDB_Pother")}
    partners = {}
    for left, right in edges:
        partners.setdefault(left, []).append(right)
    assert filter_valid_2_clusters({"cPDB_Ptrain"}, {"cPDB_Pnew", "cPDB_Pother"}, partners) == set()
