from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.sequence import (
    antibody_cluster,
    merge_cluster,
    protein_cluster,
    separate_sequences,
)

"""Build a sequence clustering Cooker.

Clustering runs once per DB (docs/seq-id-and-cluster-scheme.md), so the DB code travels
through to the writer, which stamps it into every cluster id: ``c{DB}_{rep seq_id}``.
"""

seq_cluster_recipe = RecipeBook()

seq_cluster_recipe.step(
    outputs=(("fasta_path_dict", dict),),
    instruction=separate_sequences,
    kwargs={
        "tmp_dir": ("tmp_dir", Path),
        "seq_id_map_path": ("seq_id_map_path", Path),
        "fasta_path": ("fasta_path", Path),
        "sabdab_summary_path": ("sabdab_summary_path", Path),
    },
)


seq_cluster_recipe.step(
    outputs=(
        ("antibody_cluster_dict", dict),
        ("failed_ids", set),
    ),
    instruction=antibody_cluster,
    kwargs={
        "tmp_dir": ("tmp_dir", Path),
        "fasta_path_dict": ("fasta_path_dict", Path),
    },
)

seq_cluster_recipe.step(
    outputs=(
        ("protein_cluster_dict", dict),
        ("protein_d_cluster_dict", dict),
    ),
    instruction=protein_cluster,
    kwargs={
        "tmp_dir": ("tmp_dir", Path),
        "fasta_path_dict": ("fasta_path_dict", Path),
        "failed_ids": ("failed_ids", set),
        "seq_id_map_path": ("seq_id_map_path", Path),
        "mmseqs2_seq_id": ("mmseqs2_seq_id", float | str),
        "mmseqs2_cov": ("mmseqs2_cov", float | str),
        "mmseqs2_covmode": ("mmseqs2_covmode", str),
        "mmseqs2_clustermode": ("mmseqs2_clustermode", str),
    },
)

seq_cluster_recipe.step(
    outputs=(("db_code", str),),
    instruction=lambda db_code: db_code,
    kwargs={"db_code": ("db_code", str)},
)

seq_cluster_recipe.step(
    outputs=(("cluster_dict", dict),),
    instruction=merge_cluster,
    kwargs={
        "fasta_path_dict": ("fasta_path_dict", dict),
        "protein_cluster_dict": ("protein_cluster_dict", dict),
        "protein_d_cluster_dict": ("protein_d_cluster_dict", dict),
        "antibody_cluster_dict": ("antibody_cluster_dict", dict),
    },
)

RECIPE = seq_cluster_recipe
TARGETS = ["cluster_dict", "db_code"]
