from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.seq_cluster_instructions import (
    antibody_cluster,
    merge_cluster,
    protein_cluster,
    separate_sequences,
)

"""Build a sequence clustering Cooker."""

seq_cluster_recipe = RecipeBook()

seq_cluster_recipe.add(
    targets=[
        (("fasta_path_dict", dict),),
    ],
    instruction=separate_sequences,
    inputs=[
        {
            "kwargs": {
                "tmp_dir": ("tmp_dir", Path),
                "seq_id_map_path": ("seq_id_map_path", Path),
                "fasta_path": ("fasta_path", Path),
                "sabdab_summary_path": ("sabdab_summary_path", Path),
            },
        },
    ],
)


seq_cluster_recipe.add(
    targets=[
        (
            ("antibody_cluster_dict", dict),
            ("failed_ids", set),
        ),
    ],
    instruction=antibody_cluster,
    inputs=[
        {
            "kwargs": {
                "tmp_dir": ("tmp_dir", Path),
                "fasta_path_dict": ("fasta_path_dict", Path),
            },
        },
    ],
)

seq_cluster_recipe.add(
    targets=[
        (
            ("protein_cluster_dict", dict),
            ("protein_d_cluster_dict", dict),
        ),
    ],
    instruction=protein_cluster,
    inputs=[
        {
            "kwargs": {
                "tmp_dir": ("tmp_dir", Path),
                "fasta_path_dict": ("fasta_path_dict", Path),
                "failed_ids": ("failed_ids", set),
                "seq_id_map_path": ("seq_id_map_path", Path),
                "mmseqs2_seq_id": ("mmseqs2_seq_id", float | str),
                "mmseqs2_cov": ("mmseqs2_cov", float | str),
                "mmseqs2_covmode": ("mmseqs2_covmode", str),
                "mmseqs2_clustermode": ("mmseqs2_clustermode", str),
            },
        },
    ],
)

seq_cluster_recipe.add(
    targets=[
        (("cluster_dict", dict),),
    ],
    instruction=merge_cluster,
    inputs=[
        {
            "kwargs": {
                "fasta_path_dict": ("fasta_path_dict", dict),
                "protein_cluster_dict": ("protein_cluster_dict", dict),
                "protein_d_cluster_dict": ("protein_d_cluster_dict", dict),
                "antibody_cluster_dict": ("antibody_cluster_dict", dict),
            },
        },
    ],
)

RECIPE = seq_cluster_recipe
TARGETS = ["cluster_dict"]
