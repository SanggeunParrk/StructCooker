from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.convert import load_fasta
from pipelines.instructions.filter_instructions import (
    filter_valid_2_clusters,
)
from pipelines.instructions.metadata_instructions import (
    build_seqid_map,
    classify_seq_clusters,
    load_tsv,
)

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()

recipe.add(
    targets=[
        (("raw_fasta_dict", dict),),
        (("train_fasta_dict", dict),),
        (("valid_1_fasta_dict", dict),),
    ],
    instruction=load_fasta,
    inputs=[
        {
            "kwargs": {
                "fasta_path": ("raw_fasta_path", str | Path),
            },
        },
        {
            "kwargs": {
                "fasta_path": ("train_fasta_path", str | Path),
            },
        },
        {
            "kwargs": {
                "fasta_path": ("valid_1_fasta_path", str | Path),
            },
        },
    ],
)


recipe.add(
    targets=[
        (("seqid2seq", dict),),
        (("seqclusters2seqids", dict),),
        (("interacting_seq_clusters", dict),),
    ],
    instruction=load_tsv,
    inputs=[
        {
            "kwargs": {
                "tsv_file_path": ("seqid2seq_path", Path),
            },
            "params": {
                "split_by_comma": False,
            },
        },
        {
            "kwargs": {
                "tsv_file_path": ("seqcluster_path", Path),
            },
            "params": {
                "split_by_comma": True,
            },
        },
        {
            "kwargs": {
                "tsv_file_path": ("interacting_seq_clusters_path", Path),
            },
            "params": {
                "split_by_comma": False,
            },
        },
    ],
)

recipe.add(
    targets=[
        (("seqid_map", dict),),  # moltype+seq -> seqid
    ],
    instruction=build_seqid_map,
    inputs=[
        {
            "kwargs": {
                "seqid2seq": ("seqid2seq", dict),
            },
        },
    ],
)

recipe.add(
    targets=[
        (("train_clusters", set),),
        (("valid_1_clusters", set),),
    ],
    instruction=classify_seq_clusters,
    inputs=[
        {
            "kwargs": {
                "raw_fasta_dict": ("raw_fasta_dict", dict),
                "fasta_dict": ("train_fasta_dict", dict),
                "seqid_map": ("seqid_map", dict),
                "seqclusters2seqids": ("seqclusters2seqids", dict),
            },
        },
        {
            "kwargs": {
                "raw_fasta_dict": ("raw_fasta_dict", dict),
                "fasta_dict": ("valid_1_fasta_dict", dict),
                "seqid_map": ("seqid_map", dict),
                "seqclusters2seqids": ("seqclusters2seqids", dict),
            },
        },
    ],
)

recipe.add(
    targets=[
        (("filtered_valid_2_clusters", set),),
    ],
    instruction=filter_valid_2_clusters,
    inputs=[
        {
            "kwargs": {
                "train_clusters": ("train_clusters", set),
                "valid_1_clusters": ("valid_1_clusters", set),
                "interacting_seq_clusters": ("interacting_seq_clusters", dict),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["filtered_valid_2_clusters"]
