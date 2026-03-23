from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.convert import load_fasta
from pipelines.instructions.metadata_instructions import (
    build_seq_metadata_map,
    build_template_metadata_map,
    load_signalp,
    load_tsv,
    parse_metadata,
)

"""Rebuild a CIF lmdb to train AF3"""

recipe = RecipeBook()

recipe.add(
    targets=[
        (("pdb_id2deposition_date", dict),),
    ],
    instruction=parse_metadata,
    inputs=[
        {
            "kwargs": {
                "metadata_path": ("metadata_path", Path),
            },
        },
    ],
)


recipe.add(
    targets=[
        (("raw_fasta_dict", dict),),
    ],
    instruction=load_fasta,
    inputs=[
        {
            "kwargs": {
                "fasta_path": ("raw_fasta_path", str | Path),
            },
        },
    ],
)

recipe.add(
    targets=[
        (("seqid2seq", dict),),
        (("seqclusters2seqids", dict),),
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
    ],
)

recipe.add(
    targets=[
        (("seq_metadata_map", dict),),  # cif_id -> seq id + seq cluster
    ],
    instruction=build_seq_metadata_map,
    inputs=[
        {
            "kwargs": {
                "raw_fasta_dict": ("raw_fasta_dict", dict),
                "seqid2seq": ("seqid2seq", dict),
                "seqclusters2seqids": ("seqclusters2seqids", dict),
            },
        },
    ],
)


recipe.add(
    targets=[(("signalp_dict", dict),)],
    instruction=load_signalp,
    inputs=[
        {
            "kwargs": {
                "signalp_dir": ("signalp_dir", Path | None),
            },
        },
    ],
)


recipe.add(
    targets=[
        (
            ("template_metadata_map", dict),
            ("seqid2earliest_date", dict),
            ("filtered_seqid2seq", dict),
        ),
    ],
    instruction=build_template_metadata_map,
    inputs=[
        {
            "kwargs": {
                "pdb_id2deposition_date": ("pdb_id2deposition_date", dict),
                "seq_metadata_map": ("seq_metadata_map", dict),
                "seqid2seq": ("seqid2seq", dict),
                "signalp_dict": ("signalp_dict", dict | None),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["template_metadata_map", "seqid2earliest_date", "filtered_seqid2seq"]
