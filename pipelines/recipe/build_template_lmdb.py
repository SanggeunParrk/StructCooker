import time
from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.template_instructions import (
    extract_sequences,
    filter_and_align_template_chain_ids,
    load_templates,
)

"""Build a CIFMol->fasta Cooker."""

recipe = RecipeBook()

recipe.add(
    targets=[
        (
            ("template_chain_ids", list[str]),
            ("earliest_query_date", time.struct_time),
            ("query_seq", str),
        ),
    ],
    instruction=extract_sequences,
    inputs=[
        {
            "kwargs": {
                "hmmsearch_output_path": ("file_path", Path),
                "seqid2earliest_date": ("seqid2earliest_date", dict),
                "seqid2seq": ("filtered_seqid2seq", dict),
            },
        },
    ],
)


recipe.add(
    targets=[
        (("align_results", dict),),
    ],
    instruction=filter_and_align_template_chain_ids(
        date_cutoff="2021-09-30",
        day_diff_cutoff=60,
        min_seq_len=10,
        topk=20,
    ),
    inputs=[
        {
            "kwargs": {
                "query_seq": ("query_seq", str),
                "earliest_query_date": ("earliest_query_date", time.struct_time),
                "template_chain_ids": ("template_chain_ids", list[str]),
                "metadata_dict": ("template_metadata_map", dict),
            },
        },
    ],
)

recipe.add(
    targets=[
        (("template_mols", list),),
    ],
    instruction=load_templates,
    inputs=[
        {
            "kwargs": {
                "cif_db_path": ("cif_db_path", dict),
                "align_results": ("align_results", dict),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["template_mols"]
