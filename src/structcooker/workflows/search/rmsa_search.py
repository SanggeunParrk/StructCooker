from pathlib import Path

from datacooker import RecipeBook

from structcooker.instructions.transforms.msa import (
    make_input_fasta,
    run_rna_msa_search,
)

"""Build a FASTA->RNA MSA (AlphaFold3-style nhmmer + hmmalign) Cooker."""

recipe = RecipeBook()

recipe.step(
    outputs=(
        ("input_fasta", Path),
        ("out_dir", Path),
    ),
    instruction=make_input_fasta,
    kwargs={
        "seqid": ("seqid", str),
        "sequence": ("sequence", str),
        "output_dir": ("output_dir", Path),
    },
)

recipe.step(
    outputs=(("msa_results", str),),
    instruction=run_rna_msa_search,
    kwargs={
        "input_fasta": ("input_fasta", Path),
        "out_dir": ("out_dir", Path),
        "cpu": ("cpu_per_job", int),
        "db_rfam": ("db_rfam", str | Path),
        "db_rnacentral": ("db_rnacentral", str | Path),
        "db_nt": ("db_nt", str | Path),
    },
)

RECIPE = recipe
TARGETS = ["msa_results"]
