"""Recipe: ANARCI chain type of every unique PDB protein sequence (the interface DB's antibody cache)."""

from datacooker import RecipeBook

from structcooker.instructions.interface.antibody import annotate_antibodies

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=annotate_antibodies,
    kwargs={
        "fasta_path": ("fasta_path", str),
        "out_path": ("out_path", str),
        "n_jobs": ("n_jobs", int),
    },
)
RECIPE = recipe
TARGETS = ["status"]
