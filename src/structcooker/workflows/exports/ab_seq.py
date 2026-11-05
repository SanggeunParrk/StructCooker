"""Recipe: the antibody records of a PDB fasta (``ab_seq.fasta``), headers as in the input."""

from datacooker import RecipeBook

from structcooker.instructions.transforms.antibody_fasta import filter_fasta_antibody

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=filter_fasta_antibody,
    kwargs={
        "fasta_path": ("fasta_path", str),
        "sabdab_summary_path": ("sabdab_summary_path", str),
        "out_path": ("out_path", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
