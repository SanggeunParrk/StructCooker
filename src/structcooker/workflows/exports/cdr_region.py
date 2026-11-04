"""Recipe: the CDR residues of the antibody chains of a fasta (``cdr_region.fasta``)."""
from datacooker import RecipeBook

from structcooker.instructions.transforms.antibody_fasta import extract_cdr_regions

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=extract_cdr_regions,
    kwargs={
        "fasta_path": ("fasta_path", str),
        "out_path": ("out_path", str),
        "cdrs": ("cdrs", list),
        "scheme": ("scheme", str),
        "batch_size": ("batch_size", int),
        "ncpu": ("ncpu", int),
        "failed_path": ("failed_path", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
