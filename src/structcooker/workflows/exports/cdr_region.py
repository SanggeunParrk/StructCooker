"""Recipe: the CDR positions of the antibody sequences of a fasta (``cdr_region.csv``)."""

from datacooker import RecipeBook

from structcooker.instructions.transforms.antibody_fasta import extract_cdr_table

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=extract_cdr_table,
    kwargs={
        "fasta_path": ("fasta_path", str),
        "seq_id_map_path": ("seq_id_map_path", str),
        "out_path": ("out_path", str),
        "scheme": ("scheme", str),
        "batch_size": ("batch_size", int),
        "ncpu": ("ncpu", int),
        "failed_path": ("out_failed", str),
        "extra_path": ("out_extra", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
