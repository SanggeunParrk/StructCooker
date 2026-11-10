"""Recipe: interface manifest + seq_ids + 40% clusters (``manifest_full.tsv``) and ``report.md``."""

from datacooker import RecipeBook

from structcooker.instructions.interface.tables import write_manifest_full

recipe = RecipeBook()
recipe.step(
    outputs=(("status", str),),
    instruction=write_manifest_full,
    kwargs={
        "manifest_path": ("manifest_path", str),
        "fasta_path": ("fasta_path", str),
        "seq_id_map_path": ("seq_id_map_path", str),
        "cluster_path": ("cluster_path", str),
        "cif_metadata_path": ("cif_metadata_path", str),
        "out_path": ("out_path", str),
        "out_report": ("out_report", str),
    },
)
RECIPE = recipe
TARGETS = ["status"]
