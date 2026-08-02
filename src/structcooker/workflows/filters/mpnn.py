from datacooker import RecipeBook

from structcooker.instructions.transforms.pruning import select_features
from structcooker.mols import CIFMolAttached

"""Prune a CIFMol DB down to the features an MPNN reads.

One step, because pruning is one decision: which features survive. The list
lives in the db config's ``parameters`` block so the config, not this file,
documents what a given model-input DB contains.
"""

recipe = RecipeBook()


recipe.step(
    outputs=(("mpnn_dict", dict),),
    instruction=select_features,
    kwargs={
        "cifmol": ("cifmol", CIFMolAttached),
        "atom_features": ("atom_features", list),
        "residue_features": ("residue_features", list),
        "chain_features": ("chain_features", list),
        "atom_edge_features": ("atom_edge_features", list),
        "residue_edge_features": ("residue_edge_features", list),
        "chain_edge_features": ("chain_edge_features", list),
        "compact_index_table": ("compact_index_table", bool),
    },
)


RECIPE = recipe
TARGETS = ["mpnn_dict"]
