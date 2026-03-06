from pathlib import Path

from datacooker import RecipeBook

from pipelines.instructions.msa_instructions import make_input_fasta, run_msa_search, run_signalp

"""Build a FASTA->MSA (SignalP + HHblits/HHfilter) Cooker."""

recipe = RecipeBook()

recipe.add(
	targets=[
		(
            ("input_fasta", Path),
            ("out_dir", Path),
		),
    ],
	instruction=make_input_fasta,
    inputs=[
        {
            "kwargs": {
                "seqid": ("seqid", str),
                "sequence": ("sequence", str),
                "output_dir": ("output_dir", Path),
            },
        },
    ],
)

recipe.add(
	targets=[
		(("trim_fasta", Path),),
    ],
	instruction=run_signalp,
    inputs=[
        {
            "kwargs": {
                "input_fasta": ("input_fasta", Path),
                "out_dir": ("out_dir", Path),
            },
        },
    ],
)

recipe.add(
    targets=[
        (("msa_results", dict),),
    ],
    instruction=run_msa_search,
    inputs=[
        {
            "kwargs": {
                "input_fasta": ("trim_fasta", Path),
                "out_dir": ("out_dir", Path),
                "cpu": ("cpu_per_job", int),
                "mem": ("mem_per_job", int),
                "db_ur30": ("db_ur30", str | Path),
                "db_bfd": ("db_bfd", str | Path),
                "hhsuite_bin_dir": ("hhsuite_bin_dir", str | Path),
            },
        },
    ],
)

RECIPE = recipe
TARGETS = ["msa_results"]
