# Release transfers (2026-10-28/29 refresh)

Copies of the scripts that moved the 40%-cluster refresh to H100 (`psk6950@147.47.135.34:7777`, `~/data/BioMolDB_20261014`)
and cssb (`psk6950@147.47.218.30:7994`, `/public_data02/BioMolDB_20261014`). They were run from `logs/transfer/` (git-ignored),
which also holds the generated `units_<target>_rank<r>.txt`, `files_<target>.tsv` and `prefixes_<target>.txt`; paths inside the
scripts point there.

* `plan2.py h100|cssb|h100x|h100c|cssbc K` -- file lists (ROOT, relative path, destination prefix) balanced over K ranks, group 1 first.
* `transfer2.sbatch` + `send_one.sh` -- rsync a rank's list, 3-4 streams, 3 tries per file; rerun resumes.
* `verify2.sh <target>` -- every planned file on the remote with the same size.
* `chunk_one.sh`, `chunk_afm.sbatch`, `finalize_afm.sh` -- a file too big for one stream (1 TB) written in 4 GiB ranges in place,
  then size plus 24 random 64 MiB md5 regions.
* `send_complex.sbatch` -- after the PDB complex design DB is validated.
* `send_foldbench.sbatch`, `foldbench_native.sbatch`, `send_meta_h100.sbatch` -- the FoldBench DB and the extra metadata.

Lessons: a `shards.json` may point outside its collection (the long `cif_attached` did) -- follow it before planning; `ssh` inside a
`while read` loop needs `-n`; relative remote paths break after `cd`.
