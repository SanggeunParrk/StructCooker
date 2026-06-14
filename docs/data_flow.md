# Data Processing Flow

This note captures the original end-to-end data processing flow that was previously stored in `pipelines/README.md`.

## Main Flow

1. Download CCD information from wwPDB.
2. Download CIF files from RCSB.
3. Build the CCD LMDB.
4. Build the CIF LMDB.
   Note: some CIF files are broken and may need to be replaced with fixed versions.
5. Extract FASTA files.
6. Build the sequence ID map.
   If an older sequence ID map already exists, it can be reused.
7. Download `sabdab_summary_all.tsv`.
   Note: SAbDab is updated weekly.
8. Run sequence clustering.
9. Run MSA search.
   This requires the sequence databases to be available.
10. Run template search.
   This requires the template database to be built first.
11. Build the MSA database.
   Key: sequence ID.
12. Build the template database.
   Key: sequence ID.

## Optional Flow

1. Filter the dataset.
2. Extract interacting sequence IDs and interacting sequence clusters.
3. Extract metadata.
4. Create the train/valid split.
5. Run valid filtering.
