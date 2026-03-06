1. download ccd infos from wwpdb
2. download cif files from rcsb
3. build ccd lmdb.
4. cif lmdb (Note that there are broken cif files and corresponding fixed file so you have to replace it.)
5. extract fasta files
6. build seq id map (if there exists old seq id map you can use it.)
7. download sabdab_summary.tsv
8. sequence cluster
9. msa search (You have to download seq db)
10. template search (You have to build template db)
11. build msa db (key = seq ID)
12. build template db (key = seq ID)

==== (Optional) ====
1. filter out
2. extract interacting seq ids & interacting seq clusters
3. extract metadata
4. train valid split
5. valid filter