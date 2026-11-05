"""Number of templates stored for each seq_id of a template LMDB -> TSV seq_id, n_templates."""

import sys

import lmdb

from structcooker.instructions.transforms.codecs import from_bytes

env = lmdb.open(sys.argv[1], readonly=True, lock=False, readahead=False)
n = 0
with env.begin() as txn, open(sys.argv[2], "w") as out:
    for key, value in txn.cursor():
        out.write(f"{key.decode()}\t{len(from_bytes(value)['template_ids'])}\n")
        n += 1
print(n)
