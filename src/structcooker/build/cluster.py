"""Settings that depend on the cluster rather than on the data, all from the environment.

The node size uses DataCooker's own variables, which its planner reads too
(``DATACOOKER_CPUS_PER_NODE``, ``DATACOOKER_NODE_MEM_GB``); the defaults are the CSSB cluster's
(112-core, ~500 GB nodes). The rest lives in DataCooker's executor: ``DATACOOKER_PARTITION``,
``DATACOOKER_QOS`` and ``DATACOOKER_EXCLUDE_NODES`` (nodes to keep jobs off). The CSSB values
are collected in scripts/cluster/cssb.env.
"""
from __future__ import annotations

import os

# A full node, as one job asks for it (parallel searches, materializers fanning out, verification).
NODE_CORES = int(os.environ.get("DATACOOKER_CPUS_PER_NODE", "112"))
NODE_MEM_GB = int(os.environ.get("DATACOOKER_NODE_MEM_GB", "490"))
