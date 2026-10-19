"""AF3-style vs MiniWorld-style sampling over the same PDB pool, plus the distillation sources.

Input: the tables ``extract_items.py`` wrote (chains.*.tsv, interfaces.*.tsv). Output: a markdown
report of exact per-item sampling probabilities, summarised. Nothing is simulated by drawing;
every number is the expectation under the sampler.

AF3 (supplementary information 2.5.1, eq. 1): an item is a chain or an interface (two chains whose
minimum heavy-atom distance is < 5 A), sampled with
        w  ~  beta_r / N_clust * (alpha_prot * n_prot + alpha_nuc * n_nuc + alpha_ligand * n_ligand)
beta_chain 0.5, beta_interface 1, alpha_prot 3, alpha_nuc 3, alpha_ligand 1; N_clust is the number
of chains (interfaces) in the item's cluster; an interface cluster is the pair of chain clusters.
Chains with fewer than 4 resolved residues are dropped (2.5.4).

MiniWorld (src/miniworld/data/dataloader/sources.py, ``read_train_items``): the pool is the
legacy edge_node file -- every valid chain whose cluster never occurs in an interface ("sole"),
and every stored contact edge (any atom pair within 6 A). Weight of a row:
        W[type] / (clusters of that type) / (rows of the cluster)
with the fixed per-type weights of the training config (``sampler:``), type from the edge id.

    python compare_sampling.py --items <dir> --out report.md [--cluster-tsv ... --seqid-tsv ...]
"""
from __future__ import annotations

import argparse
import glob
import re
from array import array
from collections import defaultdict
from pathlib import Path

import numpy as np

# ---- MiniWorld, copied from src/miniworld/data/dataloader/sources.py and the training config ----
MW_WEIGHTS = {
    "protein_protein": 25.0, "protein_ligand": 25.0, "protein_dna": 10.0, "protein_rna": 10.0,
    "antibody_protein": 15.0, "dna_dna": 5.0, "rna_rna": 5.0, "dna_rna": 0.5,
    "antibody_antibody": 0.5, "antibody_ligand": 1.0, "na_ligand": 1.0, "etc_interface": 1.0, "sole": 1.0,
}


def mw_edge_type_from_letters(parse: set[str]) -> str:
    """``_pdb_edge_type`` with the type letters given directly (it read them out of the edge id)."""
    if parse == {"P"}:
        return "protein_protein"
    if parse <= {"P", "L", "B"} and "P" in parse:
        return "protein_ligand"
    if parse == {"P", "D"}:
        return "protein_dna"
    if parse == {"P", "R"}:
        return "protein_rna"
    if parse == {"P", "A"}:
        return "antibody_protein"
    if parse == {"D"}:
        return "dna_dna"
    if parse == {"R"}:
        return "rna_rna"
    if parse == {"D", "R"}:
        return "dna_rna"
    if parse == {"A"}:
        return "antibody_antibody"
    if parse <= {"A", "L", "B"} and "A" in parse:
        return "antibody_ligand"
    if parse <= {"N", "L", "B"} and "N" in parse:
        return "na_ligand"
    return "etc_interface"


def mw_edge_type_verbatim(edge_id: str) -> str:
    """``_pdb_edge_type`` exactly as MiniWorld runs it (regex on the edge id string)."""
    if "_" not in edge_id:
        return "sole"
    return mw_edge_type_from_letters(set(re.findall(r"c([A-Z])", edge_id)))


# ---- AF3 ----
BETA = {"chain": 0.5, "interface": 1.0}
ALPHA = {"prot": 3.0, "nuc": 3.0, "lig": 1.0}
POLYMER = set("PAQRDN")
KIND = {**dict.fromkeys("PAQ", "prot"), **dict.fromkeys("RDN", "nuc")}     # everything else: ligand


def kind_of(letter: str) -> str:
    return KIND.get(letter, "lig")


def letter_of(cluster: str) -> str:
    """Type letter of a cluster id ``c<DB>_<letter><number>``."""
    return cluster.split("_", 1)[1][0]


# ---------------------------------------------------------------------------------------------
def load_pool(items_dir: Path) -> dict:
    """Chain rows and interface rows as numpy arrays (cluster and entry strings become int codes)."""
    entries: dict[str, int] = {}
    clusters: dict[str, int] = {}
    names: list[str] = []
    cl_letter: list[str] = []

    def ecode(s: str) -> int:
        v = entries.get(s)
        if v is None:
            v = entries[s] = len(entries)
        return v

    def ccode(s: str) -> int:
        v = clusters.get(s)
        if v is None:
            v = clusters[s] = len(clusters)
            names.append(s)
            cl_letter.append(letter_of(s))
        return v

    seqs: dict[str, int] = {}
    seq_names: list[str] = []

    def scode(x: str) -> int:
        v = seqs.get(x)
        if v is None:
            v = seqs[x] = len(seqs)
            seq_names.append(x)
        return v

    ch_entry, ch_cl, ch_ok4, ch_seq = array("i"), array("i"), array("b"), array("i")
    dropped: dict[str, set[str]] = defaultdict(set)         # "entry|asm" -> chains with < 4 resolved residues
    for path in sorted(glob.glob(str(items_dir / "chains.*.tsv"))):
        with open(path) as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                entry, asm, chain, _etype, cluster, nres, _atoms = parts[:7]
                if len(parts) > 7:
                    ch_seq.append(scode(parts[7]))
                c = ccode(cluster)
                ok4 = 1
                if cl_letter[c] in POLYMER and int(nres) < 4:
                    ok4 = 0
                    dropped[f"{entry}|{asm}"].add(chain)
                ch_entry.append(ecode(entry))
                ch_cl.append(c)
                ch_ok4.append(ok4)

    if_entry, if_c1, if_c2, if_in5, if_ok4 = array("i"), array("i"), array("i"), array("b"), array("b")
    if_s1, if_s2 = array("i"), array("i")
    for path in sorted(glob.glob(str(items_dir / "interfaces.*.tsv"))):
        with open(path) as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                entry, asm, ch1, ch2, cl1, cl2, in5 = parts[:7]
                if len(parts) > 8:
                    if_s1.append(scode(parts[7]))
                    if_s2.append(scode(parts[8]))
                ok = 1
                bad = dropped.get(f"{entry}|{asm}")
                if bad and (ch1 in bad or ch2 in bad):
                    ok = 0
                if_entry.append(ecode(entry))
                if_c1.append(ccode(cl1))
                if_c2.append(ccode(cl2))
                if_in5.append(int(in5))
                if_ok4.append(ok)

    def arr(a, dt):
        return np.frombuffer(a, dtype=dt) if len(a) else np.zeros(0, dtype=dt)

    return {
        "n_entries": len(entries), "n_clusters": len(clusters), "letters": np.array(cl_letter), "names": names,
        "ch_entry": arr(ch_entry, np.int32), "ch_cl": arr(ch_cl, np.int32), "ch_ok4": arr(ch_ok4, np.int8).astype(bool),
        "if_entry": arr(if_entry, np.int32), "if_c1": arr(if_c1, np.int32), "if_c2": arr(if_c2, np.int32),
        "if_in5": arr(if_in5, np.int8).astype(bool), "if_ok4": arr(if_ok4, np.int8).astype(bool),
        "ch_seq": arr(ch_seq, np.int32), "if_s1": arr(if_s1, np.int32), "if_s2": arr(if_s2, np.int32), "seq_names": seq_names,
    }


def edge_codes(a: np.ndarray, b: np.ndarray, k: int) -> np.ndarray:
    """One int64 per item cluster: the sorted pair of chain clusters, or the chain cluster alone (b = -1)."""
    lo = np.where(b >= 0, np.minimum(a, b), a).astype(np.int64)
    hi = np.where(b >= 0, np.maximum(a, b), -1).astype(np.int64)
    return lo * (k + 1) + (hi + 1)


def letter_codes(letters: np.ndarray) -> tuple[np.ndarray, list[str]]:
    alphabet = sorted(set(letters.tolist()))
    table = {ch: i for i, ch in enumerate(alphabet)}
    return np.array([table[x] for x in letters.tolist()], dtype=np.int64), alphabet


def category(lcode: np.ndarray, alphabet: list[str], c1: np.ndarray, c2: np.ndarray, is_iface: np.ndarray) -> np.ndarray:
    """Common label of an item (chain:<kind> or the MiniWorld interface bucket); vectorised via unique combos."""
    a = lcode[c1]
    b = np.where(c2 >= 0, lcode[np.maximum(c2, 0)], -1)
    key = (a * 64 + (b + 1)) * 2 + is_iface.astype(np.int64)
    uniq, inv = np.unique(key, return_inverse=True)
    labels = []
    for u in uniq.tolist():
        f = u % 2
        rest = u // 2
        la, lb = rest // 64, rest % 64 - 1
        x, y = alphabet[la], (alphabet[lb] if lb >= 0 else "")
        labels.append(mw_edge_type_from_letters({x, y}) if f else f"chain:{kind_of(x)}")
    return np.array(labels, dtype=object)[inv]


SIZE_BINS = [(1, 1), (2, 3), (4, 10), (11, 100), (101, 1000), (1001, 10**12)]


def bin_name(lo: int, hi: int) -> str:
    return f"{lo}-{hi}" if hi < 10**12 else f">={lo}"


def summarize(p, cat, size, c1, c2, entry, n_clusters, n_entries, ns=(128_000, 1_000_000, 10_000_000)) -> dict:
    res: dict = {}
    cats, cinv = np.unique(cat, return_inverse=True)
    mass = np.bincount(cinv, weights=p, minlength=len(cats))
    res["mass_by_category"] = {c: float(m) for c, m in zip(cats.tolist(), mass.tolist(), strict=True)}
    res["mass_by_size"] = {bin_name(lo, hi): float(p[(size >= lo) & (size <= hi)].sum()) for lo, hi in SIZE_BINS}
    m_cl = np.bincount(c1, weights=p, minlength=n_clusters)
    both = (c2 >= 0) & (c2 != c1)
    m_cl += np.bincount(c2[both], weights=p[both], minlength=n_clusters)
    m_en = np.bincount(entry, weights=p, minlength=n_entries)
    res["expected_unique"] = {}
    for n in ns:
        mc, me = m_cl[m_cl > 0], m_en[m_en > 0]
        res["expected_unique"][n] = (
            float(np.sum(-np.expm1(n * np.log1p(-np.minimum(mc, 1 - 1e-12))))),
            float(np.sum(-np.expm1(n * np.log1p(-np.minimum(me, 1 - 1e-12))))),
        )
    q = m_cl[m_cl > 0] / m_cl.sum()
    res["eff_clusters"] = float(np.exp(-(q * np.log(q)).sum()))
    r = m_en[m_en > 0] / m_en.sum()
    res["eff_entries"] = float(np.exp(-(r * np.log(r)).sum()))
    res["top100_cluster_mass"] = float(np.sort(m_cl)[::-1][:100].sum() / m_cl.sum())
    res["max_entry_prob"] = float(m_en.max())
    res["n_clusters_touched"] = int((m_cl > 0).sum())
    return res


def pdb_section(pool: dict, out: list[str]) -> dict:
    k = pool["n_clusters"]
    letters = pool["letters"]
    lcode, alphabet = letter_codes(letters)
    kinds = np.array([kind_of(x) for x in letters.tolist()])
    n_ch, n_if = len(pool["ch_cl"]), len(pool["if_c1"])
    out.append(f"Pool read: {n_ch:,} chain rows, {n_if:,} stored contact (6 A) interface rows, "
               f"{pool['n_entries']:,} entries, {k:,} chain clusters. Type letters seen: {', '.join(alphabet)}.\n")

    # ----- MiniWorld pool: sole chains + every stored contact edge
    iface_clusters = np.zeros(k, dtype=bool)
    iface_clusters[pool["if_c1"]] = True
    iface_clusters[pool["if_c2"]] = True
    sole = ~iface_clusters[pool["ch_cl"]]
    mw_c1 = np.concatenate([pool["ch_cl"][sole], pool["if_c1"]])
    mw_c2 = np.concatenate([np.full(int(sole.sum()), -1, np.int32), pool["if_c2"]])
    mw_entry = np.concatenate([pool["ch_entry"][sole], pool["if_entry"]])
    mw_is_if = np.concatenate([np.zeros(int(sole.sum()), bool), np.ones(n_if, bool)])
    mw_edge = edge_codes(mw_c1, mw_c2, k)
    uniq, first, inv, cnt = np.unique(mw_edge, return_index=True, return_inverse=True, return_counts=True)
    btype = np.empty(len(uniq), dtype=object)
    pair_cache: dict[tuple[str, str], str] = {}
    for j, row in enumerate(first.tolist()):
        if mw_c2[row] < 0:
            btype[j] = "sole"
            continue
        key = (letters[mw_c1[row]], letters[mw_c2[row]])
        v = pair_cache.get(key)
        if v is None:
            v = pair_cache[key] = mw_edge_type_from_letters(set(key))
        btype[j] = v
    n_edges_in_type: dict[str, int] = defaultdict(int)
    for t in btype.tolist():
        n_edges_in_type[t] += 1
    w_edge = np.array([MW_WEIGHTS[t] / n_edges_in_type[t] for t in btype.tolist()])
    mw_w = (w_edge / cnt)[inv]
    mw_p = mw_w / mw_w.sum()
    mw_size = cnt[inv]

    # ----- AF3 pool: every chain with >= 4 resolved residues, every 5 A interface
    ch_keep = pool["ch_ok4"]
    if_keep = pool["if_in5"] & pool["if_ok4"]
    af_c1 = np.concatenate([pool["ch_cl"][ch_keep], pool["if_c1"][if_keep]])
    af_c2 = np.concatenate([np.full(int(ch_keep.sum()), -1, np.int32), pool["if_c2"][if_keep]])
    af_entry = np.concatenate([pool["ch_entry"][ch_keep], pool["if_entry"][if_keep]])
    af_is_if = np.concatenate([np.zeros(int(ch_keep.sum()), bool), np.ones(int(if_keep.sum()), bool)])
    af_key = edge_codes(af_c1, af_c2, k) * 2 + af_is_if.astype(np.int64)      # item kind is part of the cluster
    uniq2, first2, inv2, cnt2 = np.unique(af_key, return_index=True, return_inverse=True, return_counts=True)
    n_prot = (kinds[af_c1] == "prot").astype(float) + ((af_c2 >= 0) & (kinds[np.maximum(af_c2, 0)] == "prot"))
    n_nuc = (kinds[af_c1] == "nuc").astype(float) + ((af_c2 >= 0) & (kinds[np.maximum(af_c2, 0)] == "nuc"))
    n_lig = (kinds[af_c1] == "lig").astype(float) + ((af_c2 >= 0) & (kinds[np.maximum(af_c2, 0)] == "lig"))
    beta = np.where(af_is_if, BETA["interface"], BETA["chain"])
    af_w = beta / cnt2[inv2] * (ALPHA["prot"] * n_prot + ALPHA["nuc"] * n_nuc + ALPHA["lig"] * n_lig)
    af_p = af_w / af_w.sum()
    af_size = cnt2[inv2]

    out.append("### Pools\n")
    out.append("| | MiniWorld pool | AF3 pool |\n|---|---:|---:|")
    out.append(f"| rows | {len(mw_p):,} | {len(af_p):,} |")
    out.append(f"| chain rows | {int((~mw_is_if).sum()):,} (clusters never in an interface) | {int((~af_is_if).sum()):,} (every chain) |")
    out.append(f"| interface rows | {int(mw_is_if.sum()):,} (6 A, any atom) | {int(af_is_if.sum()):,} (5 A) |")
    out.append(f"| distinct item clusters | {len(uniq):,} | {len(uniq2):,} |\n")

    mw_cat = category(lcode, alphabet, mw_c1, mw_c2, mw_is_if)
    af_cat = category(lcode, alphabet, af_c1, af_c2, af_is_if)
    af_cat_counts: dict[str, int] = defaultdict(int)
    for c in af_cat[first2].tolist():
        af_cat_counts[c] += 1
    af_cat_rows = {c: int((af_cat == c).sum()) for c in af_cat_counts}
    mw_cat_rows = {c: int((mw_cat == c).sum()) for c in set(mw_cat.tolist())}
    s_mw = summarize(mw_p, mw_cat, mw_size, mw_c1, mw_c2, mw_entry, k, pool["n_entries"])
    s_af = summarize(af_p, af_cat, af_size, af_c1, af_c2, af_entry, k, pool["n_entries"])

    out.append("### Probability of an item class (within the PDB source)\n")
    out.append("| item class | MiniWorld | AF3 | MiniWorld / AF3 |\n|---|---:|---:|---:|")
    # MiniWorld's own `sole` rows are single chains: show them under chain:<kind> for comparison
    mw_cat_chain = np.where(~mw_is_if, np.char.add("chain:", kinds[mw_c1]).astype(object), mw_cat)
    s_mw2 = summarize(mw_p, mw_cat_chain, mw_size, mw_c1, mw_c2, mw_entry, k, pool["n_entries"])
    cats = sorted(set(s_mw2["mass_by_category"]) | set(s_af["mass_by_category"]),
                  key=lambda c: -max(s_mw2["mass_by_category"].get(c, 0), s_af["mass_by_category"].get(c, 0)))
    for c in cats:
        a, b = s_mw2["mass_by_category"].get(c, 0.0), s_af["mass_by_category"].get(c, 0.0)
        out.append(f"| {c} | {a:.2%} | {b:.2%} | {f'{a / b:.2f}' if b > 0 else '—'} |")
    chain_mw = sum(v for c, v in s_mw2["mass_by_category"].items() if c.startswith("chain:"))
    chain_af = sum(v for c, v in s_af["mass_by_category"].items() if c.startswith("chain:"))
    out.append(f"\nSingle-chain items together: MiniWorld {chain_mw:.2%}, AF3 {chain_af:.2%}.\n")

    out.append("### Probability by cluster size (rows in the item's own cluster, in that sampler's pool)\n")
    out.append("| cluster size | MiniWorld | AF3 |\n|---|---:|---:|")
    for key in s_mw["mass_by_size"]:
        out.append(f"| {key} | {s_mw['mass_by_size'][key]:.2%} | {s_af['mass_by_size'][key]:.2%} |")
    out.append("")

    out.append("### Coverage and concentration\n")
    out.append("| | MiniWorld | AF3 |\n|---|---:|---:|")
    out.append(f"| chain clusters reachable | {s_mw['n_clusters_touched']:,} | {s_af['n_clusters_touched']:,} |")
    out.append(f"| effective number of chain clusters (exp entropy of the touch probabilities) | {s_mw['eff_clusters']:,.0f} | {s_af['eff_clusters']:,.0f} |")
    out.append(f"| effective number of entries | {s_mw['eff_entries']:,.0f} | {s_af['eff_entries']:,.0f} |")
    out.append(f"| touch probability of the 100 most-touched chain clusters | {s_mw['top100_cluster_mass']:.2%} | {s_af['top100_cluster_mass']:.2%} |")
    out.append(f"| largest single-entry probability | {s_mw['max_entry_prob']:.2e} | {s_af['max_entry_prob']:.2e} |")
    for n in (128_000, 1_000_000, 10_000_000):
        a, b = s_mw["expected_unique"][n], s_af["expected_unique"][n]
        out.append(f"| expected distinct chain clusters touched in {n:,} draws | {a[0]:,.0f} | {b[0]:,.0f} |")
        out.append(f"| expected distinct entries drawn in {n:,} draws | {a[1]:,.0f} | {b[1]:,.0f} |")
    out.append("")

    def contains(p, c1, c2):
        k1 = kinds[c1]
        k2 = np.where(c2 >= 0, kinds[np.maximum(c2, 0)], "")
        has_nuc = (k1 == "nuc") | (k2 == "nuc")
        has_lig = (k1 == "lig") | (k2 == "lig")
        return float(p[has_nuc].sum()), float(p[has_lig].sum()), float(p[~has_nuc & ~has_lig].sum())
    a, b = contains(mw_p, mw_c1, mw_c2), contains(af_p, af_c1, af_c2)
    out.append("### What a draw contains\n")
    out.append("| | MiniWorld | AF3 |\n|---|---:|---:|")
    out.append(f"| contains a nucleic-acid chain | {a[0]:.2%} | {b[0]:.2%} |")
    out.append(f"| contains a ligand | {a[1]:.2%} | {b[1]:.2%} |")
    out.append(f"| protein only | {a[2]:.2%} | {b[2]:.2%} |\n")

    # ----- MiniWorld's type regex applied to the new cluster-id format
    names = pool["names"]
    counts: dict[str, int] = defaultdict(int)
    rows = np.nonzero(mw_is_if)[0][:200000].tolist()
    for row in rows:
        counts[mw_edge_type_verbatim(f"{names[mw_c1[row]]}_{names[mw_c2[row]]}")] += 1
    sole_ids = {mw_edge_type_verbatim(names[mw_c1[row]]) for row in np.nonzero(~mw_is_if)[0][:1000].tolist()}
    out.append("### MiniWorld's `_pdb_edge_type`, unchanged, on the new cluster ids (`cPDB_<letter><number>`)\n")
    out.append(f"{len(rows):,} interface rows classified: "
               + ", ".join(f"{t} {n:,}" for t, n in sorted(counts.items(), key=lambda x: -x[1]))
               + f". Single-chain ids classified as: {sorted(sole_ids)}.\n")
    return {
        "mw_mass": s_mw2["mass_by_category"], "af_mass": s_af["mass_by_category"],
        "mw_clusters": dict(n_edges_in_type), "af_clusters": dict(af_cat_counts),
        "mw_rows": mw_cat_rows, "af_rows": af_cat_rows,
        "mw_pool_rows": len(mw_p), "af_pool_rows": len(af_p),
        "mw_eff_entries": s_mw["eff_entries"], "af_eff_entries": s_af["eff_entries"],
        "mw_eff_clusters": s_mw["eff_clusters"], "af_eff_clusters": s_af["eff_clusters"],
        "n_chain_clusters": int(len(set(np.concatenate([af_c1, af_c2[af_c2 >= 0]]).tolist()))),
    }


def distillation_section(cluster_tsv: Path, seqid_tsv: Path, out: list[str]) -> None:
    """Within-source weighting of the distillation sets: AF3 uniform vs MiniWorld cluster-uniform x row-uniform."""
    cl_of: dict[str, int] = {}
    with open(cluster_tsv) as handle:
        for i, line in enumerate(handle):
            _rep, members = line.rstrip("\n").split("\t")
            for m in members.split(","):
                cl_of[m] = i
    codes: dict[str, array] = {"P": array("i"), "R": array("i")}
    missing = 0
    with open(seqid_tsv) as handle:
        for line in handle:
            _entry, sid = line.rstrip("\n").split("\t")
            if sid[0] not in codes:
                continue
            c = cl_of.get(sid)
            if c is None:
                missing += 1
                continue
            codes[sid[0]].append(c)
    out.append("## Distillation sources: within-source weighting\n")
    out.append("AF3 draws a distillation example uniformly (SI table 3, 'uniform'). MiniWorld gives every cluster the same "
               "mass and splits it evenly over the cluster's entries (`1 / rowcount`, then source-normalised).\n")
    out.append(f"Entries whose sequence is missing from the cluster table: {missing:,}.\n")
    for kind, label in (("P", "protein monomers (long and short sets together)"), ("R", "RNA monomers")):
        arr = np.frombuffer(codes[kind], dtype=np.int32)
        if len(arr) == 0:
            continue
        _u, inv, cnt = np.unique(arr, return_inverse=True, return_counts=True)
        n, ncl = len(arr), len(cnt)
        size = cnt[inv]
        p_af = np.full(n, 1.0 / n)
        p_mw = (1.0 / ncl) / size
        out.append(f"### {label}\n")
        out.append(f"{n:,} entries in {ncl:,} clusters at 30% identity (mean {n / ncl:.2f} entries per cluster, "
                   f"largest {int(cnt.max()):,}).\n")
        out.append("| cluster size | share of entries | AF3 uniform | MiniWorld cluster-uniform | MiniWorld / AF3 |\n|---|---:|---:|---:|---:|")
        for lo, hi in SIZE_BINS:
            m = (size >= lo) & (size <= hi)
            a, b = float(p_af[m].sum()), float(p_mw[m].sum())
            out.append(f"| {bin_name(lo, hi)} | {m.mean():.2%} | {a:.2%} | {b:.2%} | {b / a if a else 0:.2f} |")
        eff_mw = float(np.exp(-(p_mw * np.log(p_mw)).sum()))
        out.append(f"\nEffective number of entries: AF3 {n:,} (uniform), MiniWorld {eff_mw:,.0f}.\n")


def mixture_section(out: list[str]) -> None:
    af3 = {"weighted PDB": 0.5, "disordered PDB": 0.02, "protein monomer (long)": 0.495,
           "short protein monomer": 0.005, "RNA distillation": 0.05}
    mw = {"weighted PDB": 0.5, "disordered PDB": 0.0, "protein monomer (long)": 0.495,
          "short protein monomer": 0.005, "RNA distillation": 0.0}
    sa, sm = sum(af3.values()), sum(mw.values())
    out.append("## Dataset mixture\n")
    out.append("AF3 SI table 3 (pre-training rows; the transcription-factor sets are fine-tuning only) against the "
               "`source_weights` of `H100_d512.yaml` (RNA distillation not declared, disordered PDB masked out of the "
               "sampler). Both are normalised over the active sources.\n")
    out.append("| source | AF3 weight | AF3 share | MiniWorld weight | MiniWorld share |\n|---|---:|---:|---:|---:|")
    for name in af3:
        out.append(f"| {name} | {af3[name]} | {af3[name] / sa:.2%} | {mw[name]} | {mw[name] / sm:.2%} |")
    out.append(f"| sum | {sa:.3f} | | {sm:.3f} | |\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--items", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--cluster-tsv", type=Path)
    ap.add_argument("--seqid-tsv", type=Path)
    args = ap.parse_args()
    out: list[str] = ["# Sampling: AF3 eq. 1 vs MiniWorld\n"]
    mixture_section(out)
    out.append("## PDB source\n")
    pdb_section(load_pool(args.items), out)
    if args.cluster_tsv and args.seqid_tsv:
        distillation_section(args.cluster_tsv, args.seqid_tsv, out)
    args.out.write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
