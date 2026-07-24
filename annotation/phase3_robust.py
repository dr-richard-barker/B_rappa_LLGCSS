#!/usr/bin/env python3
"""
Phase 3 (hardened) — correlation-aware, multiple-testing-corrected gene-set test.

Upgrades over phase3_geneset_test.py:
  * effect measured from the 39-sample log-normalised expression matrix as a
    genotype×preservative-STRATIFIED dose contrast (controls the two confounds directly),
    instead of the marginal DESeq2 log2FC;
  * analytic competitive p = scipy Mann-Whitney U (exact-ish), plus
  * a CORRELATION-AWARE p from sample-label permutation stratified by genotype×preservative
    (shuffling samples preserves gene-gene correlation, which the MWU ignores and which
    otherwise inflates significance for co-expressed gene sets);
  * Benjamini-Hochberg FDR across the routes/tiers tested.

Writes annotation/phase3_dose_stats.tsv (per-gene stratified dose effect) for the figures,
and annotation/phase3_robust_results.tsv (the test table).
"""
import csv, os, re
import numpy as np
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
RNG = np.random.default_rng(20260724)
NPERM = 5000

def norm(g):
    m = re.fullmatch(r'(?i)bra(\d+)', g.strip().strip('"')); return f"Bra{m.group(1)}" if m else None

# ---- scent set ----
tier, routes = {}, {}
with open(os.path.join(HERE, "scent_geneset.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        tier[r["bra_id"]] = int(r["tier"]); routes[r["bra_id"]] = r["routes"].split(";")

# ---- expression matrix + sample design ----
path = os.path.join(ROOT, "Radiation model/Radi_model_with_interactions_deg_values_DESeq2.csv")
with open(path) as fh:
    rows = list(csv.reader(fh))
hdr = rows[0]; i_id = hdr.index("ensembl_ID")
start = hdr.index("Processed data:") + 1
samp_cols = [(i, hdr[i]) for i in range(start, len(hdr)) if hdr[i].strip()]
pat = re.compile(r'^(AL|W)(0|40)(D|R)_')
geno, dose, pres = [], [], []
keep_idx = []
for i, name in samp_cols:
    m = pat.match(name)
    if not m: continue
    keep_idx.append(i); geno.append(m.group(1)); dose.append(int(m.group(2))); pres.append(m.group(3))
geno = np.array(geno); dose = np.array(dose); pres = np.array(pres)
strata = np.array([g+p for g, p in zip(geno, pres)])   # AL D / AL R / W D / W R
print(f"samples: {len(keep_idx)}  strata: {dict(zip(*np.unique(strata, return_counts=True)))}")
print(f"  dose split: {dict(zip(*np.unique(dose, return_counts=True)))}")

genes, M = [], []
for r in rows[1:]:
    g = norm(r[i_id])
    if not g: continue
    try:
        vals = [float(r[i]) for i in keep_idx]
    except (ValueError, IndexError):
        continue
    genes.append(g); M.append(vals)
M = np.array(M); genes = np.array(genes)
print(f"genes with full matrix: {M.shape[0]} × {M.shape[1]}")

# ---- stratified dose effect per gene ----
def dose_effect(mat, dose_labels):
    """mean over strata of (mean_40 - mean_0) within stratum -> confound-adjusted dose effect."""
    eff = np.zeros(mat.shape[0])
    ns = 0
    for s in np.unique(strata):
        sel = strata == s
        d = dose_labels[sel]
        if (d == 40).sum() == 0 or (d == 0).sum() == 0: continue
        eff += mat[:, sel][:, d == 40].mean(1) - mat[:, sel][:, d == 0].mean(1)
        ns += 1
    return eff / ns
obs = dose_effect(M, dose)

# save per-gene stat for figures
with open(os.path.join(HERE, "phase3_dose_stats.tsv"), "w") as fo:
    fo.write("bra_id\tdose_effect\ttier\troutes\n")
    for g, e in zip(genes, obs):
        fo.write(f"{g}\t{e:.5f}\t{tier.get(g,'')}\t{';'.join(routes.get(g,[]))}\n")

# ---- permutation null (shuffle dose within each stratum) ----
gene_idx = {g: k for k, g in enumerate(genes)}
def idx_for(members):
    return np.array([gene_idx[g] for g in members if g in gene_idx])

def perm_null(set_idx):
    bg_mask = np.ones(M.shape[0], bool); bg_mask[set_idx] = False
    T = np.empty(NPERM)
    for b in range(NPERM):
        dp = dose.copy()
        for s in np.unique(strata):
            sel = np.where(strata == s)[0]
            dp[sel] = RNG.permutation(dp[sel])
        eff = dose_effect(M, dp)
        T[b] = eff[set_idx].mean() - eff[bg_mask].mean()
    return T

def test(name, members):
    idx = idx_for(members)
    if len(idx) < 5: return None
    bg = np.ones(M.shape[0], bool); bg[idx] = False
    s_eff, b_eff = obs[idx], obs[bg]
    Tobs = s_eff.mean() - b_eff.mean()
    # analytic competitive MWU (two-sided) on signed effect
    U, p_mwu = mannwhitneyu(s_eff, b_eff, alternative="two-sided")
    # correlation-aware permutation p
    Tnull = perm_null(idx)
    p_perm = (1 + np.sum(np.abs(Tnull) >= abs(Tobs))) / (NPERM + 1)
    return dict(name=name, n=len(idx), med_set=float(np.median(s_eff)),
                med_bg=float(np.median(b_eff)), Tobs=float(Tobs),
                p_mwu=float(p_mwu), p_perm=float(p_perm))

SETS = [("Tier1_core", [g for g in tier if tier[g] == 1]),
        ("Tier1+2_all", list(tier))]
for rt in ["terpenoid", "benzenoid_phenylpropanoid", "ester", "fatty_acid_GLV", "apocarotenoid"]:
    SETS.append((f"route:{rt}", [g for g in tier if rt in routes.get(g, [])]))

res = [r for r in (test(n, m) for n, m in SETS) if r]
# BH across all tests (use permutation p)
ps = np.array([r["p_perm"] for r in res]); order = np.argsort(ps); m = len(ps)
q = np.empty(m); prev = 1.0
for rank, i in enumerate(order[::-1]):
    prev = min(prev, ps[i] * m / (m - rank)); q[i] = prev
for r, qi in zip(res, q): r["q_perm_BH"] = float(qi)

with open(os.path.join(HERE, "phase3_robust_results.tsv"), "w") as fo:
    cols = ["name", "n", "med_set", "med_bg", "Tobs", "p_mwu", "p_perm", "q_perm_BH"]
    fo.write("\t".join(cols) + "\n")
    for r in res: fo.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

print(f"\n{'set':28s} {'n':>4s} {'med_set':>8s} {'med_bg':>7s} {'T':>7s} {'p_MWU':>8s} {'p_perm':>8s} {'q_BH':>7s}")
for r in res:
    print(f"{r['name']:28s} {r['n']:4d} {r['med_set']:+8.3f} {r['med_bg']:+7.3f} {r['Tobs']:+7.3f} "
          f"{r['p_mwu']:8.3g} {r['p_perm']:8.3g} {r['q_perm_BH']:7.3g}")
print(f"\n(permutation: {NPERM} sample-label shuffles stratified by genotype×preservative; "
      "q = BH-FDR across the {} tests)".format(len(res)))
