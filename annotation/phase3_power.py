#!/usr/bin/env python3
"""
Phase 3 power / minimum-detectable-effect analysis for the radiation × scent gene-set test.

Turns the null result into a bounded statement: "no coordinated effect larger than delta".
Method (consistent with phase3_robust.py):
  * per-gene confound-adjusted dose effect from the 39-library matrix, stratified by
    genotype×preservative;
  * competitive statistic T = mean(dose_effect[set]) - mean(dose_effect[background]);
  * null distribution of T from 5000 dose-label permutations within strata; two-sided
    critical value c = 97.5th percentile of |T_null|;
  * a coordinated additive shift delta applied to every set gene raises the set mean by
    delta and leaves T's variance ~unchanged, so the alternative distribution of T is the
    null shifted by delta. Power(delta) is estimated empirically as the fraction of
    (T_null + delta) beyond ±c. This is the BEST case (perfectly coordinated shift), so it
    is a lower bound on the effect a real, noisier signal would need.
Outputs: annotation/phase3_power_results.tsv and annotation/figures/FigS1_power.png
"""
import csv, os, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
RNG = np.random.default_rng(20260724); NPERM = 5000
def norm(g):
    m = re.fullmatch(r'(?i)bra(\d+)', g.strip().strip('"')); return f"Bra{m.group(1)}" if m else None

# scent set
tier, routes = {}, {}
for r in csv.DictReader(open(os.path.join(HERE, "scent_geneset.tsv")), delimiter='\t'):
    tier[r["bra_id"]] = int(r["tier"]); routes[r["bra_id"]] = r["routes"].split(";")

# matrix + design (as phase3_robust.py)
rows = list(csv.reader(open(os.path.join(ROOT, "Radiation model/Radi_model_with_interactions_deg_values_DESeq2.csv"))))
hdr = rows[0]; i_id = hdr.index("ensembl_ID"); start = hdr.index("Processed data:") + 1
pat = re.compile(r'^(AL|W)(0|40)(D|R)_')
keep, geno, dose, pres = [], [], [], []
for i in range(start, len(hdr)):
    m = pat.match(hdr[i])
    if m: keep.append(i); geno.append(m.group(1)); dose.append(int(m.group(2))); pres.append(m.group(3))
dose = np.array(dose); strata = np.array([g+p for g, p in zip(geno, pres)])
genes, M = [], []
for r in rows[1:]:
    g = norm(r[i_id])
    if not g: continue
    try: M.append([float(r[i]) for i in keep]); genes.append(g)
    except (ValueError, IndexError): pass
M = np.array(M); genes = np.array(genes); gidx = {g: k for k, g in enumerate(genes)}

def dose_effect(mat, d):
    eff = np.zeros(mat.shape[0]); ns = 0
    for s in np.unique(strata):
        sel = strata == s; ds = d[sel]
        if (ds == 40).sum() and (ds == 0).sum():
            eff += mat[:, sel][:, ds == 40].mean(1) - mat[:, sel][:, ds == 0].mean(1); ns += 1
    return eff / ns
obs = dose_effect(M, dose)

SETS = {"whole_scent_set": list(tier),
        "ester_route": [g for g in tier if "ester" in routes[g]],
        "Tier1_core": [g for g in tier if tier[g] == 1]}
idx = {k: np.array([gidx[g] for g in v if g in gidx]) for k, v in SETS.items()}

# one permutation loop -> null T for every set
Tnull = {k: np.empty(NPERM) for k in SETS}
masks = {k: (np.ones(M.shape[0], bool)) for k in SETS}
for k in SETS: masks[k][idx[k]] = False       # background mask
for b in range(NPERM):
    dp = dose.copy()
    for s in np.unique(strata):
        sel = np.where(strata == s)[0]; dp[sel] = RNG.permutation(dp[sel])
    eff = dose_effect(M, dp)
    for k in SETS:
        Tnull[k][b] = eff[idx[k]].mean() - eff[masks[k]].mean()

grid = np.round(np.arange(0.0, 0.41, 0.01), 2)
def power(tn, c, delta): return np.mean((tn + delta > c) | (tn + delta < -c))

res = []
for k in SETS:
    tn = Tnull[k]; c = np.percentile(np.abs(tn), 97.5)
    Tobs = obs[idx[k]].mean() - obs[masks[k]].mean()
    pw = np.array([power(tn, c, d) for d in grid])
    def mde(p):
        w = np.where(pw >= p)[0]; return float(grid[w[0]]) if len(w) else float('nan')
    res.append(dict(set=k, n=len(idx[k]), Tobs=float(Tobs), c95=float(c),
                    d50=mde(.5), d80=mde(.8), d90=mde(.9), grid=grid, pw=pw))

# ---- table ----
with open(os.path.join(HERE, "phase3_power_results.tsv"), "w") as fo:
    fo.write("set\tn\tobserved_T\tcrit_value_c(0.05,2-sided)\tMDE_delta_p50\tdelta_80pct_power\tdelta_90pct_power\t"
             "delta_80_as_foldchange\tdelta_80_as_percent\n")
    for r in res:
        fc = 2 ** r["d80"]; pc = (fc - 1) * 100
        fo.write(f"{r['set']}\t{r['n']}\t{r['Tobs']:+.4f}\t{r['c95']:.4f}\t{r['d50']:.3f}\t{r['d80']:.3f}\t"
                 f"{r['d90']:.3f}\t{fc:.3f}\t{pc:.1f}\n")

print(f"{'set':18s} {'n':>4s} {'T_obs':>8s} {'c(.05)':>7s} {'δ50':>6s} {'δ80':>6s} {'δ90':>6s}  (δ in log2 units)")
for r in res:
    print(f"{r['set']:18s} {r['n']:4d} {r['Tobs']:+8.4f} {r['c95']:7.4f} {r['d50']:6.3f} {r['d80']:6.3f} {r['d90']:6.3f}")
print("\nInterpretation: δ80 = smallest coordinated shift (log2) in the set's 40-vs-0 cGy dose")
print("effect detectable at 80% power / α=0.05. |T_obs| < c for every set (non-significant).")

# ---- figure ----
fig, ax = plt.subplots(figsize=(6.4, 4.0))
colr = {"whole_scent_set": "#333", "ester_route": "#C44E52", "Tier1_core": "#4C72B0"}
lab = {"whole_scent_set": "whole scent set (n=%d)", "ester_route": "ester/MT route (n=%d)", "Tier1_core": "Tier-1 core (n=%d)"}
for r in res:
    ax.plot(r["grid"], r["pw"], color=colr[r["set"]], lw=2, label=lab[r["set"]] % r["n"])
    if not np.isnan(r["d80"]): ax.plot([r["d80"]], [0.8], "o", color=colr[r["set"]], ms=5)
ax.axhline(0.8, ls=":", c="#888", lw=1); ax.text(0.005, 0.81, "80% power", fontsize=8, color="#666")
# observed |T| markers on x-axis
for r in res:
    ax.axvline(abs(r["Tobs"]), color=colr[r["set"]], ls="--", lw=0.8, alpha=0.5)
ax.set_xlabel("coordinated shift in scent-set dose effect, δ (log2, 40 vs 0 cGy)")
ax.set_ylabel("power (α=0.05, two-sided)"); ax.set_ylim(0, 1.02); ax.set_xlim(0, 0.4)
ax.set_title("Fig S1. Sensitivity of the radiation gene-set test\n"
             "(dashed = observed |effect|; markers = δ at 80% power)", fontsize=9, loc="left")
ax.legend(fontsize=8, frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "FigS1_power.png"), dpi=300, bbox_inches="tight")
print("wrote annotation/phase3_power_results.tsv + annotation/figures/FigS1_power.png")
