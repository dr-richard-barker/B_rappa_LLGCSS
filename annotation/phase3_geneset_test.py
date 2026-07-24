#!/usr/bin/env python3
"""
Phase 3 — properly powered test: does radiation perturb the floral-scent gene set?

The Phase-1 preview used a Fisher test on the 13 hard DEGs — badly under-powered. This
uses a RANK-BASED competitive gene-set test (Mann-Whitney U of the scent set's radiation
effect vs the genome-wide background), which uses every gene's effect size, not a DEG
cutoff, so it detects a coordinated subtle shift if one exists.

Radiation effects (in-house counts, iDEP interaction model):
  main dose      40_CGY-0_CGY_log2FC                         (col 4/5)  -- all samples
  dose×genotype  I:condition_40_CGY.genotype_ANTHO_LESS      (col 8/9)  -- confound
  dose×preserv.  I:condition_40_CGY.preservative_RL          (col 6/7)  -- confound
Tests, per tier (Tier1 core / Tier1+2) and per route:
  1. two-sided MWU on signed log2FC  -> directional shift of scent genes
  2. two-sided MWU on |log2FC|       -> are scent genes more perturbed (any direction)
  3. sign test on log2FC             -> up/down imbalance
Also: how many scent genes have a significant dose×genotype / dose×preservative interaction
(does the scent radiation-response depend on the confounds?), and a ranked candidate list.
Caveat: in-house re-analysis, 40 cGy; a null here is 'no detectable coordinated effect'.
"""
import csv, os, math, re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
FDR = 0.10
def norm(g):
    m = re.fullmatch(r'(?i)bra(\d+)', g.strip().strip('"')); return f"Bra{m.group(1)}" if m else None
def num(x):
    try: return float(x)
    except: return None

# ---- scent set (tiered + routes) ----
tier, routes = {}, {}
with open(os.path.join(HERE, "scent_geneset.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        tier[r["bra_id"]] = int(r["tier"]); routes[r["bra_id"]] = r["routes"].split(";")
SET_T1 = {g for g,t in tier.items() if t==1}; SET_ALL = set(tier)

# ---- radiation effects (interaction model) ----
dose, dose_p, ixG, ixG_p, ixP, ixP_p = {}, {}, {}, {}, {}, {}
with open(os.path.join(ROOT, "Radiation model/Radi_model_with_interactions_deg_values_DESeq2.csv")) as fh:
    r = csv.reader(fh); h = next(r)
    i_id = h.index("ensembl_ID")
    c = {k: h.index(k) for k in ["40_CGY-0_CGY_log2FC","40_CGY-0_CGY_adjPval",
         "I:condition_40_CGY.genotype_ANTHO_LESS_log2FC","I:condition_40_CGY.genotype_ANTHO_LESS_adjPval",
         "I:condition_40_CGY.preservative_RL_log2FC","I:condition_40_CGY.preservative_RL_adjPval"]}
    for row in r:
        g = norm(row[i_id])
        if not g: continue
        v = num(row[c["40_CGY-0_CGY_log2FC"]])
        if v is not None:
            dose[g]=v; dose_p[g]=num(row[c["40_CGY-0_CGY_adjPval"]])
            ixG[g]=num(row[c["I:condition_40_CGY.genotype_ANTHO_LESS_log2FC"]]); ixG_p[g]=num(row[c["I:condition_40_CGY.genotype_ANTHO_LESS_adjPval"]])
            ixP[g]=num(row[c["I:condition_40_CGY.preservative_RL_log2FC"]]); ixP_p[g]=num(row[c["I:condition_40_CGY.preservative_RL_adjPval"]])

# ---- scent DE (High vs Low) for cross-tab ----
scentDE = {}
with open(os.path.join(ROOT, "NewTest/differential_expression_GLbulkRNAseq (1).csv")) as fh:
    r = csv.reader(fh); h = next(r); gi=h.index("gene_id"); li=h.index("Log2fc_(High)v(Low)"); pi=h.index("Adj.p.value_(High)v(Low)")
    for row in r:
        g=norm(row[gi])
        if g: scentDE[g]=(num(row[li]),num(row[pi]))

# ---- Mann-Whitney U (normal approx, tie-corrected), two-sided ----
def mwu(a, b):
    na, nb = len(a), len(b)
    if na == 0 or nb == 0: return None
    allv = sorted([(v,0) for v in a] + [(v,1) for v in b])
    # rank with ties (average ranks)
    ranks = [0.0]*len(allv); i = 0; tie_terms = 0.0
    while i < len(allv):
        j = i
        while j+1 < len(allv) and allv[j+1][0] == allv[i][0]: j += 1
        avg = (i + j)/2.0 + 1
        for k in range(i, j+1): ranks[k] = avg
        t = j - i + 1
        if t > 1: tie_terms += t**3 - t
        i = j + 1
    Ra = sum(ranks[k] for k in range(len(allv)) if allv[k][1]==0)
    Ua = Ra - na*(na+1)/2.0
    mu = na*nb/2.0
    N = na + nb
    sigma = math.sqrt(na*nb/12.0 * ((N+1) - tie_terms/(N*(N-1))))
    if sigma == 0: return (Ua, mu, 1.0)
    z = (Ua - mu)/sigma
    p = math.erfc(abs(z)/math.sqrt(2))          # two-sided
    return (Ua, mu, p, z)

def median(v):
    s=sorted(v); n=len(s); return (s[n//2] if n%2 else (s[n//2-1]+s[n//2])/2) if n else float('nan')

def run(setname, S):
    S = [g for g in S if g in dose]
    bg = [g for g in dose if g not in set(S)]
    sd = [dose[g] for g in S]; bd = [dose[g] for g in bg]
    print(f"\n[{setname}] n_scent={len(S)}  n_background={len(bg)}")
    m = mwu(sd, bd)
    print(f"  signed log2FC : median scent={median(sd):+.3f} vs bg={median(bd):+.3f}  MWU p={m[2]:.3g} (z={m[3]:+.2f})")
    ma = mwu([abs(x) for x in sd], [abs(x) for x in bd])
    print(f"  |log2FC|      : median scent={median([abs(x) for x in sd]):.3f} vs bg={median([abs(x) for x in bd]):.3f}  MWU p={ma[2]:.3g} (z={ma[3]:+.2f})")
    up = sum(1 for x in sd if x>0); dn = sum(1 for x in sd if x<0)
    # binomial two-sided sign test p (normal approx)
    n=up+dn; z=(up-n/2)/math.sqrt(n/4) if n else 0; ps=math.erfc(abs(z)/math.sqrt(2))
    print(f"  direction     : {up} up / {dn} down  sign-test p={ps:.3g}")
    return m, ma

print("="*70)
print("PHASE 3 — radiation main dose effect (40 vs 0 cGy) on the scent gene set")
print("="*70)
print(f"genes with a testable dose effect: {len(dose)}")
run("Tier1 core", SET_T1); run("Tier1+2", SET_ALL)
# per route (Tier1+2)
print("\n--- per route (Tier1+2, signed log2FC MWU vs background) ---")
for rt in ["terpenoid","benzenoid_phenylpropanoid","ester","fatty_acid_GLV","apocarotenoid"]:
    Srt = [g for g in SET_ALL if rt in routes.get(g,[]) and g in dose]
    bg = [dose[g] for g in dose if g not in set(Srt)]
    if len(Srt) < 3: print(f"  {rt:26s} n={len(Srt)} (too few)"); continue
    m = mwu([dose[g] for g in Srt], bg)
    print(f"  {rt:26s} n={len(Srt):3d}  median={median([dose[g] for g in Srt]):+.3f}  MWU p={m[2]:.3g}")

# ---- confounds: do scent genes' radiation response depend on genotype / preservative? ----
print("\n--- confound interactions among scent genes (adjP<{:.2f}) ---".format(FDR))
sg_ixG = [g for g in SET_ALL if ixG_p.get(g) is not None and ixG_p[g]<FDR]
sg_ixP = [g for g in SET_ALL if ixP_p.get(g) is not None and ixP_p[g]<FDR]
print(f"  scent genes with sig dose×genotype interaction:    {len(sg_ixG)}")
print(f"  scent genes with sig dose×preservative interaction: {len(sg_ixP)}")

# ---- candidate list: scent genes most responsive to dose, cross-ref High/Low ----
print("\n--- top scent genes by |dose log2FC| (with scent High/Low DE) ---")
ranked = sorted(SET_ALL, key=lambda g: abs(dose.get(g,0)), reverse=True)
print(f"  {'gene':11s} {'route':18s} {'doseLFC':>8s} {'doseAdjP':>9s} {'scentLFC':>8s} {'scentAdjP':>9s}")
for g in ranked[:15]:
    sl, sp = scentDE.get(g,(None,None))
    print(f"  {g:11s} {routes.get(g,['?'])[0][:18]:18s} {dose[g]:+8.2f} "
          f"{(dose_p.get(g) if dose_p.get(g) is not None else float('nan')):9.2g} "
          f"{('' if sl is None else f'{sl:+.2f}'):>8s} {('' if sp is None else f'{sp:.2g}'):>9s}")
