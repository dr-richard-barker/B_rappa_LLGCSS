#!/usr/bin/env python3
"""
Preliminary test of the core question:
  Does GCR radiation preferentially perturb floral-scent biosynthesis genes in B. rapa?

Uses the Phase-1 annotation (annotation/scent_geneset.tsv) to intersect:
  - scent DE      : NewTest/differential_expression... Log2fc/adjP (High vs Low)
  - radiation DE  : Radiation model/Results... WT_40cGy vs WT_0cGy (clean radiation contrast)
                    and antho_less 40 vs 0 cGy.
Reports overlap counts + a Fisher exact enrichment test (scent genes among radiation DEGs).
This is EXPLORATORY: our radiation counts are an in-house re-analysis; treat as a signal to
follow up, not a final result.
"""
import csv, os, math
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
FDR = 0.10

def num(x):
    try: return float(x)
    except: return None
def norm(g):
    g = g.strip().strip('"'); import re
    m = re.fullmatch(r'(?i)bra(\d+)', g); return f"Bra{m.group(1)}" if m else None

# scent-candidate gene set (Phase 2 tiered); keep genes present in both experiments
tier_of = {}; route_of = {}
with open(os.path.join(HERE, "scent_geneset.tsv")) as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        if row["in_scent_data"] == "1" and row["in_radiation_data"] == "1":
            tier_of[row["bra_id"]] = int(row["tier"]); route_of[row["bra_id"]] = row["routes"]
SET_T1 = {g for g, t in tier_of.items() if t == 1}          # specific
SET_ALL = set(tier_of)                                       # sensitive (Tier 1+2)
scent_set = SET_ALL

# radiation DE: WT 40 vs 0 cGy, and antho_less 40 vs 0 cGy
rad = {}
with open(os.path.join(ROOT, "Radiation model/Results_LFC_Pval_DESeq2.csv")) as fh:
    r = csv.reader(fh); h = next(r)
    iid = h.index("ensembl_ID")
    iWc = h.index("WT_&_40_CGY-WT_&_0_CGY_log2FC");  iWp = h.index("WT_&_40_CGY-WT_&_0_CGY_adjPval")
    iAc = h.index("ANTHO_LESS_&_40_CGY-ANTHO_LESS_&_0_CGY_log2FC"); iAp = h.index("ANTHO_LESS_&_40_CGY-ANTHO_LESS_&_0_CGY_adjPval")
    for row in r:
        n = norm(row[iid])
        if n: rad[n] = (num(row[iWc]), num(row[iWp]), num(row[iAc]), num(row[iAp]))

# scent DE: High vs Low
scentDE = {}
with open(os.path.join(ROOT, "NewTest/differential_expression_GLbulkRNAseq (1).csv")) as fh:
    r = csv.reader(fh); h = next(r)
    iid = h.index("gene_id"); ilfc = h.index("Log2fc_(High)v(Low)"); iap = h.index("Adj.p.value_(High)v(Low)")
    for row in r:
        n = norm(row[iid])
        if n: scentDE[n] = (num(row[ilfc]), num(row[iap]))

def sig(p): return p is not None and p < FDR

# background = genes with a testable WT-radiation adjP
bg = [g for g, v in rad.items() if v[1] is not None]
rad_deg_W = {g for g in bg if sig(rad[g][1])}
rad_deg_A = {g for g, v in rad.items() if v[3] is not None and sig(v[3])}
print(f"Radiation background (WT contrast, testable): {len(bg)}")
print(f"Radiation DEGs (WT 40 vs 0 cGy, adjP<{FDR}): {len(rad_deg_W)}")
print(f"Radiation DEGs (antho_less 40 vs 0 cGy):     {len(rad_deg_A)}")
print(f"Scent set (both experiments): Tier1 core={len(SET_T1)}  Tier1+2={len(SET_ALL)}")

# --- validate the gene set on the SCENT axis (should move in High vs Low) ---
def scent_axis(S):
    tested = [g for g in S if g in scentDE and scentDE[g][1] is not None]
    return sum(1 for g in tested if sig(scentDE[g][1])), len(tested)
for label, S in [("Tier1 core", SET_T1), ("Tier1+2", SET_ALL)]:
    s, n = scent_axis(S)
    print(f"[scent-axis validation] {label}: DE High vs Low (adjP<{FDR}) = {s}/{n}")
from collections import defaultdict as _dd
route_hit = _dd(lambda: [0, 0])
for g in SET_ALL:
    if g in scentDE and scentDE[g][1] is not None:
        for rt in route_of[g].split(";"):
            route_hit[rt][1] += 1
            if sig(scentDE[g][1]): route_hit[rt][0] += 1
print("  by route (DE/tested on scent axis):", {k: f"{v[0]}/{v[1]}" for k, v in route_hit.items()})

# --- enrichment: scent genes among radiation DEGs (Fisher exact) ---
def fisher(a, b, c, d):
    # 2x2: a=scent&DEG, b=scent&notDEG, c=notScent&DEG, d=notScent&notDEG; one-sided (enrichment)
    def logf(n): return math.lgamma(n + 1)
    def hyp(a, b, c, d):
        n = a + b + c + d
        return math.exp(logf(a+b)+logf(c+d)+logf(a+c)+logf(b+d)-logf(a)-logf(b)-logf(c)-logf(d)-logf(n))
    p = 0.0; amin = max(0, (a+c) - (c+d) if False else 0)
    r1 = a + b; c1 = a + c; n = a + b + c + d
    lo = max(0, c1 - (n - r1)); hi = min(r1, c1)
    p0 = hyp(a, b, c, d)
    for x in range(lo, hi + 1):
        px = hyp(x, r1 - x, c1 - x, n - r1 - (c1 - x))
        if px <= p0 * (1 + 1e-9): p += px
    return p

print("\n[enrichment] scent genes among radiation DEGs (Fisher exact, one-sided):")
for clabel, degset, use_wt in [("WT 40 vs 0 cGy", rad_deg_W, True), ("antho_less 40 vs 0 cGy", rad_deg_A, False)]:
    bgset = set(bg) if use_wt else {g for g, v in rad.items() if v[3] is not None}
    for slabel, Sfull in [("Tier1 core", SET_T1), ("Tier1+2", SET_ALL)]:
        S = Sfull & bgset
        a = len(S & degset); b = len(S) - a
        c = len(degset & bgset) - a; d = len(bgset) - a - b - c
        exp = len(S) * len(degset & bgset) / len(bgset) if bgset else 0
        p = fisher(a, b, c, d)
        fold = (a / len(S)) / (len(degset & bgset) / len(bgset)) if S and (degset & bgset) else float('nan')
        print(f"  [{clabel} | {slabel}] {a}/{len(S)} scent DEGs (exp ~{exp:.1f})  fold={fold:.2f}  Fisher p={p:.3g}")

# --- headline candidates: scent-associated AND radiation-responsive ---
print("\n[headline] genes that are BOTH scent-candidate AND radiation-DE (WT contrast):")
hits = []
for g in sorted(scent_set & rad_deg_W):
    wlfc, wp, _, _ = rad[g]
    slfc, sp = scentDE.get(g, (None, None))
    hits.append((g, wlfc, wp, slfc, sp))
print(f"  n = {len(hits)}")
for g, wlfc, wp, slfc, sp in hits[:25]:
    sflag = "scentDE" if sig(sp) else ""
    print(f"  {g}  radLFC={wlfc:+.2f} radAdjP={wp:.1e}  scentLFC={slfc if slfc is None else round(slfc,2)} {sflag}")
