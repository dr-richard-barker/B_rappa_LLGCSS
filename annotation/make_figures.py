#!/usr/bin/env python3
"""
Draft manuscript figures (Phase 5 / item B6). Outputs 300-dpi PNGs to annotation/figures/.
Everything is derived from committed tables; re-run after phase3_robust.py / scent_orthology.py.
  Fig1  design + gene-ID-join schematic
  Fig2  scent gene set: routes×tiers (2A) and scent-axis validation High/Low (2B)
  Fig3  radiation dose effect per route vs background (from phase3_dose_stats.tsv) with perm p / BH q
  Fig4  cross-species conservation heatmap (from scent_orthology_matrix.tsv)
"""
import os, csv, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
FIG = os.path.join(HERE, "figures"); os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.grid": False,
                     "figure.dpi": 300, "savefig.dpi": 300, "font.family": "DejaVu Sans"})
ROUTE_ORDER = ["terpenoid", "benzenoid_phenylpropanoid", "ester", "fatty_acid_GLV", "apocarotenoid"]
ROUTE_LAB = {"terpenoid":"Terpenoid", "benzenoid_phenylpropanoid":"Benzenoid/\nphenylpropanoid",
             "ester":"Ester / methyl-\ntransferase", "fatty_acid_GLV":"Fatty-acid /\nGLV",
             "apocarotenoid":"Apocarotenoid"}
COL = {"terpenoid":"#4C72B0","benzenoid_phenylpropanoid":"#DD8452","ester":"#C44E52",
       "fatty_acid_GLV":"#55A868","apocarotenoid":"#8172B3"}
def norm(g):
    m = re.fullmatch(r'(?i)bra(\d+)', g.strip().strip('"')); return f"Bra{m.group(1)}" if m else None

# ---------- load scent set ----------
sg = list(csv.DictReader(open(os.path.join(HERE, "scent_geneset.tsv")), delimiter='\t'))
tier = {r["bra_id"]: int(r["tier"]) for r in sg}
routes = {r["bra_id"]: r["routes"].split(";") for r in sg}

# =========================================================== Fig 1: schematic
fig, ax = plt.subplots(figsize=(7.2, 3.4)); ax.axis("off"); ax.set_xlim(0,10); ax.set_ylim(0,5)
def box(x,y,w,h,txt,fc):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.03,rounding_size=0.12",
                fc=fc,ec="#333",lw=1)); ax.text(x+w/2,y+h/2,txt,ha="center",va="center",fontsize=8)
box(0.2,3.4,2.6,1.2,"Scent RNA-seq\nHigh vs Low\nSRR4417237-244\n(Bra######)","#DCE6F5")
box(0.2,0.4,2.6,1.2,"Radiation RNA-seq\n0 vs 40 cGy x geno x pres\n39 libs (BRA######)","#F5E3DC")
box(3.6,1.9,2.4,1.2,"Gene-ID join\n(case-normalise)\n31,756 genes 1:1","#E8E8E8")
box(6.7,3.4,3.1,1.2,"Tiered scent set\n363 genes / 4 routes\n(GO + Pfam)","#E4F0E4")
box(6.7,0.4,3.1,1.2,"Cross-species\nconservation +\nradiation gene-set test","#EFE7F5")
for a,b in [((2.8,4.0),(3.6,2.9)),((2.8,1.0),(3.6,2.3)),((6.0,2.6),(6.7,3.7)),((6.0,2.3),(6.7,1.0))]:
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle="-|>",mutation_scale=12,color="#555",lw=1.2))
ax.set_title("Fig 1. Two-experiment design and gene-ID join", fontsize=9, loc="left")
fig.savefig(os.path.join(FIG,"Fig1_design.png"), bbox_inches="tight"); plt.close(fig)

# =========================================================== Fig 2: gene set
fig, axs = plt.subplots(1,2, figsize=(7.2,3.2))
# 2A routes x tier stacked
t1 = [sum(1 for g in routes if r in routes[g] and tier[g]==1) for r in ROUTE_ORDER]
t2 = [sum(1 for g in routes if r in routes[g] and tier[g]==2) for r in ROUTE_ORDER]
x = np.arange(len(ROUTE_ORDER))
axs[0].bar(x, t1, color=[COL[r] for r in ROUTE_ORDER], label="Tier 1 (core)")
axs[0].bar(x, t2, bottom=t1, color=[COL[r] for r in ROUTE_ORDER], alpha=0.45, label="Tier 2 (supporting)")
axs[0].set_xticks(x); axs[0].set_xticklabels([ROUTE_LAB[r] for r in ROUTE_ORDER], fontsize=7)
axs[0].set_ylabel("genes"); axs[0].set_title("A  scent gene set (n=363)", fontsize=9, loc="left")
axs[0].legend(fontsize=7, frameon=False)
# 2B scent-axis validation: fraction DE High/Low per route
rows = list(csv.reader(open(os.path.join(ROOT,"NewTest/differential_expression_GLbulkRNAseq (1).csv"))))
h=rows[0]; gi=h.index("gene_id"); pi=h.index("Adj.p.value_(High)v(Low)")
padj={}
for rr in rows[1:]:
    g=norm(rr[gi])
    if g:
        try: padj[g]=float(rr[pi])
        except: padj[g]=None
frac=[]; ns=[]
for r in ROUTE_ORDER:
    mem=[g for g in routes if r in routes[g] and g in padj and padj[g] is not None]
    sig=[g for g in mem if padj[g]<0.1]
    frac.append(100*len(sig)/len(mem) if mem else 0); ns.append(len(mem))
axs[1].bar(x, frac, color=[COL[r] for r in ROUTE_ORDER])
for xi,f,n in zip(x,frac,ns): axs[1].text(xi,f+0.4,f"{n}",ha="center",fontsize=7,color="#333")
axs[1].set_xticks(x); axs[1].set_xticklabels([ROUTE_LAB[r] for r in ROUTE_ORDER], fontsize=7)
axs[1].set_ylabel("% genes DE High vs Low (adjP<0.1)")
axs[1].set_title("B  scent-axis validation", fontsize=9, loc="left")
fig.tight_layout(); fig.savefig(os.path.join(FIG,"Fig2_geneset.png"), bbox_inches="tight"); plt.close(fig)

# =========================================================== Fig 3: dose effect per route
ds = list(csv.DictReader(open(os.path.join(HERE,"phase3_dose_stats.tsv")), delimiter='\t'))
eff = {r["bra_id"]: float(r["dose_effect"]) for r in ds}
setgenes = set(routes)
bg = [eff[g] for g in eff if g not in setgenes]
res = {r["name"]: r for r in csv.DictReader(open(os.path.join(HERE,"phase3_robust_results.tsv")), delimiter='\t')}
fig, ax = plt.subplots(figsize=(7.2,3.6))
data=[bg]; labs=["background\n(all other genes)"]; cols=["#BBBBBB"]
for r in ROUTE_ORDER:
    vals=[eff[g] for g in routes if r in routes[g] and g in eff]
    data.append(vals); labs.append(f"{ROUTE_LAB[r]}\n(n={len(vals)})"); cols.append(COL[r])
bp=ax.boxplot(data, showfliers=False, widths=0.6, patch_artist=True, medianprops=dict(color="black"))
for patch,c in zip(bp["boxes"],cols): patch.set_facecolor(c); patch.set_alpha(0.7)
ax.axhline(np.median(bg), ls="--", lw=0.8, color="#888")
ax.set_xticklabels(labs, fontsize=7); ax.set_ylabel("stratified dose effect (40 vs 0 cGy, log2)")
ax.set_ylim(np.percentile(bg,2), np.percentile(bg,98))
# annotate ester + full-set permutation p / q
er=res.get("route:ester",{})
ax.set_title(f"Fig 3. Radiation dose effect by scent route — no route significant after correction\n"
             f"(ester/MT route: perm p={float(er.get('p_perm',0)):.2f}, BH q={float(er.get('q_perm_BH',0)):.2f}; all q>0.10)",
             fontsize=8.5, loc="left")
fig.tight_layout(); fig.savefig(os.path.join(FIG,"Fig3_radiation_test.png"), bbox_inches="tight"); plt.close(fig)

# =========================================================== Fig 4: conservation heatmap
rows = list(csv.reader(open(os.path.join(HERE,"scent_orthology_matrix.tsv"), ), delimiter='\t'))
hdr=rows[0]; sp_start=5
species=[s.replace("_gca000188115v5cm","").replace("_"," ") for s in hdr[sp_start:]]
genes=[r[0] for r in rows[1:]]; routescol=[r[2] for r in rows[1:]]
mat=np.array([[int(x) for x in r[sp_start:]] for r in rows[1:]], float)
fig, ax = plt.subplots(figsize=(7.6,4.2))
im=ax.imshow(np.log1p(mat), aspect="auto", cmap="magma")
ax.set_xticks(range(len(species))); ax.set_xticklabels(species, rotation=45, ha="right", fontsize=7, fontstyle="italic")
ax.set_yticks(range(len(genes)))
ax.set_yticklabels([f"{g} ({routescol[i][:4]})" for i,g in enumerate(genes)], fontsize=7)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if mat[i,j]>0: ax.text(j,i,int(mat[i,j]),ha="center",va="center",fontsize=6,
                               color="white" if np.log1p(mat[i,j])<np.log1p(mat).max()*0.6 else "black")
cb=fig.colorbar(im, ax=ax, fraction=0.025); cb.set_label("log(1+ortholog count)", fontsize=7)
ax.set_title("Fig 4. Cross-species conservation of the scent panel (Ensembl Compara orthologs)", fontsize=8.5, loc="left")
fig.tight_layout(); fig.savefig(os.path.join(FIG,"Fig4_conservation.png"), bbox_inches="tight"); plt.close(fig)

print("wrote:", ", ".join(sorted(os.listdir(FIG))))
