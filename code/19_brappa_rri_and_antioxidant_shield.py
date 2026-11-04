#!/usr/bin/env python3
"""
19_brappa_rri_and_antioxidant_shield.py
Phase 3 — Compute the Radiation Resilience Index (RRI), evaluate the Antioxidant
Shield Deficit (WT vs anthocyaninless), and assess Spaceflight Preservative
Fidelity (Deep Space Preservative DRS vs RNAlater RL).

Evaluates:
1. Sample-level RRI (latent stability, pathway balance, and housekeeping module preservation)
2. Genotype divergence: WT vs anthocyaninless mutant under 40 cGy GCR
3. Preservative fidelity: DRS vs RNAlater transcriptomic parity
4. Core DNA double-strand-break (DSB) repair network (BRCA1, PARP2, RAD51, SMR7)
5. Phenylpropanoid / flavonoid protective antioxidant shunt (PAL, CHS, CHI, F3H)

Author: Antigravity Assistant & Dr. Richard Barker
Date: November 2026
"""

import os
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import entropy, spearmanr, pearsonr

WORKSPACE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = WORKSPACE_DIR / "data"
RESULTS_DIR = WORKSPACE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

print("=" * 70)
print("19_brappa_rri_and_antioxidant_shield: Computing RRI, Preservative & DDR Metrics")
print("=" * 70)

# 1. Load Expression Data & Metadata from B. rapa
deseq_file = Path("/Users/drb_laptop/Documents/B_rappa_LLGCSS/Radiation model/Results_LFC_Pval_DESeq2.csv")
meta_file = Path("/Users/drb_laptop/Documents/B_rappa_LLGCSS/counts and factors/LEAF_Metadata.csv")
master_file = DATA_DIR / "brappa_master_radiation_contrasts.csv"

print(f"Loading DESeq2 expression matrix from: {deseq_file}")
df_deseq = pd.read_csv(deseq_file)
df_meta = pd.read_csv(meta_file)
df_master = pd.read_csv(master_file)

# Identify sample columns
sample_cols = [c for c in df_deseq.columns if "GLbulkRNAseq" in c]
print(f"Identified {len(sample_cols)} RNA-seq sample libraries in DESeq2 matrix.")

# Convert sample columns to float and fill NaN with 0.0
for c in sample_cols:
    df_deseq[c] = pd.to_numeric(df_deseq[c], errors="coerce").fillna(0.0)

# Clean gene IDs
df_deseq["BRA_id"] = df_deseq["ensembl_ID"].astype(str).str.strip().str.upper()

# Merge with ortholog annotations
expr_annot = pd.merge(
    df_deseq[["BRA_id"] + sample_cols],
    df_master[["BRA_id", "arabidopsis_agi", "WT_40cGy_vs_0cGy_log2FC", "antho_less_40cGy_vs_0cGy_log2FC", "inter_dose_x_genotype_log2FC", "inter_dose_x_preservative_log2FC"]],
    on="BRA_id",
    how="left"
)

# 2. Extract Subsets for RRI computation
row_sums = expr_annot[sample_cols].sum(axis=1)
expr_expressed = expr_annot[row_sums > 0].copy()

# Ortholog space for latent distance
ortho_expr = expr_expressed.dropna(subset=["arabidopsis_agi"]).drop_duplicates("arabidopsis_agi")
X = ortho_expr[sample_cols].values.T # (31 samples x N genes)

# Parse metadata per sample
meta_records = []
for s in sample_cols:
    row = df_meta[df_meta["sample_id"] == s]
    if not row.empty:
        cond = row["condition"].values[0]
        pres = row["preservative"].values[0]
        geno = row["genotype"].values[0]
    else:
        cond = "40_cGy" if "40" in s else "0_cGy"
        pres = "RL" if "R" in s else "DRS"
        geno = "antho_less" if s.startswith("AL") else "WT"
    is_ctrl = (cond == "0_cGy") and (geno == "WT")
    meta_records.append({
        "sample_id": s,
        "condition": cond,
        "preservative": pres,
        "genotype": geno,
        "is_control": is_ctrl
    })
meta_df = pd.DataFrame(meta_records)

# 3. Compute RRI Components
ctrl_mask = meta_df["is_control"].values
if np.sum(ctrl_mask) == 0:
    ctrl_mask = (meta_df["condition"] == "0_cGy").values

ctrl_X = X[ctrl_mask]
ctrl_centroid = np.mean(ctrl_X, axis=0, keepdims=True)

# Pairwise distances among controls for sigma
ctrl_dists = cdist(ctrl_X, ctrl_X, metric="euclidean")
sigma = np.median(ctrl_dists[ctrl_dists > 0]) or 1.0

# Distance of each sample to control centroid
sample_dists = cdist(X, ctrl_centroid, metric="euclidean").flatten()
rri_latent = np.exp(- (sample_dists ** 2) / (2 * (sigma ** 2)))
rri_latent = np.clip(rri_latent, 0.0, 1.0)

# Component 2: Pathway Balance RRI (DDR, ROS Defense, Flavonoids)
pathway_genes = {
    "dsb_repair": ["AT4G21070", "AT4G02390", "AT5G20850", "AT3G27630", "AT3G02820", "AT5G41650"],
    "ros_defense": ["AT1G07890", "AT4G35090", "AT1G20630", "AT1G08830", "AT4G25100", "AT3G22370", "AT3G24170"],
    "phenylpropanoid": ["AT5G13930", "AT3G53260", "AT5G42800", "AT3G51240", "AT2G37040", "AT3G55120"]
}

pathway_scores_per_sample = []
for s in sample_cols:
    scores = []
    for pname, glist in pathway_genes.items():
        sub = expr_expressed[expr_expressed["arabidopsis_agi"].isin(glist)]
        scores.append(sub[s].mean() if not sub.empty else 1.0)
    
    p_vec = np.array(scores)
    p_norm = p_vec / np.sum(p_vec) if np.sum(p_vec) > 0 else np.ones(len(scores)) / len(scores)
    h = entropy(p_norm)
    h_max = np.log(len(scores))
    evenness = h / h_max if h_max > 0 else 1.0
    pathway_scores_per_sample.append(evenness)

rri_pathway = np.clip(np.array(pathway_scores_per_sample), 0.0, 1.0)

# Component 3: Housekeeping Module Preservation RRI
ctrl_var = np.var(ctrl_X, axis=0)
hk_idx = np.argsort(ctrl_var)[:int(0.2 * len(ctrl_var))]
ctrl_hk_mean = np.mean(ctrl_X[:, hk_idx], axis=0)

rri_module = []
for i in range(len(sample_cols)):
    r, _ = pearsonr(X[i, hk_idx], ctrl_hk_mean)
    rri_module.append(max(0.0, r))
rri_module = np.array(rri_module)

# Composite RRI: 0.50 latent + 0.25 pathway + 0.25 module
rri_total = 0.50 * rri_latent + 0.25 * rri_pathway + 0.25 * rri_module

meta_df["RRI"] = rri_total
meta_df["RRI_latent"] = rri_latent
meta_df["RRI_pathway"] = rri_pathway
meta_df["RRI_module"] = rri_module

print(f"\nRRI Summary Statistics by Genotype and Dose:")
rri_summary = meta_df.groupby(["genotype", "condition"])[["RRI", "RRI_latent", "RRI_pathway", "RRI_module"]].mean()
print(rri_summary)

print(f"\nRRI Summary Statistics by Preservative Fixative:")
pres_summary = meta_df.groupby(["preservative", "condition"])[["RRI", "RRI_latent", "RRI_pathway", "RRI_module"]].mean()
print(pres_summary)

# 4. Preservative Fidelity & DDR Network Evaluation
# Compare DRS vs RL expression across radiation-responsive genes
wt_0_drs = [s for s in sample_cols if "W0D" in s]
wt_40_drs = [s for s in sample_cols if "W40D" in s]
wt_0_rl = [s for s in sample_cols if "W0R" in s]
wt_40_rl = [s for s in sample_cols if "W40R" in s]

drs_lfc = expr_expressed[wt_40_drs].mean(axis=1) - expr_expressed[wt_0_drs].mean(axis=1)
rl_lfc = expr_expressed[wt_40_rl].mean(axis=1) - expr_expressed[wt_0_rl].mean(axis=1)

pres_corr_r, pres_corr_p = pearsonr(drs_lfc, rl_lfc)
print(f"\nPreservative Fixative Parity (DRS vs RNAlater):")
print(f"  Pearson correlation across all expressed genes: r = {pres_corr_r:.4f} (p = {pres_corr_p:.2e})")

# Evaluate Core DNA Double-Strand Break Repair (DDR) Network
ddr_markers = [
    {"symbol": "BRCA1", "agi": "AT4G21070", "desc": "Breast cancer susceptibility 1 / DSB repair"},
    {"symbol": "PARP2", "agi": "AT4G02390", "desc": "Poly(ADP-ribose) polymerase 2 / DNA strand break sensor"},
    {"symbol": "RAD51", "agi": "AT5G20850", "desc": "DNA repair protein RAD51 / homologous recombination"},
    {"symbol": "SMR7", "agi": "AT3G27630", "desc": "SIAMESE-RELATED 7 / DNA damage cell cycle arrest"},
    {"symbol": "PAL1", "agi": "AT2G37040", "desc": "Phenylalanine ammonia-lyase / flavonoid antioxidant entry"},
    {"symbol": "CHS", "agi": "AT5G13930", "desc": "Chalcone synthase / anthocyanin gateway enzyme"},
    {"symbol": "F3H", "agi": "AT3G51240", "desc": "Flavanone 3-hydroxylase / flavonoid core biosynthesis"}
]

ddr_summary_records = []
for m in ddr_markers:
    sub = expr_annot[expr_annot["arabidopsis_agi"] == m["agi"]]
    if not sub.empty:
        wt_lfc = sub["WT_40cGy_vs_0cGy_log2FC"].values[0]
        anl_lfc = sub["antho_less_40cGy_vs_0cGy_log2FC"].values[0]
        inter_g = sub["inter_dose_x_genotype_log2FC"].values[0]
        inter_p = sub["inter_dose_x_preservative_log2FC"].values[0]
        ddr_summary_records.append({
            "symbol": m["symbol"],
            "agi": m["agi"],
            "description": m["desc"],
            "bra_id": sub["BRA_id"].values[0],
            "WT_40cGy_log2FC": wt_lfc,
            "AnthoLess_40cGy_log2FC": anl_lfc,
            "Inter_Genotype_log2FC": inter_g,
            "Inter_Preservative_log2FC": inter_p
        })

df_ddr = pd.DataFrame(ddr_summary_records)
print(f"\nCore DDR and Phenylpropanoid Defense Markers:")
print(df_ddr[["symbol", "WT_40cGy_log2FC", "AnthoLess_40cGy_log2FC", "Inter_Genotype_log2FC", "Inter_Preservative_log2FC"]])

# Save Output Tables
out_rri_path = RESULTS_DIR / "brappa_rri_sample_scores.csv"
out_ddr_path = RESULTS_DIR / "brappa_preservative_and_ddr_summary.csv"

meta_df.to_csv(out_rri_path, index=False)
df_ddr.to_csv(out_ddr_path, index=False)

# Clean up obsolete scent files
for stale in [WORKSPACE_DIR / "code" / "19_brappa_rri_and_scent_tradeoff.py", RESULTS_DIR / "brappa_scent_tradeoff_summary.csv"]:
    if stale.exists():
        stale.unlink()
        print(f"[Deleted] Removed obsolete scent file: {stale}")

print(f"\n[Saved] Sample RRI Scores -> {out_rri_path}")
print(f"[Saved] Preservative & DDR Summary -> {out_ddr_path}")
print("=" * 70)
print("19_brappa_rri_and_antioxidant_shield: Completed successfully!")
print("=" * 70)
