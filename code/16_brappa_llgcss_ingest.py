#!/usr/bin/env python3
"""
16_brappa_llgcss_ingest.py
Ingest and harmonize the Brassica rapa (Lunar LEAF / LLGCSS) GCR radiation dataset.
Maps B. rapa gene identifiers (BRA###### / Bra######) to Arabidopsis thaliana orthologs (TAIR10 AGI),
extracts factorial radiation contrasts, and formats query signatures for cross-species decoding.

Factorial design:
  - Radiation Condition: 40 cGy GCR vs 0 cGy
  - Genotype: Wild-Type (WT) vs anthocyaninless (anl)
  - Preservative: Deep Space Preservative (DRS) vs RNAlater (RL)

Author: Antigravity Assistant & Dr. Richard Barker
Date: November 2026
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Base paths
WORKSPACE_DIR = Path(__file__).resolve().parents[1]
LLGCSS_DIR = Path("/Users/drb_laptop/Documents/B_rappa_LLGCSS")
DATA_OUT = WORKSPACE_DIR / "data"
DATA_OUT.mkdir(exist_ok=True)

print("=" * 70)
print("16_brappa_llgcss_ingest: Starting B. rapa LLGCSS pure radiation ingestion")
print("=" * 70)

# 1. Load Orthology Map (18,770 RBH pairs)
ortho_file = LLGCSS_DIR / "annotation" / "brapa_to_arabidopsis_orthologs.tsv"
print(f"Loading orthology map from: {ortho_file}")
df_ortho = pd.read_csv(ortho_file, sep="\t")
df_ortho["BRA_id"] = df_ortho["brapa_gene"].str.strip().str.upper()
df_ortho["bra_id"] = df_ortho["brapa_gene"].str.strip().str.lower().str.capitalize()
df_ortho["arabidopsis_agi"] = df_ortho["arabidopsis_agi"].str.strip().str.upper()
print(f"Loaded {len(df_ortho):,} 1:1 ortholog pairs between B. rapa and A. thaliana.")

# 2. Load Radiation Model Results (DESeq2 main contrasts)
rad_results_file = LLGCSS_DIR / "Radiation model" / "Results_LFC_Pval_DESeq2.csv"
print(f"Loading radiation contrasts from: {rad_results_file}")
df_rad = pd.read_csv(rad_results_file)
df_rad["BRA_id"] = df_rad["ensembl_ID"].astype(str).str.strip().str.upper()

# 3. Load Interaction Terms Model
rad_inter_file = LLGCSS_DIR / "Radiation model" / "Radi_model_with_interactions_deg_values_DESeq2.csv"
print(f"Loading interaction terms from: {rad_inter_file}")
df_inter = pd.read_csv(rad_inter_file)
df_inter["BRA_id"] = df_inter["ensembl_ID"].astype(str).str.strip().str.upper()

# Select essential columns from radiation model
rad_cols = {
    "BRA_id": "BRA_id",
    "baseMean": "baseMean",
    "WT_&_40_CGY-WT_&_0_CGY_log2FC": "WT_40cGy_vs_0cGy_log2FC",
    "WT_&_40_CGY-WT_&_0_CGY_adjPval": "WT_40cGy_vs_0cGy_padj",
    "ANTHO_LESS_&_40_CGY-ANTHO_LESS_&_0_CGY_log2FC": "antho_less_40cGy_vs_0cGy_log2FC",
    "ANTHO_LESS_&_40_CGY-ANTHO_LESS_&_0_CGY_adjPval": "antho_less_40cGy_vs_0cGy_padj",
    "ANTHO_LESS_&_0_CGY-WT_&_0_CGY_log2FC": "antho_less_vs_WT_0cGy_log2FC",
    "ANTHO_LESS_&_0_CGY-WT_&_0_CGY_adjPval": "antho_less_vs_WT_0cGy_padj",
    "ANTHO_LESS_&_40_CGY-WT_&_40_CGY_log2FC": "antho_less_vs_WT_40cGy_log2FC",
    "ANTHO_LESS_&_40_CGY-WT_&_40_CGY_adjPval": "antho_less_vs_WT_40cGy_padj"
}
rad_extracted = df_rad[[c for c in rad_cols.keys() if c in df_rad.columns]].rename(columns=rad_cols)

# Select interaction columns
inter_cols = {
    "BRA_id": "BRA_id",
    "I:condition_40_CGY.genotype_ANTHO_LESS_log2FC": "inter_dose_x_genotype_log2FC",
    "I:condition_40_CGY.genotype_ANTHO_LESS_adjPval": "inter_dose_x_genotype_padj",
    "I:condition_40_CGY.preservative_RL_log2FC": "inter_dose_x_preservative_log2FC",
    "I:condition_40_CGY.preservative_RL_adjPval": "inter_dose_x_preservative_padj"
}
inter_extracted = df_inter[[c for c in inter_cols.keys() if c in df_inter.columns]].rename(columns=inter_cols)

# Merge radiation contrasts
rad_merged = pd.merge(rad_extracted, inter_extracted, on="BRA_id", how="outer")
print(f"Radiation dataset contains {len(rad_merged):,} genes.")

# 4. Master Integration: Orthology + Radiation Contrasts
master = pd.merge(rad_merged, df_ortho, on="BRA_id", how="left")
master["has_arabidopsis_ortholog"] = ~master["arabidopsis_agi"].isna()

print(f"\nMaster radiation table summary:")
print(f"  Total genes: {len(master):,}")
print(f"  Genes with Arabidopsis TAIR10 ortholog: {master['has_arabidopsis_ortholog'].sum():,}")

# Export master dataset
out_master_path = DATA_OUT / "brappa_master_radiation_contrasts.csv"
master.to_csv(out_master_path, index=False)
print(f"\n[Saved] Master radiation table -> {out_master_path}")

# Export specific query signatures for DSRS and Redox CVAE
# 1. WT 40 cGy vs 0 cGy (Arabidopsis AGI keyed)
sig_wt = master.dropna(subset=["arabidopsis_agi", "WT_40cGy_vs_0cGy_log2FC"])[
    ["arabidopsis_agi", "WT_40cGy_vs_0cGy_log2FC", "WT_40cGy_vs_0cGy_padj"]
].drop_duplicates("arabidopsis_agi").rename(columns={
    "arabidopsis_agi": "gene_id",
    "WT_40cGy_vs_0cGy_log2FC": "log2FC",
    "WT_40cGy_vs_0cGy_padj": "padj"
})
sig_wt_path = DATA_OUT / "brappa_sig_wt_gcr_40cGy.csv"
sig_wt.to_csv(sig_wt_path, index=False)
print(f"[Saved] WT 40 cGy GCR signature ({len(sig_wt):,} genes) -> {sig_wt_path}")

# 2. Anthocyaninless 40 cGy vs 0 cGy
sig_anl = master.dropna(subset=["arabidopsis_agi", "antho_less_40cGy_vs_0cGy_log2FC"])[
    ["arabidopsis_agi", "antho_less_40cGy_vs_0cGy_log2FC", "antho_less_40cGy_vs_0cGy_padj"]
].drop_duplicates("arabidopsis_agi").rename(columns={
    "arabidopsis_agi": "gene_id",
    "antho_less_40cGy_vs_0cGy_log2FC": "log2FC",
    "antho_less_40cGy_vs_0cGy_padj": "padj"
})
sig_anl_path = DATA_OUT / "brappa_sig_antho_less_gcr_40cGy.csv"
sig_anl.to_csv(sig_anl_path, index=False)
print(f"[Saved] Antho_less 40 cGy GCR signature ({len(sig_anl):,} genes) -> {sig_anl_path}")

# 3. Interaction Condition 40 cGy x Genotype Antho_less
sig_inter_geno = master.dropna(subset=["arabidopsis_agi", "inter_dose_x_genotype_log2FC"])[
    ["arabidopsis_agi", "inter_dose_x_genotype_log2FC", "inter_dose_x_genotype_padj"]
].drop_duplicates("arabidopsis_agi").rename(columns={
    "arabidopsis_agi": "gene_id",
    "inter_dose_x_genotype_log2FC": "log2FC",
    "inter_dose_x_genotype_padj": "padj"
})
sig_inter_geno_path = DATA_OUT / "brappa_sig_inter_dose_genotype.csv"
sig_inter_geno.to_csv(sig_inter_geno_path, index=False)
print(f"[Saved] Interaction dose x genotype signature ({len(sig_inter_geno):,} genes) -> {sig_inter_geno_path}")

# 4. Interaction Condition 40 cGy x Preservative RL (vs DRS)
sig_inter_pres = master.dropna(subset=["arabidopsis_agi", "inter_dose_x_preservative_log2FC"])[
    ["arabidopsis_agi", "inter_dose_x_preservative_log2FC", "inter_dose_x_preservative_padj"]
].drop_duplicates("arabidopsis_agi").rename(columns={
    "arabidopsis_agi": "gene_id",
    "inter_dose_x_preservative_log2FC": "log2FC",
    "inter_dose_x_preservative_padj": "padj"
})
sig_inter_pres_path = DATA_OUT / "brappa_sig_inter_dose_preservative.csv"
sig_inter_pres.to_csv(sig_inter_pres_path, index=False)
print(f"[Saved] Interaction dose x preservative signature ({len(sig_inter_pres):,} genes) -> {sig_inter_pres_path}")

# Delete old scent files if present
for old_file in [DATA_OUT / "brappa_sig_scent_high_vs_low.csv", DATA_OUT / "brappa_scent_contrasts.csv"]:
    if old_file.exists():
        old_file.unlink()
        print(f"[Deleted] Scent file removed -> {old_file}")

print("=" * 70)
print("16_brappa_llgcss_ingest: Completed successfully without floral scent!")
print("=" * 70)
