#!/usr/bin/env python3
"""
17_run_dsrs_classifier.py
Phase 2.1 — Run the DeepSpace Stress-Recognition System (DSRS) on Brassica rapa
GCR radiation contrasts.

Evaluates 4 pure radiation query signatures:
1. WT 40 cGy vs 0 cGy GCR
2. Anthocyaninless (anl) 40 cGy vs 0 cGy GCR
3. Dose x Genotype Interaction (excess radiation effect in anl)
4. Dose x Preservative Interaction (RL vs DRS fixative)

Projects the orthology-mapped B. rapa gene signatures onto the 123 seed developmental
panels and evaluates correlation with 27 reference spaceflight and physical stressor contrasts.

Author: Antigravity Assistant & Dr. Richard Barker
Date: November 2026
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

WORKSPACE_DIR = Path(__file__).resolve().parents[1]
DEEPSPACE_DIR = Path("/Users/drb_laptop/Documents/deepspace-seed-stress-decoder")
DATA_DIR = WORKSPACE_DIR / "data"
RESULTS_DIR = WORKSPACE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

print("=" * 70)
print("17_run_dsrs_classifier: Running DSRS Pure Radiation Stress-Recognition")
print("=" * 70)

# 1. Load Reference Panel Library
panel_file = DEEPSPACE_DIR / "panels" / "panel_library.csv"
print(f"Loading seed panel library from: {panel_file}")
df_panels = pd.read_csv(panel_file)
df_panels["panel"] = df_panels["panel_source"] + "::" + df_panels["panel_group"].astype(str)
panel_gene_sets = {p: set(g["gene"].dropna().unique()) for p, g in df_panels.groupby("panel")}
print(f"Loaded {len(panel_gene_sets)} panels.")

# 2. Load Reference NES Matrix & Stressor Classes
ref_nes_file = DEEPSPACE_DIR / "results" / "tables" / "decoder_nes_matrix_v7.csv"
class_file = DEEPSPACE_DIR / "results" / "tables" / "contrast_classes.csv"
print(f"Loading reference stressor NES matrix from: {ref_nes_file}")
ref_nes = pd.read_csv(ref_nes_file, index_col=0)
ref_classes = pd.read_csv(class_file, index_col=0)["stressor_class"].to_dict()
print(f"Reference matrix shape: {ref_nes.shape} ({len(ref_classes)} classified stressors)")

# GSEA Enrichment Score calculation
def compute_gsea_es(ranked_genes, ranked_vals, gene_set, p=1.0):
    n = len(ranked_genes)
    is_hit = np.isin(ranked_genes, list(gene_set))
    n_hits = np.sum(is_hit)
    if n_hits == 0 or n_hits == n:
        return 0.0
    
    weights = np.abs(ranked_vals) ** p
    hit_weights = weights * is_hit
    nr = np.sum(hit_weights)
    if nr == 0:
        return 0.0
    
    p_hit = np.cumsum(hit_weights) / nr
    p_miss = np.cumsum(1.0 - is_hit) / (n - n_hits)
    diff = p_hit - p_miss
    max_d = np.max(diff)
    min_d = np.min(diff)
    return float(max_d if max_d > -min_d else min_d)

def project_signature(sig_df, n_perms=60, seed=42):
    """Project signature onto 123 panels with permutation normalization."""
    np.random.seed(seed)
    sorted_sig = sig_df.sort_values("log2FC", ascending=False).dropna(subset=["gene_id", "log2FC"])
    ranked_genes = sorted_sig["gene_id"].values
    ranked_vals = sorted_sig["log2FC"].values
    
    panel_names = [p for p in ref_nes.index if p in panel_gene_sets]
    raw_es = {}
    for p in panel_names:
        raw_es[p] = compute_gsea_es(ranked_genes, ranked_vals, panel_gene_sets[p])
    
    # Fast permutation for NES scaling
    null_es = {p: [] for p in panel_names}
    for _ in range(n_perms):
        perm_vals = np.random.permutation(ranked_vals)
        for p in panel_names:
            null_es[p].append(compute_gsea_es(ranked_genes, perm_vals, panel_gene_sets[p]))
    
    nes = {}
    for p in panel_names:
        obs = raw_es[p]
        perms = np.array(null_es[p])
        if obs >= 0:
            pos_null = perms[perms >= 0]
            denom = np.mean(pos_null) if len(pos_null) > 0 and np.mean(pos_null) > 0 else 1.0
            nes[p] = obs / denom
        else:
            neg_null = perms[perms < 0]
            denom = np.abs(np.mean(neg_null)) if len(neg_null) > 0 and np.abs(np.mean(neg_null)) > 0 else 1.0
            nes[p] = obs / denom
            
    return pd.Series(nes)

# Pure radiation query signatures
queries = {
    "Brapa_WT_GCR_40cGy": DATA_DIR / "brappa_sig_wt_gcr_40cGy.csv",
    "Brapa_AnthoLess_GCR_40cGy": DATA_DIR / "brappa_sig_antho_less_gcr_40cGy.csv",
    "Brapa_Interaction_Dose_x_Genotype": DATA_DIR / "brappa_sig_inter_dose_genotype.csv",
    "Brapa_Interaction_Dose_x_Preservative": DATA_DIR / "brappa_sig_inter_dose_preservative.csv"
}

projected_nes_dict = {}
all_rankings = []
all_family_rankings = []

for qname, qpath in queries.items():
    print(f"\nProcessing query signature: {qname}...")
    qdf = pd.read_csv(qpath)
    print(f"  Input genes: {len(qdf):,}")
    q_nes = project_signature(qdf, n_perms=60)
    projected_nes_dict[qname] = q_nes
    
    # Calculate similarity with reference stressors
    common = q_nes.index.intersection(ref_nes.index)
    sims = ref_nes.loc[common].corrwith(q_nes.loc[common], method="spearman").sort_values(ascending=False)
    pears = ref_nes.loc[common].corrwith(q_nes.loc[common], method="pearson")
    
    res_df = pd.DataFrame({
        "query": qname,
        "stressor": sims.index,
        "stressor_class": [ref_classes.get(s, "Unknown") for s in sims.index],
        "spearman_rho": sims.values,
        "pearson_r": [pears[s] for s in sims.index]
    })
    all_rankings.append(res_df)
    
    # Aggregate by family
    fam_df = res_df.groupby("stressor_class").agg({
        "spearman_rho": "mean",
        "pearson_r": "mean"
    }).sort_values("spearman_rho", ascending=False).reset_index()
    fam_df["query"] = qname
    all_family_rankings.append(fam_df)
    
    print(f"  Top 3 Stressors:")
    for idx, row in res_df.head(3).iterrows():
        print(f"    {row['stressor']} ({row['stressor_class']}): rho = {row['spearman_rho']:.4f}")
    print(f"  Top 3 Families:")
    for idx, row in fam_df.head(3).iterrows():
        print(f"    {row['stressor_class']}: mean rho = {row['spearman_rho']:.4f}")

# Combine and save results
df_all_rankings = pd.concat(all_rankings, ignore_index=True)
df_all_family = pd.concat(all_family_rankings, ignore_index=True)
df_projected_nes = pd.DataFrame(projected_nes_dict)

out_rank_path = RESULTS_DIR / "dsrs_brappa_stressor_rankings.csv"
out_fam_path = RESULTS_DIR / "dsrs_brappa_family_rankings.csv"
out_nes_path = RESULTS_DIR / "dsrs_brappa_panel_nes_matrix.csv"

df_all_rankings.to_csv(out_rank_path, index=False)
df_all_family.to_csv(out_fam_path, index=False)
df_projected_nes.to_csv(out_nes_path)

print(f"\n[Saved] DSRS Stressor Rankings -> {out_rank_path}")
print(f"[Saved] DSRS Family Rankings -> {out_fam_path}")
print(f"[Saved] DSRS Projected NES Matrix -> {out_nes_path}")
print("=" * 70)
print("17_run_dsrs_classifier: Completed successfully!")
print("=" * 70)
