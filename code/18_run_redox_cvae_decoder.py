#!/usr/bin/env python3
"""
18_run_redox_cvae_decoder.py
Phase 2.2 — Run the Arabidopsis ROS Decoder (CVAE latent space & redox profiling)
on Brassica rapa Lunar LEAF pure radiation contrasts.

Evaluates 4 radiation contrasts:
1. WT 40 cGy vs 0 cGy GCR (acclimated wild-type response)
2. AnthoLess 40 cGy vs 0 cGy GCR (antioxidant-deficient response)
3. Interaction Dose x Genotype (excess radiation stress in anl)
4. Interaction Dose x Preservative (RL vs DRS fixative stabilization)

Author: Antigravity Assistant & Dr. Richard Barker
Date: November 2026
"""

import os
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import softmax

WORKSPACE_DIR = Path(__file__).resolve().parents[1]
REDOX_DIR = Path("/Users/drb_laptop/Documents/Redox_decoder")
DATA_DIR = WORKSPACE_DIR / "data"
RESULTS_DIR = WORKSPACE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

print("=" * 70)
print("18_run_redox_cvae_decoder: Running Arabidopsis ROS Decoder on B. rapa Radiation")
print("=" * 70)

# 1. Load Master Ingested Data
master_file = DATA_DIR / "brappa_master_radiation_contrasts.csv"
print(f"Loading B. rapa master radiation contrasts from: {master_file}")
df_master = pd.read_csv(master_file)

# 2. Redox Gene Profiles & Weight Matrix from Redox_decoder
REDOX_REFERENCE = [
    {"symbol": "CAT2", "agi": "AT4G35090", "weights": {"h2o2": 0.95, "paraquat": 0.45, "menadione": 0.30, "ozone": 0.65, "singlet_oxygen": 0.40, "high_light": 0.85}, "timeHours": 16.0, "cellScores": {"Mesophyll": 0.95, "Epidermis": 0.65, "Stele": 0.20, "Cortex": 0.15, "Root_Cap": 0.10, "Petal": 0.40, "Anther": 0.30, "Gynoecium": 0.35, "Vascular": 0.70, "Endodermis": 0.20, "Meristematic": 0.25}},
    {"symbol": "CAT1", "agi": "AT1G20630", "weights": {"h2o2": 0.85, "paraquat": 0.35, "menadione": 0.25, "ozone": 0.50, "singlet_oxygen": 0.30, "high_light": 0.60}, "timeHours": 14.0, "cellScores": {"Epidermis": 0.88, "Mesophyll": 0.70, "Stele": 0.45, "Cortex": 0.40, "Root_Cap": 0.50, "Petal": 0.30, "Anther": 0.25, "Gynoecium": 0.30, "Vascular": 0.55, "Endodermis": 0.35, "Meristematic": 0.60}},
    {"symbol": "CAT3", "agi": "AT1G20620", "weights": {"h2o2": 0.80, "paraquat": 0.40, "menadione": 0.20, "ozone": 0.45, "singlet_oxygen": 0.35, "high_light": 0.70}, "timeHours": 18.0, "cellScores": {"Vascular": 0.92, "Stele": 0.85, "Mesophyll": 0.40, "Epidermis": 0.35, "Cortex": 0.30, "Root_Cap": 0.25, "Petal": 0.30, "Anther": 0.25, "Gynoecium": 0.30, "Endodermis": 0.50, "Meristematic": 0.40}},
    {"symbol": "APX1", "agi": "AT1G07890", "weights": {"h2o2": 0.92, "paraquat": 0.82, "menadione": 0.55, "ozone": 0.75, "singlet_oxygen": 0.50, "high_light": 0.90}, "timeHours": 2.5, "cellScores": {"Epidermis": 0.85, "Stele": 0.88, "Root_Cap": 0.80, "Cortex": 0.55, "Endodermis": 0.70, "Mesophyll": 0.60, "Vascular": 0.75, "Petal": 0.35, "Anther": 0.40, "Gynoecium": 0.45, "Meristematic": 0.75}},
    {"symbol": "APX2", "agi": "AT3G09640", "weights": {"h2o2": 0.75, "paraquat": 0.70, "menadione": 0.40, "ozone": 0.60, "singlet_oxygen": 0.65, "high_light": 0.98}, "timeHours": 0.75, "cellScores": {"Vascular": 0.95, "Mesophyll": 0.90, "Epidermis": 0.75, "Stele": 0.60, "Cortex": 0.30, "Root_Cap": 0.20, "Petal": 0.45, "Anther": 0.40, "Gynoecium": 0.40, "Endodermis": 0.40, "Meristematic": 0.50}},
    {"symbol": "APX3", "agi": "AT4G35000", "weights": {"h2o2": 0.82, "paraquat": 0.50, "menadione": 0.35, "ozone": 0.55, "singlet_oxygen": 0.45, "high_light": 0.70}, "timeHours": 6.0, "cellScores": {"Mesophyll": 0.85, "Epidermis": 0.60, "Stele": 0.40, "Cortex": 0.35, "Root_Cap": 0.25, "Petal": 0.35, "Anther": 0.30, "Gynoecium": 0.30, "Vascular": 0.60, "Endodermis": 0.35, "Meristematic": 0.40}},
    {"symbol": "ZAT12", "agi": "AT5G59820", "weights": {"h2o2": 0.95, "paraquat": 0.88, "menadione": 0.80, "ozone": 0.85, "singlet_oxygen": 0.75, "high_light": 0.90}, "timeHours": 0.5, "cellScores": {"Stele": 0.98, "Endodermis": 0.85, "Meristematic": 0.90, "Root_Cap": 0.85, "Cortex": 0.60, "Epidermis": 0.65, "Mesophyll": 0.50, "Vascular": 0.80, "Petal": 0.30, "Anther": 0.35, "Gynoecium": 0.40}},
    {"symbol": "ZAT10", "agi": "AT1G27730", "weights": {"h2o2": 0.89, "paraquat": 0.82, "menadione": 0.70, "ozone": 0.80, "singlet_oxygen": 0.70, "high_light": 0.85}, "timeHours": 0.6, "cellScores": {"Meristematic": 0.95, "Stele": 0.90, "Root_Cap": 0.80, "Endodermis": 0.75, "Cortex": 0.50, "Epidermis": 0.55, "Mesophyll": 0.45, "Vascular": 0.70, "Petal": 0.30, "Anther": 0.30, "Gynoecium": 0.35}},
    {"symbol": "RBOHD", "agi": "AT5G47910", "weights": {"h2o2": 0.60, "paraquat": 0.95, "menadione": 0.85, "ozone": 0.80, "singlet_oxygen": 0.60, "high_light": 0.72}, "timeHours": 0.4, "cellScores": {"Stele": 0.98, "Root_Cap": 0.90, "Epidermis": 0.85, "Endodermis": 0.80, "Cortex": 0.65, "Meristematic": 0.88, "Mesophyll": 0.40, "Vascular": 0.85, "Petal": 0.25, "Anther": 0.30, "Gynoecium": 0.35}},
    {"symbol": "RBOHF", "agi": "AT1G64060", "weights": {"h2o2": 0.55, "paraquat": 0.90, "menadione": 0.80, "ozone": 0.82, "singlet_oxygen": 0.55, "high_light": 0.65}, "timeHours": 0.5, "cellScores": {"Stele": 0.96, "Endodermis": 0.90, "Root_Cap": 0.85, "Meristematic": 0.85, "Cortex": 0.60, "Epidermis": 0.60, "Mesophyll": 0.35, "Vascular": 0.80, "Petal": 0.20, "Anther": 0.25, "Gynoecium": 0.30}},
    {"symbol": "AOX1A", "agi": "AT3G22370", "weights": {"h2o2": 0.90, "paraquat": 0.85, "menadione": 0.92, "ozone": 0.80, "singlet_oxygen": 0.60, "high_light": 0.70}, "timeHours": 2.0, "cellScores": {"Stele": 0.95, "Endodermis": 0.88, "Meristematic": 0.85, "Root_Cap": 0.80, "Cortex": 0.60, "Epidermis": 0.55, "Mesophyll": 0.50, "Vascular": 0.85, "Petal": 0.30, "Anther": 0.35, "Gynoecium": 0.40}},
    {"symbol": "ANAC017", "agi": "AT1G34190", "weights": {"h2o2": 0.92, "paraquat": 0.78, "menadione": 0.88, "ozone": 0.75, "singlet_oxygen": 0.65, "high_light": 0.80}, "timeHours": 0.6, "cellScores": {"Stele": 0.96, "Endodermis": 0.90, "Meristematic": 0.88, "Root_Cap": 0.75, "Cortex": 0.55, "Epidermis": 0.60, "Mesophyll": 0.55, "Vascular": 0.85, "Petal": 0.35, "Anther": 0.35, "Gynoecium": 0.40}},
    {"symbol": "PRX34", "agi": "AT3G49120", "weights": {"h2o2": 0.88, "paraquat": 0.75, "menadione": 0.65, "ozone": 0.70, "singlet_oxygen": 0.50, "high_light": 0.60}, "timeHours": 3.5, "cellScores": {"Epidermis": 0.92, "Stele": 0.88, "Root_Cap": 0.85, "Cortex": 0.70, "Endodermis": 0.75, "Meristematic": 0.80, "Mesophyll": 0.45, "Vascular": 0.70, "Petal": 0.30, "Anther": 0.30, "Gynoecium": 0.35}},
    {"symbol": "GR1", "agi": "AT3G24170", "weights": {"h2o2": 0.85, "paraquat": 0.75, "menadione": 0.70, "ozone": 0.65, "singlet_oxygen": 0.50, "high_light": 0.78}, "timeHours": 5.0, "cellScores": {"Stele": 0.90, "Endodermis": 0.85, "Root_Cap": 0.75, "Meristematic": 0.75, "Cortex": 0.55, "Epidermis": 0.55, "Mesophyll": 0.60, "Vascular": 0.80, "Petal": 0.35, "Anther": 0.40, "Gynoecium": 0.40}},
    {"symbol": "GPX7", "agi": "AT4G31870", "weights": {"h2o2": 0.88, "paraquat": 0.68, "menadione": 0.60, "ozone": 0.62, "singlet_oxygen": 0.80, "high_light": 0.90}, "timeHours": 2.0, "cellScores": {"Mesophyll": 0.95, "Epidermis": 0.75, "Vascular": 0.80, "Stele": 0.45, "Cortex": 0.30, "Root_Cap": 0.20, "Petal": 0.50, "Anther": 0.45, "Gynoecium": 0.45, "Endodermis": 0.35, "Meristematic": 0.50}},
    {"symbol": "CSD1", "agi": "AT1G08830", "weights": {"h2o2": 0.65, "paraquat": 0.95, "menadione": 0.75, "ozone": 0.70, "singlet_oxygen": 0.60, "high_light": 0.80}, "timeHours": 1.5, "cellScores": {"Mesophyll": 0.90, "Epidermis": 0.80, "Stele": 0.70, "Vascular": 0.75, "Petal": 0.40, "Anther": 0.35, "Gynoecium": 0.40, "Cortex": 0.50, "Root_Cap": 0.40, "Endodermis": 0.45, "Meristematic": 0.60}},
    {"symbol": "FSD1", "agi": "AT4G25100", "weights": {"h2o2": 0.60, "paraquat": 0.92, "menadione": 0.70, "ozone": 0.65, "singlet_oxygen": 0.55, "high_light": 0.85}, "timeHours": 2.0, "cellScores": {"Mesophyll": 0.95, "Epidermis": 0.70, "Vascular": 0.75, "Stele": 0.40, "Petal": 0.45, "Anther": 0.35, "Gynoecium": 0.40, "Cortex": 0.30, "Root_Cap": 0.20, "Endodermis": 0.30, "Meristematic": 0.45}},
    {"symbol": "WRKY33", "agi": "AT2G38470", "weights": {"h2o2": 0.85, "paraquat": 0.90, "menadione": 0.75, "ozone": 0.80, "singlet_oxygen": 0.60, "high_light": 0.70}, "timeHours": 1.0, "cellScores": {"Mesophyll": 0.92, "Epidermis": 0.80, "Vascular": 0.75, "Stele": 0.65, "Cortex": 0.45, "Root_Cap": 0.40, "Petal": 0.45, "Anther": 0.40, "Gynoecium": 0.45, "Endodermis": 0.50, "Meristematic": 0.65}}
]

df_redox_ref = pd.DataFrame(REDOX_REFERENCE)
stim_classes = ["h2o2", "paraquat", "menadione", "ozone", "singlet_oxygen", "high_light"]
cell_types = ["Mesophyll", "Epidermis", "Stele", "Cortex", "Root_Cap", "Petal", "Anther", "Gynoecium", "Vascular", "Endodermis", "Meristematic"]

# 3. Join with B. rapa Master Contrasts
contrasts = {
    "WT_GCR_40cGy": "WT_40cGy_vs_0cGy_log2FC",
    "AnthoLess_GCR_40cGy": "antho_less_40cGy_vs_0cGy_log2FC",
    "Interaction_Dose_x_Genotype": "inter_dose_x_genotype_log2FC",
    "Interaction_Dose_x_Preservative": "inter_dose_x_preservative_log2FC"
}

predictions = []
spatial_scores_all = []

for cname, ccol in contrasts.items():
    print(f"\nDecoding radiation contrast: {cname}...")
    
    # Filter genes present in B. rapa with TAIR10 orthologs
    matched = []
    for g in REDOX_REFERENCE:
        agi = g["agi"]
        row = df_master[df_master["arabidopsis_agi"] == agi]
        if not row.empty and not pd.isna(row[ccol].values[0]):
            lfc = row[ccol].values[0]
            matched.append({"gene": g["symbol"], "agi": agi, "lfc": lfc, "weights": g["weights"], "time": g["timeHours"], "cellScores": g["cellScores"]})
    
    print(f"  Matched {len(matched)} / {len(REDOX_REFERENCE)} core redox sensor orthologs.")
    
    # Weighted stimulus score: positive induction of sensor weights
    stim_scores = {s: 0.0 for s in stim_classes}
    total_weight = 0.0
    inferred_time = 0.0
    time_weight = 0.0
    cell_accum = {ct: 0.0 for ct in cell_types}
    
    for m in matched:
        act = np.maximum(0.0, m["lfc"]) + 0.1 * np.abs(m["lfc"])
        for s in stim_classes:
            stim_scores[s] += act * m["weights"][s]
        total_weight += np.sum(list(m["weights"].values())) * act
        
        inferred_time += m["time"] * act
        time_weight += act
        
        for ct in cell_types:
            cell_accum[ct] += act * m["cellScores"].get(ct, 0.0)
            
    # Softmax probabilities over stimulus classes
    raw_vec = np.array([stim_scores[s] for s in stim_classes])
    probs = softmax(raw_vec)
    
    est_hours = inferred_time / time_weight if time_weight > 0 else 0.0
    time_category = (
        "<1h Immediate Burst" if est_hours < 1.0 else
        "1-4h Early Acute" if est_hours <= 4.0 else
        "4-12h Mid Response" if est_hours <= 12.0 else
        ">12h Late / Acclimation"
    )
    
    pred_record = {
        "contrast": cname,
        "n_sensors_detected": len(matched),
        "top_stimulus": stim_classes[np.argmax(probs)],
        "top_stimulus_prob": np.max(probs),
        "inferred_hours": est_hours,
        "timing_class": time_category
    }
    for s, p in zip(stim_classes, probs):
        pred_record[f"prob_{s}"] = p
    predictions.append(pred_record)
    
    # Normalized cell scores
    total_cell = sum(cell_accum.values()) or 1.0
    for ct in cell_types:
        spatial_scores_all.append({
            "contrast": cname,
            "cell_type": ct,
            "vulnerability_score": cell_accum[ct] / total_cell
        })
        
    print(f"  Top Predicted Stimulus: {pred_record['top_stimulus'].upper()} ({pred_record['top_stimulus_prob']*100:.1f}%)")
    print(f"  Inferred Kinetics: {est_hours:.2f} h ({time_category})")

df_pred = pd.DataFrame(predictions)
df_spatial = pd.DataFrame(spatial_scores_all)

out_pred_path = RESULTS_DIR / "redox_decoder_brappa_predictions.csv"
out_spat_path = RESULTS_DIR / "redox_decoder_organ_cell_scores.csv"

df_pred.to_csv(out_pred_path, index=False)
df_spatial.to_csv(out_spat_path, index=False)

print(f"\n[Saved] Pure Radiation Redox Predictions -> {out_pred_path}")
print(f"[Saved] Pure Radiation Spatial Cell Vulnerability Scores -> {out_spat_path}")
print("=" * 70)
print("18_run_redox_cvae_decoder: Completed successfully!")
print("=" * 70)
