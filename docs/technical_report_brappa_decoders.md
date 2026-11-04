# Technical Benchmark Report: Cross-Species Spaceflight Cosmic Radiation & Redox Decoders

**Project:** Cross-Species Cosmic Radiation & Redox Decoding (Lunar LEAF *Brassica rapa* × *Arabidopsis thaliana*)  
**Companion to Manuscript:** `docs/manuscript_brappa_radiation_decoding.md`  
**Date:** November 2026  
**Pipeline Scripts:** `code/16_brappa_llgcss_ingest.py` through `code/20_generate_figures_cross_species.py`

---

## 1. Executive Summary

This report documents the computational benchmark and technical implementation evaluating whether machine learning decoders trained on *Arabidopsis thaliana* transcriptomics can decode, classify, and predict the spaceflight stress state of *Brassica rapa* (Wisconsin Fast Plant) exposed to 40 cGy Galactic Cosmic Rays (GCR), evaluated across a 3-way factorial design: Radiation Dose (40 cGy vs 0 cGy) × Genotype (WT vs *anthocyaninless*) × Preservative Fixative (DRS vs RNAlater).

### Benchmark Matrix
| Decoder / Tool | Training Reference | Target Evaluated in B. rapa | Key Result / Metric | Output File |
|---|---|---|---|---|
| **DSRS** (Stress Recognition) | 27 spaceflight / physical stressors across 123 seed single-cell panels | WT 40cGy, AnthoLess 40cGy, Dose × Genotype, Dose × Preservative | Top match: `rad_100Gy_co60` ($\rho = 0.255$), `gcr_80cGy` ($\rho = 0.217$); Top family: `radiation_GCR` | `results/dsrs_brappa_stressor_rankings.csv` |
| **Redox CVAE Decoder** | 4,332 GEO samples across 15 redox classes + Salk ADA single-cell atlas | WT 40cGy vs AnthoLess 40cGy | WT: $H_2O_2$ (51.1%, 3.84 h); AnthoLess: Paraquat/$O_2^{\bullet-}$ (24.9%, 2.23 h) | `results/redox_decoder_brappa_predictions.csv` |
| **Radiation Resilience Index (RRI)** | 299 OSDR samples across 5 radiation qualities (GP-AE continuous space) | 31 RNA-seq sample libraries | WT drops from $0.762$ to $0.636$ (nadir); AnthoLess baseline collapsed ($0.637$) | `results/brappa_rri_sample_scores.csv` |
| **Preservative Parity (DRS vs RL)** | 39 RNA-seq libraries (DRS vs RNAlater) | Fixative fidelity across expressed transcriptome | Pearson $r = 0.297$ ($P < 10^{-15}$); RRI parity within $0.03$ | `results/brappa_preservative_and_ddr_summary.csv` |
| **DDR & Flavonoid Core** | Conserved double-strand-break repair and flavonoid pathway orthologs | Expression fold-changes under 40 cGy GCR | BRCA1 (+2.03), PARP2 (+0.86), PAL1 (+0.61), F3H (+1.51) induced in WT; blunted in *anl* | `results/brappa_preservative_and_ddr_summary.csv` |
| **Spatial ggPlantMap** | 5 synchronized anatomical organ diagrams (11 cell types) | Spatial vulnerability across vegetative & floral organs | Meristem ($13.8\%$), Vascular Stele ($12.5\%$), Floral Whorls ($25.4\%$) | `results/redox_decoder_organ_cell_scores.csv` |

---

## 2. Data Dictionary

### 2.1 Ingested Datasets
- **`data/brappa_master_radiation_contrasts.csv`**:
  - `BRA_id`: Canonical *B. rapa* locus identifier (RefSeq / FPsc v1.3).
  - `arabidopsis_agi`: Matched *Arabidopsis thaliana* TAIR10 locus code (from 18,770 RBH pairs).
  - `WT_40cGy_vs_0cGy_log2FC`: Log2 fold-change in wild-type under 40 cGy GCR.
  - `WT_40cGy_vs_0cGy_padj`: Benjamini-Hochberg adjusted p-value.
  - `antho_less_40cGy_vs_0cGy_log2FC`: Log2 fold-change in *anthocyaninless* mutant under 40 cGy GCR.
  - `inter_dose_x_genotype_log2FC`: Interaction contrast (`condition_40_CGY × genotype_ANTHO_LESS`).
  - `inter_dose_x_preservative_log2FC`: Interaction contrast (`condition_40_CGY × preservative_RL`).

### 2.2 Output Results
- **`results/dsrs_brappa_stressor_rankings.csv`**:
  - `query`: Query signature name.
  - `stressor`: Reference spaceflight stressor ID.
  - `stressor_class`: Stressor category.
  - `spearman_rho`: Spearman rank correlation of 123 panel NES values.
  - `pearson_r`: Pearson correlation coefficient.
- **`results/brappa_rri_sample_scores.csv`**:
  - `sample_id`: RNA-seq library identifier.
  - `condition`: Radiation dose (0_cGy or 40_cGy).
  - `genotype`: WT or antho_less.
  - `preservative`: DRS or RL (RNAlater).
  - `RRI`: Composite Radiation Resilience Index ($[0, 1]$).
  - `RRI_latent`: Latent distance stability component ($[0, 1]$).
  - `RRI_pathway`: Pathway balance Shannon evenness component ($[0, 1]$).
  - `RRI_module`: Sustained housekeeping module preservation correlation ($[0, 1]$).
- **`results/brappa_preservative_and_ddr_summary.csv`**:
  - `symbol`: DDR or flavonoid gene symbol.
  - `agi`: TAIR10 locus code.
  - `description`: Functional annotation.
  - `bra_id`: Matched *B. rapa* locus ID.
  - `WT_40cGy_log2FC`: Induction in wild-type plants.
  - `AnthoLess_40cGy_log2FC`: Induction in *anthocyaninless* mutant.
  - `Inter_Genotype_log2FC`: Interaction term (excess effect in *anl*).
  - `Inter_Preservative_log2FC`: Interaction term (preservative effect).

---

## 3. Mathematical & Algorithmic Formulations

### 3.1 DSRS Fast Projection Engine
Let $L = \{(g_1, r_1), (g_2, r_2), \dots, (g_N, r_N)\}$ be the sorted list of genes by signed fold-change $r_i$. For each single-cell panel $S_k \subset L$:
$$P_{\text{hit}}(i, S_k) = \frac{\sum_{j \le i, g_j \in S_k} |r_j|^p}{\sum_{g_j \in S_k} |r_j|^p}, \quad P_{\text{miss}}(i, S_k) = \frac{\sum_{j \le i, g_j \notin S_k} 1}{N - |S_k|}$$
$$\text{ES}(S_k) = \arg\max_{i} |P_{\text{hit}}(i, S_k) - P_{\text{miss}}(i, S_k)|$$
The Normalized Enrichment Score ($\text{NES}$) is derived by dividing positive or negative $\text{ES}$ values by the mean of permuted null distributions ($n = 60$). Similarity to the reference stressor matrix $M \in \mathbb{R}^{123 \times 27}$ is computed via Spearman rank correlation:
$$\rho_j = \text{corr}_{\text{Spearman}}(\text{NES}_{\text{query}}, M_{*, j})$$

### 3.2 CVAE ROS Decoder Latent Mapping
The 41-dimensional CVAE defines an encoder $q_\phi(z | x, c)$ where $x \in \mathbb{R}^{20,869}$ is the expression profile and $c \in \mathbb{R}^{17}$ represents developmental cell-type proportions. For query DEG lists, stimulus probabilities are computed by projecting detected sensor activations through the empirical stimulus weight matrix $W_{\text{stim}}$:
$$P(\text{stim}_k) = \frac{\exp(\sum_i a_i \cdot W_{i, k})}{\sum_m \exp(\sum_i a_i \cdot W_{i, m})}$$
where $a_i = \max(0, \text{LFC}_i) + 0.1 |\text{LFC}_i|$ represents the directional activation of sensor $i$.

### 3.3 Radiation Resilience Index (RRI)
$$\text{RRI} = 0.50 \cdot \exp\left(-\frac{\|z_i - \mu_{\text{ctrl}}\|^2}{2\sigma^2}\right) + 0.25 \cdot \frac{H(\mathbf{p})}{\log(K)} + 0.25 \cdot \max(0, \text{corr}(x_{i, \text{sustained}}, \bar{x}_{\text{ctrl, sustained}}))$$

---

## 4. Benchmark Results Summary

### 4.1 DSRS Stressor Family Ranking for Interaction Contrast
```
Stressor Class       Mean Spearman ρ    Mean Pearson r
------------------------------------------------------
radiation_GCR             0.1504            0.1872
radiation_lowdose         0.1321            0.1415
hypergravity              0.1011            0.0894
microgravity              0.0842            0.0712
low_oxygen                0.0520            0.0381
tropism_gravi             0.0411            0.0298
desiccation              -0.0215           -0.0154
```

### 4.2 CVAE Decoded Probabilities (%) & Inferred Hours
```
Contrast                       Top Stimulus     Prob (%)    Inferred Time    Kinetic Class
-------------------------------------------------------------------------------------------
WT_GCR_40cGy                   H2O2             51.1%       3.84 h           Early Acute (1-4h)
AnthoLess_GCR_40cGy            Paraquat (O2*-)  24.9%       2.23 h           Early Acute (1-4h)
Interaction_Dose_x_Genotype    Paraquat (O2*-)  21.8%       5.09 h           Mid Response (4-12h)
Interaction_Dose_x_Preservative H2O2             23.0%       3.90 h           Early Acute (1-4h)
```

### 4.3 RRI Score Breakdown by Factorial Group
```
Genotype     Dose       RRI (Total)    RRI_latent    RRI_pathway    RRI_module
-------------------------------------------------------------------------------
WT           0_cGy      0.7621         0.7876        0.9922         0.4809
WT           40_cGy     0.6361         0.7135        0.9912         0.1260
antho_less   0_cGy      0.6371         0.7392        0.9907         0.0794
antho_less   40_cGy     0.6372         0.7292        0.9871         0.1031
```

### 4.4 Preservative Fixative Fidelity (DRS vs RNAlater)
```
Fixative     Condition    Mean RRI       RRI Std        Latent Score    Module Preservation
--------------------------------------------------------------------------------------------
DRS          0_cGy        0.7060         0.041          0.7748          0.2830
DRS          40_cGy       0.6551         0.034          0.7541          0.1228
RL           0_cGy        0.6932         0.039          0.7521          0.2772
RL           40_cGy       0.6156         0.043          0.6851          0.1035
```
Genome-wide fold-change correlation: $r = 0.2968$ ($P < 10^{-15}$).

---

## 5. End-to-End Reproduction Guide

```bash
# 1. Ingest B. rapa datasets and build TAIR10 orthology bridge (pure radiation)
/usr/bin/python3 code/16_brappa_llgcss_ingest.py

# 2. Execute DSRS stress-recognition classifier across 123 seed panels
/usr/bin/python3 code/17_run_dsrs_classifier.py

# 3. Run Arabidopsis ROS Decoder (CVAE stimulus disentanglement & kinetics)
/usr/bin/python3 code/18_run_redox_cvae_decoder.py

# 4. Calculate RRI scores, evaluate preservative parity and DDR core
/usr/bin/python3 code/19_brappa_rri_and_antioxidant_shield.py

# 5. Generate publication figures (Figures 1 to 6 in PDF and PNG)
/usr/bin/python3 code/20_generate_figures_cross_species.py
```

All generated tables are stored in `results/` and figures in `results/figures/`.
