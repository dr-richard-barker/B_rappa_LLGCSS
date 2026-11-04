# Data Files & Accessions: Brassica rapa Lunar LEAF Study

This repository tracks **reproducible code, orthology maps, and metadata** in git. Large sequencing data and count matrices are preserved on disk and accessible via NASA OSDR.

## Raw Sequencing & Study Data

| Dataset | Accessions | Source | Experimental Design |
|---|---|---|---|
| Lunar LEAF Radiation / GCR | NASA OSDR (collaborator-led) | NASA OSDR / GeneLab | 39 libraries: 0 vs 40 cGy GCR × WT vs *anl* × DRS vs RL preservative |

## Reference Genomes & Orthology

| Assembly | Annotation | Use |
|---|---|---|
| NCBI RefSeq `GCF_000309985.2` (CAAS Brap v3.01, Chiifu) | FPsc / RefSeq (`BRA######` / `Bra######`) | Reference genome for Lunar LEAF counts |
| *Arabidopsis thaliana* TAIR10 | TAIR10 (`AT#G#####`) | Orthology target for cross-species latent decoders |

## Key Tables

- `counts and factors/LEAF_Metadata.csv` — radiation sample sheet (canonical 39 libraries)
- `counts and factors/Results_LFC_Pval_DESeq2.csv` — DESeq2 contrast statistics
- `Radiation model/Radi_model_with_interactions_deg_values_DESeq2.csv` — factorial interaction model
- `annotation/brapa_to_arabidopsis_orthologs.tsv` — 18,770 1:1 DIAMOND RBH ortholog pairs
