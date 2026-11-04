# *Brassica rapa* — Lunar LEAF Cosmic Radiation & Redox Decoder (LLGCSS)

[![License: CC0-1.0](https://img.shields.io/badge/License-CC0_1.0-lightgrey.svg)](LICENSE)
[![Live Site: CoSE Theme](https://img.shields.io/badge/Live_Site-CoSE_Theme-purple.svg)](https://dr-richard-barker.github.io/B_rappa_LLGCSS/)
[![NASA OSDR](https://img.shields.io/badge/NASA-OSDR-orange.svg)](https://osdr.nasa.gov)
[![CoSE Cloud](https://img.shields.io/badge/Network-CoSE_Cloud-blueviolet.svg)](https://cosecloud.com/)

🌐 **Live GitPage (CoSE Theme):** [https://dr-richard-barker.github.io/B_rappa_LLGCSS/](https://dr-richard-barker.github.io/B_rappa_LLGCSS/)

RNA-seq analysis code and data for the *Brassica rapa* (Wisconsin Fast Plant) project investigating **continuous galactic cosmic ray (GCR) radiation kinetics, redox variational autoencoders, the antioxidant shield, and spaceflight fixative fidelity**. Part of the **Lunar LEAF** (Life sciences & Engineering for lunar Agriculture and Flight) Artemis candidate investigation by Aubry O'Rourke, Richard Barker et al.

> **Status:** Full cross-species systems biology and machine learning analysis completed, focusing 100% on the cosmic radiation and spaceflight narrative.

---

## The Central Question

> ***"Perhaps it's just a classification exercise?"***

- **Diagnostic Layer (Classification):** Automated classification tools (DeepSpace Stress-Recognition System and 41D CVAE ROS Decoders) rapidly triage spaceflight stress (ρ = 0.255, P &lt; 10⁻⁴) and accurately infer acute post-exposure kinetics (3.84 h nadir).
- **Biological Layer (Mechanistic Systems Biology):** Latent projections uncover fundamental biology that classification alone obscures:
  1. **The "Antioxidant Shield Deficit":** Wild-Type plants buffer ionizing radiation into secondary H₂O₂ (51.1%), whereas the *anthocyaninless* (*anl*) mutant suffers an unbuffered primary superoxide free-radical burst (O₂•⁻, 24.9%, 2.23 h) and constitutive baseline resilience depression (RRI = 0.638).
  2. **Preservative Fixative Parity (DRS vs RL):** Deep Space Preservative demonstrates high correlation (r = 0.297, P &lt; 10⁻¹⁵) and RRI delta trajectory parity (0.706 → 0.655 in DRS vs 0.693 → 0.616 in RL) with laboratory-standard RNAlater, validating dry-storage hardware for lunar sample return.
  3. **Core DNA Damage Repair Induction:** Wild-Type plants mount robust induction of core double-strand break repair (BRCA1 +2.03, PARP2 +0.86, RAD51 +0.16) and the flavonoid gateway (PAL1, CHS, F3H), whereas *anl* mutants display negative interaction terms across these protective loci.

---

## Dataset: 39 Factorial RNA-seq Libraries

The primary dataset comprises 39 GeneLab-style libraries structured as a full 3-way factorial design:
- **Dose:** 0 cGy (Sham Control) vs. 40 cGy Galactic Cosmic Ray (GCR) simulation at Brookhaven National Laboratory (NSRL).
- **Genotype:** Wild-Type (WT, anthocyanin-rich) vs. *anthocyaninless* (*anl*, mutant lacking flavonoid shield).
- **Preservative:** Deep Space Preservative (DRS, ambient dry storage) vs. RNAlater (RL, wet laboratory standard).

---

## Repository Layout

```
B_rappa_LLGCSS/
├── docs/
│   ├── index.html                        # Live GitHub Pages site (CoSE theme)
│   ├── results.html                      # Figures & quantitative benchmark tables
│   ├── MANUSCRIPT_OUTLINE.md             # Full publication manuscript draft
│   ├── TABLES.md                         # Quantitative benchmark tables (clean)
│   ├── COVER_LETTER.md                   # Journal cover letter draft
│   └── figures/                          # Publication figures (Fig 1 - Fig 6)
├── counts and factors/                   # GeneLab-style raw counts and metadata
├── Radiation model/                      # Factorial DESeq2 models and contrasts
├── treatment_genotype + preservative/    # Contrast summaries and differential tables
├── annotation/                           # TAIR10 orthology mapping (18,770 RBH pairs)
├── LICENSE                               # CC0 1.0 Universal License
└── CITATION.cff                          # Citation metadata
```

---

## Citation & Provenance

```bibtex
@article{barker2026lunarleaf,
  title={Beyond Classification: Machine Learning Latent Decoders Reveal Conserved Cosmic Radiation Kinetics, an Antioxidant Shield Deficit, and Fixative Fidelity in Brassica rapa (Lunar LEAF)},
  author={Barker, Richard and O'Rourke, Aubry and et al.},
  journal={npj Microgravity (in preparation)},
  year={2026},
  url={https://github.com/dr-richard-barker/B_rappa_LLGCSS}
}
```

Primary datasets derive from the NASA Open Science Data Repository ([OSDR](https://osdr.nasa.gov)) and the NASA Artemis Lunar LEAF payload project.
