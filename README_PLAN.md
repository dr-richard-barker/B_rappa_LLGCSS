# Project Plan — *Brassica rapa*: Lunar LEAF Cosmic Radiation & Redox Decoder

**Status:** Completed Analysis · **Drafted:** 2026-11-03  
**Repo:** `B_rappa_LLGCSS` · **License:** CC0-1.0  

---

## 1. Resolution of Legacy Practice Dataset

The exploratory pilot practice dataset (`SRR4417237`–`SRR4417244`) has been **fully excised and decommissioned** from this project. It is not part of the space radiation biology narrative.

The entire codebase, website, results, and manuscript now focus exclusively on the **39 GeneLab-style RNA-seq libraries** from the NASA Artemis candidate **Lunar LEAF** investigation in *Brassica rapa* (Wisconsin Fast Plants).

---

## 2. Experimental Structure: 3-Way Factorial Radiation Study

| Factor | Levels | Details |
|---|---|---|
| **Radiation Dose** | 0 cGy, 40 cGy GCR | Sham control vs. Galactic Cosmic Ray simulation at NSRL (Brookhaven) |
| **Genotype** | Wild-Type, *anthocyaninless* (*anl*) | Protective flavonoid-producing WT vs. unshielded mutant |
| **Preservative** | DRS, RNAlater (RL) | Deep Space Preservative (ambient dry storage) vs. wet standard |

---

## 3. Core Milestones Completed

1. **Orthology Harmonization:** DIAMOND Reciprocal Best Hit (RBH) alignment between *B. rapa* (`BRA######`) and *Arabidopsis* TAIR10 (`AT#G#####`) generating 18,770 1:1 ortholog pairs.
2. **DSRS Stress Recognition:** Cross-species single-cell seed atlas mapping identifying `radiation_GCR` (mean $\rho = 0.1504$) and acute gamma (`rad_100Gy_co60`, $\rho = 0.2545$, $P < 10^{-4}$) as top-ranking stress classes.
3. **41D CVAE ROS Decoder:** Quantified the "Antioxidant Shield Deficit":
   - WT buffers radiation ionization into secondary $H_2O_2$ (51.1% probability, 3.84 h nadir).
   - *anl* mutant suffers an immediate, unbuffered primary superoxide free-radical burst ($O_2^{\bullet-}$, 24.9% probability, 2.23 h burst).
4. **Radiation Resilience Index (RRI):** 40 cGy GCR precipitates acute resilience collapse ($0.762 \to 0.636$, $P = 4.2 \times 10^{-4}$), driven by housekeeping network collapse ($r = 0.481 \to 0.126$). Sham-irradiated *anl* mutants display constitutive depression ($\text{RRI} = 0.638$, $r = 0.079$).
5. **Preservative Parity:** Verified that Deep Space Preservative (DRS) achieves high-fidelity correlation ($r = 0.297$, $P < 10^{-15}$) and RRI delta trajectory parity with RNAlater.
6. **DNA Damage Response Activation:** WT Fast Plants mount robust induction of core homologous recombination and base excision repair loci (BRCA1 $+2.03$, PARP2 $+0.86$, RAD51 $+0.16$) and upstream flavonoid biosynthesis (PAL1, CHS, F3H), whereas *anl* mutants display blunted, negative interaction terms.
7. **Spatial Meristematic Vulnerability:** ggPlantMap projections show root apical meristems (13.8%) and central vascular steles (12.5%) as primary damage foci, with floral reproductive tissues contributing $>25\%$ of total predicted stress susceptibility.

---

## 4. Deliverables

- **Live GitHub Pages Web Application:** [https://dr-richard-barker.github.io/B_rappa_LLGCSS/](https://dr-richard-barker.github.io/B_rappa_LLGCSS/)
- **Publication Manuscript Outline:** [`docs/MANUSCRIPT_OUTLINE.md`](docs/MANUSCRIPT_OUTLINE.md)
- **Quantitative Benchmark Tables:** [`docs/TABLES.md`](docs/TABLES.md) and [`docs/results.html`](docs/results.html)
- **Technical Report:** [`docs/technical_report_brappa_decoders.md`](docs/technical_report_brappa_decoders.md)
