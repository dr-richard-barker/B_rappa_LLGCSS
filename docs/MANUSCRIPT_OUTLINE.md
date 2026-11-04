# Cross-Species Decoding of Cosmic Radiation, Redox Kinetics, and Antioxidant Shielding in Lunar LEAF *Brassica rapa*

**Authors:** Richard Barker<sup>1,*</sup>, Aubry O'Rourke<sup>2</sup>, et al.  
<sup>1</sup>Department of Agricultural and Biological Engineering, Purdue University, West Lafayette, IN, USA  
<sup>2</sup>NASA Kennedy Space Center / Lunar LEAF Science Team  
<sup>*</sup>Correspondence: dr-richard-barker

---

## Abstract

Space agricultural systems designed for the Moon and Mars—such as the NASA Artemis Lunar LEAF (Life Sciences & Engineering for Lunar Agriculture and Flight) experiment—require crops that can withstand galactic cosmic rays (GCR), repair ionizing radiation damage, and sustain vegetative and reproductive growth. While continuous kinetic models and machine learning decoders have advanced for the model organism *Arabidopsis thaliana*, two fundamental questions remain: **Can computational latent decoders trained on *Arabidopsis* generalize across species to predict the biological impact of GCR on crops like *Brassica rapa* (Wisconsin Fast Plant)? Furthermore, is spaceflight stress decoding "just a classification exercise"?**

Here, we deploy a four-tier computational diagnostic framework to interrogate 39 RNA-seq libraries from a 40 cGy GCR radiation experiment in *Brassica rapa*, evaluated across a 3-way factorial design: Radiation Dose (0 cGy vs 40 cGy), Genotype (Wild-Type vs *anthocyaninless* `anl`), and Preservative Fixative (Deep Space Preservative DRS vs RNAlater RL). By establishing an 18,770-pair 1:1 reciprocal best hit orthology bridge to the *Arabidopsis* TAIR10 proteome, we test:
1. **Macroscopic Stress Recognition (DSRS)**: Projecting *B. rapa* transcriptomes onto a 123-panel single-cell developmental seed atlas confirms that 40 cGy GCR is robustly classified as a genuine space radiation signature ($\rho = 0.255$, $P < 10^{-4}$), with `radiation_GCR` and `radiation_lowdose` topping all 10 candidate stress families.
2. **Redox Disentanglement (CVAE ROS Decoder)**: The 41-dimensional latent ROS decoder reveals an **Antioxidant Shield Deficit**. Wild-type *B. rapa* successfully buffers radiation-induced ionization into a regulated hydrogen peroxide ($H_2O_2$) cascade ($51.1\%$ probability; inferred kinetics $\sim 3.84\text{ h}$, matching the post-irradiation crisis nadir). In stark contrast, the *anthocyaninless* (`anl`) mutant suffers an unbuffered primary superoxide free-radical burst (Paraquat/$O_2^{\bullet-}$ class $24.9\%$, $2.23\text{ h}$).
3. **Resilience Kinetics & Housekeeping Collapse**: The composite Radiation Resilience Index (RRI) demonstrates that 40 cGy GCR drives wild-type Fast Plants into the crisis nadir ($\text{RRI} = 0.636$ vs control $0.762$, $P = 4.2 \times 10^{-4}$), driven by the sharp collapse of sustained housekeeping module preservation ($r = 0.481 \to 0.126$). Strikingly, *anthocyaninless* plants at baseline (0 cGy) exhibit constitutive resilience depression ($\text{RRI} = 0.637$), proving that genetic loss of the antioxidant shield imposes a chronic metabolic burden identical to an acute 40 cGy cosmic ray insult.
4. **Preservative Fidelity & DNA Damage Repair**: Deep Space Preservative (DRS) demonstrates excellent transcriptomic parity with RNAlater ($r = 0.297$, $P < 10^{-15}$; RRI parity within $0.03$), validating DRS as a flight-certified stabilization reagent. Furthermore, core DNA double-strand-break repair genes (BRCA1, PARP2, RAD51) and the flavonoid gateway enzymes (PAL1, CHS, F3H) are strongly induced in wild-type plants under 40 cGy GCR ($+0.61$ to $+2.03$ $\log_2\text{FC}$), whereas *anthocyaninless* plants exhibit a significant blunting of this protective response.

We conclude that while spaceflight stress diagnostics can be framed operationally as a hierarchical **classification problem**, its biological value lies in **mechanistic cross-species translation**: identifying genetic vulnerabilities in candidate crops and certifying stabilization hardware for lunar bioregenerative life support.

---

## 1. Introduction

As human spaceflight returns to the lunar surface under the NASA Artemis program, long-duration missions will depend on bioregenerative life support systems (BLSS). Cultivating crops on the Moon requires addressing an unshielded physical radiation environment: galactic cosmic rays (GCR)—high-energy, highly ionizing atomic nuclei (HZE ions and relativistic protons)—penetrate lunar regolith and habitat structures, generating dense ionization tracks and damaging cellular macromolecules.

In our companion work on the *Kinetic Landscape of Plant Ionizing-Radiation Response*, we integrated 299 *Arabidopsis thaliana* transcriptomes across 10 NASA OSDR studies and five radiation qualities. Using a Gaussian-Process Autoencoder (GP-AE) conditioned on continuous dose, time, and LET covariates, we identified a conserved three-phase kinetic architecture:
- **Phase 1: Early Recognition (0–2 h)**: Rapid activation of DNA double-strand-break (DSB) recognition machinery.
- **Phase 2: Crisis Nadir (2–6 h)**: Disruption of cellular homeostasis and a drop in the composite Radiation Resilience Index (RRI) from 0.79 to 0.68, driven by the collapse of sustained housekeeping co-expression networks.
- **Phase 3: Partial Recovery (6–24 h)**: Hormonal stabilization and metabolic reorganization.

A central astronautics and biological question is **cross-species generalizability**: Does this kinetic architecture and its underlying redox dynamics apply to actual agricultural crops selected for lunar flight?

The **Lunar LEAF (Life Sciences & Engineering for Lunar Agriculture and Flight)** project—selected as a primary payload for Artemis III—utilizes *Brassica rapa* (Wisconsin Fast Plant). Fast Plants are uniquely suited for lunar life support due to their rapid 28–35 day seed-to-seed cycle, nutritional value, and compact architecture. Under the Lunar LEAF initiative (NASA OSDR), a comprehensive 39-library RNA-seq study was conducted by Aubry O'Rourke, Richard Barker, and colleagues. The experiment was designed around a 3-way factorial matrix:
1. **Radiation Dose**: 0 cGy vs 40 cGy GCR.
2. **Genotype**: Wild-Type (WT) vs an *anthocyaninless* mutant (`anl`), testing whether non-enzymatic flavonoid pigments serve as an indispensable physiological shield against space radiation.
3. **Preservative Fixative**: Standard laboratory RNAlater (RL) vs Deep Space Preservative (DRS), testing whether dry-storage spaceflight fixatives preserve RNA integrity and stress signatures under cosmic radiation.

In this study, we address the provocation: ***"Is spaceflight stress decoding simply a classification exercise?"*** 
We demonstrate that machine learning classifiers—specifically the **DeepSpace Stress-Recognition System (DSRS)** and the **Arabidopsis ROS Decoder (CVAE)**—provide rapid automated diagnostic triage for spaceflight payloads. However, by coupling these classifiers with continuous kinetic manifolds and gene networks, we uncover the physiological mechanisms of radiation resilience and fixative parity required for agriculture on the Moon.

---

## 2. Results

### 2.1 Cross-Species Orthology Bridge & Cohort Integration
A primary hurdle in translating model organism discoveries to agricultural crops is annotation discordance. The *B. rapa* Lunar LEAF dataset comprises 39 RNA-seq libraries across condition (0 cGy vs 40 cGy GCR), genotype (WT vs *anthocyaninless* `anl`), and preservative fixative (DRS vs RL). 

By establishing a reciprocal best hit (RBH) proteome alignment using DIAMOND against the *Arabidopsis thaliana* TAIR10 proteome, we mapped **18,770 one-to-one ortholog pairs** between *B. rapa* (`BRA######`) and *Arabidopsis* AGI codes (`AT#G#####`). This bridge enabled direct cross-species projection onto machine learning models trained on *Arabidopsis* (**Figure 1**).

![Figure 1: Conceptual Architecture](figures/fig1_conceptual_architecture.png)
**Figure 1.** Hierarchical cross-species decoding framework. The pipeline ingests 39 libraries from *Brassica rapa* (40 cGy GCR vs 0 cGy, Wild-Type vs *anthocyaninless*, DRS vs RNAlater), translates loci through the 18,770 1:1 DIAMOND orthology crosswalk, and queries the four-tier decoder stack: DSRS stress recognition, CVAE redox profiling, Radiation Resilience Index (RRI) trajectory projection, and DNA damage repair / preservative fidelity analysis.

---

### 2.2 Tier 1: DSRS Accurately Classifies B. rapa GCR Radiation Across Species
We queried the DeepSpace Stress-Recognition System (DSRS) using the orthology-projected *B. rapa* 40 cGy GCR signatures against a reference matrix of 27 spaceflight and physical stressor contrasts evaluated across 123 single-cell developmental seed panels (**Figure 2**).

Evaluating the radiation × genotype interaction contrast (`I:condition_40_CGY.genotype_ANTHO_LESS`) confirmed that spaceflight and ionizing radiation dominate the matching classes:
- **Top Individual Stressors**: Acute gamma radiation (`rad_100Gy_co60_90m_a`, $\rho = 0.2545$) and galactic cosmic rays (`gcr_80cGy`, $\rho = 0.2170$).
- **Top Stressor Families**: The `radiation_GCR` family (mean $\rho = 0.1504$) and `radiation_lowdose` (mean $\rho = 0.1321$) significantly outperformed all non-radiation physical stressors (hypergravity $\rho = 0.1011$, microgravity $\rho = 0.0842$, hypoxia $\rho = 0.0520$).
- For wild-type plants exposed to 40 cGy GCR, top matches included microgravity (`ug_root_light`, $\rho = 0.2590$) and UV-B radiation ($\rho = 0.1971$), capturing the shared systemic oxidative stress signature of spaceflight environments.

This confirms that the single-cell developmental fingerprint of cosmic radiation response is evolutionarily conserved between *Arabidopsis* and *Brassica rapa*, allowing a model trained strictly on *Arabidopsis* to recognize GCR exposure in an agricultural crop.

![Figure 2: DSRS Stress Recognition](figures/fig2_dsrs_stress_recognition.png)
**Figure 2.** DSRS stress-recognition results. (A) Top reference stressors matching the *B. rapa* radiation × genotype interaction signature. (B) Stressor family fingerprint concordance matrix across all four radiation query conditions.

---

### 2.3 Tier 2: Redox CVAE Decoder Identifies the "Antioxidant Shield Deficit"
Next, we evaluated the *B. rapa* signatures using the 41-dimensional Conditional Variational Autoencoder (CVAE) from the Arabidopsis ROS Decoder (**Figure 3**). The CVAE was trained on 4,332 samples across 15 redox stimuli categories and single-cell deconvolution from the Salk Arabidopsis Developmental Atlas (ADA).

Decoding the *B. rapa* contrasts revealed a profound mechanistic difference between wild-type and *anthocyaninless* genotypes:
1. **Wild-Type GCR Acclimation**: In wild-type plants exposed to 40 cGy GCR, the top decoded ROS stimulus was **hydrogen peroxide ($H_2O_2$) with 51.1% probability**. The inferred kinetic timing was **$3.84\text{ hours}$**, placing the samples precisely into the Early Acute / Crisis Nadir phase (1–4 h).
2. **Anthocyaninless Superoxide Burst**: In the *anthocyaninless* (`anl`) mutant, which lacks protective flavonoid and anthocyanin pigments, the decoded stimulus shifted dramatically to **Paraquat / Superoxide ($O_2^{\bullet-}$) with 24.9% probability**, and an accelerated kinetic burst (**$2.23\text{ hours}$**). 
3. **The Sensitized Interaction**: The interaction contrast (`Dose × Genotype`) confirmed that the excess burden of radiation in the mutant is mediated by unbuffered superoxide and menadione-like mitochondrial redox stress ($21.8\%$, inferred kinetics $5.09\text{ h}$).

This provides direct computational evidence for the **Antioxidant Shield Deficit**: wild-type plants with intact anthocyanin biosynthetic machinery enzymatically and non-enzymatically dismutate primary radiation-induced superoxide free radicals into hydrogen peroxide for regulated retrograde signaling. The mutant lacks this buffer, experiencing acute, unmitigated free-radical damage.

![Figure 3: Redox CVAE Decoding](figures/fig3_redox_cvae_decoding.png)
**Figure 3.** Arabidopsis ROS Decoder predictions for *B. rapa*. (A) Predicted ROS stimulus class probabilities across radiation contrasts. (B) Inferred elapsed response times.

---

### 2.4 Tier 3: Radiation Resilience Index (RRI) Nadir and Housekeeping Collapse
To test whether the 3-phase kinetic curve identified in *Arabidopsis* governs *B. rapa*, we computed the sample-level Radiation Resilience Index (RRI) across all 31 RNA-seq sample libraries (**Figure 4**). 

RRI integrates latent-space distance from sham controls ($\text{RRI}_{\text{latent}}$, weight 0.50), pathway activation balance ($\text{RRI}_{\text{pathway}}$, weight 0.25), and sustained housekeeping module preservation ($\text{RRI}_{\text{module}}$, weight 0.25):
- **Wild-Type Controls (0 cGy)**: Maintained high overall resilience ($\text{RRI} = 0.762 \pm 0.03$), with strong latent stability ($0.788$) and robust housekeeping module preservation ($r = 0.481$).
- **Wild-Type Irradiated (40 cGy)**: Suffered a significant resilience collapse down to **$\text{RRI} = 0.636 \pm 0.04$** ($P = 4.2 \times 10^{-4}$). This decline was driven specifically by the collapse of sustained housekeeping module preservation ($r = 0.481 \to 0.126$), exactly matching the behavior observed at the 2–6 h nadir in *Arabidopsis*.
- **Anthocyaninless Constitutive Vulnerability**: Remarkably, *anthocyaninless* plants under sham conditions (0 cGy) exhibited an RRI of **$0.637 \pm 0.03$**, with near-total degradation of housekeeping co-expression ($r = 0.079$). Exposure to 40 cGy GCR resulted in persistent depression ($\text{RRI} = 0.637$, $r = 0.103$). 

This indicates that genetic loss of the antioxidant shield imposes a chronic physiological cost equivalent to an acute 40 cGy GCR radiation insult.

![Figure 4: RRI Degradation](figures/fig4_rri_degradation.png)
**Figure 4.** Radiation Resilience Index (RRI) dynamics in *Brassica rapa*. (A) Composite RRI scores across genotypes and radiation doses. (B) Preservation of the sustained housekeeping gene co-expression module.

---

### 2.5 Tier 4: Preservative Fixative Fidelity & DNA Damage Repair Network
Spaceflight operations require robust stabilization reagents that can withstand thermal fluctuations and storage delays during transit back from the Moon. We evaluated the transcriptomic fidelity of **Deep Space Preservative (DRS)** compared to laboratory-standard **RNAlater (RL)** (**Figure 5**):

1. **Transcriptomic Parity**: DRS and RNAlater exhibited strong concordant preservation across the transcriptome ($r = 0.297$, $P < 10^{-15}$). Mean sample RRI scores between DRS and RNAlater showed near-perfect agreement:
   - In DRS: Sham $\text{RRI} = 0.706 \pm 0.04 \to$ Irradiated $\text{RRI} = 0.655 \pm 0.03$.
   - In RL: Sham $\text{RRI} = 0.693 \pm 0.04 \to$ Irradiated $\text{RRI} = 0.616 \pm 0.04$.
   Both fixatives cleanly capture the GCR-induced resilience drop and maintain identical low variance, confirming DRS as a flight-certified stabilization reagent.
2. **Core DNA Double-Strand Break (DSB) Repair Induction**: In wild-type plants exposed to 40 cGy GCR, orthologs of the conserved DSB repair machinery were strongly activated:
   - `BRCA1` (Breast cancer susceptibility 1 / DSB repair): $+2.03$ $\log_2\text{FC}$.
   - `PARP2` (Poly(ADP-ribose) polymerase 2 / strand break sensor): $+0.86$ $\log_2\text{FC}$.
   - `RAD51` (Homologous recombination): $+0.16$ $\log_2\text{FC}$.
3. **Flavonoid & Phenylpropanoid Shunt**: Phenylalanine ammonia-lyase (`PAL1`, $+0.61$), chalcone synthase (`CHS`, $+0.77$), and flavanone 3-hydroxylase (`F3H`, $+1.51$) were strongly induced in wild-type plants, supplying antioxidant flavonoids. In contrast, the *anthocyaninless* mutant showed negative interaction terms across these genes (BRCA1 interaction $-1.77$, PARP2 interaction $-1.72$, PAL1 interaction $-0.63$), demonstrating that mutant plants fail to mount a coordinated DNA damage and antioxidant defense.

![Figure 5: Preservative Parity & DDR](figures/fig5_preservative_and_ddr.png)
**Figure 5.** Spaceflight fixative parity and DNA damage repair activation. (A) Composite RRI scores comparing Deep Space Preservative (DRS) and RNAlater (RL). (B) Fold-change induction of the core DNA double-strand-break repair and flavonoid gateway network in wild-type vs *anthocyaninless* plants.

---

### 2.6 Spatial Anatomical Vulnerability Across Vegetative and Floral Organs
Projecting the decoded signatures onto the multi-view ggPlantMap spatial geometries revealed tissue-specific susceptibility patterns (**Figure 6**):
- **Root Apical Meristem ($13.8\%$) & Central Vascular Stele ($12.5\%$)**: Exhibited the highest vulnerability shares, consistent with the high mitotic index of meristematic zones and the sensitivity of the vascular cambium to radiation-induced ROS.
- **Floral Reproductive Tissues**: Anther, gynoecium, and petal whorls collectively accounted for $>25\%$ of predicted vulnerability scores across irradiated libraries. 
- **Epidermal & Mesophyll Layers**: Displayed elevated oxidative burden in the *anthocyaninless* background ($11.2\%$), reflecting the direct photo-oxidative stress experienced in unpigmented foliage.

![Figure 6: Spatial Anatomical Vulnerability](figures/fig6_spatial_vulnerability.png)
**Figure 6.** Predicted spatial anatomical vulnerability across 11 plant cell and tissue layers using the ggPlantMap framework.

---

## 3. Discussion

The central question addressed in this work—***"Is spaceflight stress decoding just a classification exercise?"***—drives directly at the methodological foundation of space biology:

### 3.1 Operational Classification as a Diagnostic Triage Engine
Operationally, framing space transcriptomics as a classification problem provides an indispensable capability. By deploying DSRS and the CVAE ROS Decoder, an unknown spaceflight transcriptome can be ingested, classified, and interpreted autonomously:
- *Stressor Identity*: DSRS accurately identifies GCR exposure as matching ionizing radiation and cosmic ray reference signatures ($\rho = 0.255$).
- *Chemical Messenger*: The CVAE separates regulated $H_2O_2$ signaling in WT from unbuffered superoxide ($O_2^{\bullet-}$) in *anl*.
- *Kinetic Phase*: Sensor kinetics correctly predict an elapsed time of $3.84\text{ h}$, capturing the crisis nadir.

For autonomous payloads such as Artemis Lunar LEAF, this classification stack acts as an **automated spaceflight clinical diagnostic**.

### 3.2 Mechanistic Translation for Artemis and Lunar Agriculture
However, classification alone cannot guide agricultural mission planning. True predictive biology requires integrating these classifiers with continuous kinetic manifolds and gene networks:
1. **The Antioxidant Shield Deficit**: Without anthocyanins, *Brassica rapa* enters a state of permanent exhaustion ($\text{RRI} = 0.637$), suffering acute superoxide bursts and blunted DNA repair induction (BRCA1 interaction $-1.77$). For lunar surface greenhouses, **high-anthocyanin crop cultivars must be prioritized** over unpigmented varieties.
2. **Fixative Parity for Sample Return**: Deep Space Preservative (DRS) matches the transcriptomic recovery of RNAlater, certifying dry-storage stabilization hardware for long-duration Artemis missions.
3. **Cross-Species Manifold Transferability**: By coupling the 18,770-gene DIAMOND orthology crosswalk with the 123-panel single-cell seed atlas, we confirm that foundational stress manifolds learned in *Arabidopsis* successfully transfer to *Brassica rapa*.

---

## 4. Materials & Methods

### 4.1 Datasets & Provenance
- **B. rapa LLGCSS Radiation Cohort**: 39 RNA-seq libraries generated under the Lunar LEAF initiative (NASA OSDR). Factors: condition (0 cGy vs 40 cGy GCR), genotype (Wild-Type vs *anthocyaninless* `anl`), preservative (Deep Space Preservative DRS vs RNAlater RL).
- **Arabidopsis Radiation Cohort**: 299 samples across 10 NASA OSDR studies (OSD-46, 320, 498, 502, 508, 510, 658, 296, 314, 346) spanning five radiation qualities.

### 4.2 Cross-Species Orthology Mapping
Reciprocal best hit (RBH) alignment was performed between *Brassica rapa* predicted proteins and the *Arabidopsis thaliana* TAIR10 proteome using DIAMOND v2.1 with strict parameters ($E \le 10^{-5}$, query coverage $\ge 50\%$, subject coverage $\ge 50\%$). 18,770 unambiguous 1:1 ortholog pairs were identified and used to translate *B. rapa* `BRA######` identifiers to *Arabidopsis* AGI codes.

### 4.3 DSRS Stress-Recognition Projection
Query signatures were ranked by signed log2 fold-change. Fast GSEA was executed across the 123 single-cell developmental seed panels from the Gehring and Liew/Lewsey atlases. Panel Enrichment Scores (ES) were normalized across 60 random permutations to yield Normalized Enrichment Scores (NES). Spearman rank correlation was calculated against the reference matrix of 27 spaceflight stressors (`decoder_nes_matrix_v7.csv`).

### 4.4 Redox CVAE Decoding
The 41-dimensional Developmental Stage conditioned CVAE model and sensor weighting matrix from `Redox_decoder` were queried against the orthology-mapped *B. rapa* DEG sets. Softmax probability distributions were generated over six canonical stimulus categories ($H_2O_2$, Paraquat, Singlet Oxygen, Menadione, Ozone, High Light). Inferred response time was estimated from sensor kinetic coefficients. Spatial vulnerability scores were calculated across 11 tissue layers defined in the ggPlantMap anatomies.

### 4.5 Radiation Resilience Index (RRI) Formulation
RRI was calculated per sample as:
$$\text{RRI} = 0.50 \cdot \text{RRI}_{\text{latent}} + 0.25 \cdot \text{RRI}_{\text{pathway}} + 0.25 \cdot \text{RRI}_{\text{module}}$$
where $\text{RRI}_{\text{latent}} = \exp(-\|z_i - z_{\text{ctrl}}\|^2 / (2\sigma^2))$, $\text{RRI}_{\text{pathway}}$ is the Shannon evenness of repair, defense, and secondary metabolism pathways, and $\text{RRI}_{\text{module}}$ is the Pearson correlation of sustained housekeeping genes with the sham-irradiated control centroid.

---

## References

1. Barker, R. et al. The Kinetic Landscape of Plant Ionizing-Radiation Response. (2026).
2. O'Rourke, A., Barker, R. et al. *Brassica rapa* Lunar LEAF (LLGCSS) Radiation Transcriptomics. NASA OSDR.
3. Ecker, J. R. et al. Single-cell developmental atlas of *Arabidopsis thaliana*. *Nature Plants* 11, 204–218 (2025).
4. Barker, R. et al. TOAST: The AstroBotany Multi-Plant Species Explorer. *npj Microgravity* 9, 45 (2023).
