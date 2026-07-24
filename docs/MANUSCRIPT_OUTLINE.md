# Manuscript outline — *Floral scent under space radiation in Brassica rapa*

**Status:** draft (2026-07). Introduction, Methods, Results and Discussion are full prose;
Abstract, figure/table legends and reference list are drafted. Remaining placeholders are marked
**[brackets]** — mostly citations to add ([ref]) and specifics only the authors can supply
(tool versions, reference build, radiation OSDR accession, Ensembl release/date, Zenodo DOI).
Framing is deliberately modest and honest: a **resource + analytical framework + a statistically
bounded (non-significant) result** at 40 cGy — not a claimed discovery. The set-level test is
well powered (~3–4% coordinated shift detectable); a positive finding, if any, awaits the
collaborators' higher-power counts.

## Working titles (pick/adapt)
1. *A transcriptomic framework for testing whether space radiation alters floral scent, applied to* Brassica rapa
2. *Floral-scent biosynthesis genes and their response to low-dose radiation in* Brassica rapa*: a cross-species-anchored analysis*
3. *Is floral scent altered in space? Building and testing the question in* Brassica rapa

## Article type
Resource / methods-and-preliminary-results (e.g. *Frontiers in Plant Science / Astrobotany*,
*Life (Basel)*, *npj Microgravity*, or a data-descriptor venue). Not a high-claim discovery paper.

## Authors
[R. Barker (lead)]; [radiation-experiment collaborators — confirm contribution/authorship];
[others]. ORCID 0000-0001-5681-9857.

## Abstract (structured skeleton)
- **Background:** Lunar/Martian agriculture needs to know which plant traits change under
  galactic cosmic radiation (GCR); floral scent underlies pollination and is biochemically
  well-mapped, but its radiation response is untested.
- **Approach:** We combine two *B. rapa* (Wisconsin Fast Plant) transcriptomes — a High/Low
  floral-scent contrast and a 0/40 cGy radiation dose-response — reconcile their gene IDs,
  define a tiered floral-volatile gene set, and test whether radiation perturbs it, then place
  the result in a cross-species conservation framework.
- **Results:** [scent set of 363 genes across 4 routes, validated]; at 40 cGy there is **no
  statistically significant effect** on the scent gene set after correlation-aware and
  multiple-testing correction. The core volatile enzymes are unaffected; the only sub-threshold
  trend (permutation p≈0.06, BH q≈0.16) is a mild relative under-response of the
  ester/methyltransferase "tailoring" route — the same family that is lineage-specifically
  expanded and scent-defining across flowering plants.
- **Conclusions:** We provide a reusable framework and a specific, falsifiable prediction — if
  radiation alters scent, it acts on the tailoring methyltransferases, not the core synthases —
  which the current low-dose data cannot confirm and a higher-powered dataset should test.

## 1. Introduction

Sustained human presence beyond low-Earth orbit will depend on bioregenerative life-support
systems in which plants recycle atmosphere, water and waste and provide fresh food [ref]. Realising
this on the Moon or Mars requires understanding how the space environment — in particular chronic
exposure to ionising radiation from galactic cosmic rays (GCR), against which planetary surfaces
offer little shielding — reshapes plant biology [ref]. NASA's GeneLab / Open Science Data Repository
(OSDR) has catalogued transcriptomic and phenotypic responses of plants to spaceflight and to
ground-based radiation analogues [ref], but this work has concentrated on growth, development and
canonical stress pathways; reproductive traits, and pollination-related traits in particular, remain
comparatively unexamined.

Floral scent is an attractive trait through which to ask whether radiation perturbs reproduction.
The volatile organic compounds (VOCs) that constitute scent mediate plant–pollinator communication
and thereby fruit and seed set [ref], a consideration that becomes acute for the pollinator-dependent
crops likely to be grown in controlled off-Earth agriculture. Floral scent is also one of the
best-characterised outputs of plant specialised metabolism: its VOCs derive from a small number of
well-mapped routes — terpenoid, benzenoid/phenylpropanoid, fatty-acid–derived/green-leaf-volatile,
and apocarotenoid — whose core biosynthetic enzymes (terpene synthases, phenylalanine ammonia-lyase,
lipoxygenases, carotenoid-cleavage dioxygenases) are conserved across flowering plants, while the
final "tailoring" steps that define species-specific scent are carried out by rapidly evolving
methyltransferase and acyltransferase families (Dudareva et al. 2013; Pichersky & Gershenzon 2002).
This layered architecture — a conserved core with a labile decorating layer — makes scent both
biochemically tractable and well suited to asking *where* in a pathway an environmental perturbation
would act.

Whether ionising or space radiation alters floral scent is essentially untested. Addressing it
requires linking a radiation dose-response to a scent-focused gene framework and testing, with
adequate statistical power, whether scent-biosynthesis genes are preferentially affected — and,
if so, whether the effect falls on the conserved core or the evolutionarily labile tailoring layer.

Here we use *Brassica rapa* (Wisconsin Fast Plant), a fast-cycling Brassicaceae model, for which two
complementary bulk RNA-seq datasets are available: a High- versus Low-scent contrast and a 0/40 cGy
radiation dose-response. We reconcile the two datasets, curate a tiered floral-volatile gene set,
test radiation perturbation of that set with a correlation-aware, power-characterised gene-set
method, and interpret the result within a cross-species conservation framework spanning Brassicaceae,
Solanaceae and the classic scent-model species. Our aim is deliberately modest: a reusable resource
and analytical framework, an honest and statistically bounded first answer, and a specific,
falsifiable prediction to guide properly powered follow-up.

## 2. Materials and Methods

All analysis code is available at github.com/dr-richard-barker/B_rappa_LLGCSS (archived at Zenodo,
[DOI]); the scripts named below reside in the repository's `annotation/` directory, and each result
regenerates from the committed input tables by re-running the corresponding script.

### 2.1 Plant RNA-seq datasets
Two independent *Brassica rapa* (Wisconsin Fast Plant) bulk RNA-seq datasets were analysed. The
**floral-scent dataset** comprised eight libraries contrasting High- and Low-scent lines (four
biological replicates each), retrieved from the NCBI Sequence Read Archive (accessions
`SRR4417237`–`SRR4417244`) and processed through NASA Open Science Data Repository (OSDR) /
GeneLab-style tooling. The **radiation dataset** comprised 39 GeneLab-style libraries in a
three-factor design — radiation dose (0 or 40 cGy of simulated galactic cosmic radiation),
genotype (wild-type or *anthocyaninless*) and tissue preservative (DRS or RL) — obtained from NASA
OSDR ([accession/DOI to be added]). The radiation experiment is led by collaborators; the counts and
differential-expression results analysed here are an in-house re-analysis, with OSDR as the source of
record.

### 2.2 Read processing and differential expression
Reads were processed with a GeneLab RNA-seq workflow (quality control with FastQC/MultiQC, adapter
and quality trimming with Trim Galore!, alignment with STAR and transcript quantification with RSEM;
[confirm exact tool versions and *B. rapa* reference build for each dataset before submission]).
Differential expression was computed in iDEP/DESeq2 [ref] under a negative-binomial generalised
linear model. For the scent dataset, genes were tested between High- and Low-scent groups (reported
as log2 fold-change and Benjamini–Hochberg-adjusted *p*-value). For the radiation dataset, a combined
genotype×dose factor plus preservative was modelled (`~ Treatment_genotype + preservative`) to yield
the four pairwise contrasts, and a separate model with explicit two-way interaction terms
(dose×preservative, dose×genotype, preservative×genotype) provided the marginal 40-vs-0 cGy dose
effect and the interaction estimates used below. Per-gene log-normalised expression values for all 39
radiation libraries were exported from iDEP for the gene-set analysis (§2.5).

### 2.3 Gene-identifier reconciliation
The scent dataset uses Chiifu/Ensembl Plants `Bra######` gene identifiers (assembly `Brapa_1.0`) and
the radiation dataset uses `BRA######` identifiers. These denote the same loci differing only in
letter case; after case-normalisation all 31,756 scent genes mapped one-to-one onto radiation genes,
allowing the two datasets to be joined directly on gene identifier without cross-assembly mapping
(`build_annotation.py`).

### 2.4 Functional annotation and floral-scent gene-set definition
Gene Ontology (GO) terms and Pfam protein-domain assignments for every *B. rapa* gene were retrieved
from Ensembl Plants BioMart (schema `plants_mart`, dataset `brapa_eg_gene`, assembly `Brapa_1.0`;
release [xx], accessed [date]) via `build_annotation.py`. A floral-volatile biosynthesis gene set was
then curated in *B. rapa* space (`curate_scent_geneset.py`), because this Ensembl release exposes no
*Arabidopsis* orthologues or functional descriptions for *B. rapa*. Genes were assigned to
biosynthetic routes and to two confidence tiers on the basis of specific GO terms and Pfam domains,
following canonical accounts of plant volatile biosynthesis [Dudareva et al. 2013; Chen et al. 2011;
D'Auria 2006; Effmert et al. 2005]. **Tier 1** (core volatile-forming enzyme families) comprised
terpene synthases (Pfam PF01397/PF03936), SABATH methyltransferases (PF03492), lipoxygenases
(PF00305) and carotenoid-cleavage dioxygenases (PF03055). **Tier 2** (supporting/route-level genes)
comprised BAHD acyltransferases (PF02458), O-methyltransferases (PF00891; GO:0008171), phenylalanine
ammonia-lyases (PF00221) and genes carrying route-level GO terms for terpenoid, phenylpropanoid,
jasmonate and apocarotenoid metabolism; the full term/domain lists are defined in the script.
Deliberately broad terms (e.g. GO:0008299, "isoprenoid biosynthetic process," which also captures
sterol and photosynthetic-carotenoid genes) were excluded. Because a gene may satisfy criteria for
more than one route, a gene can appear under multiple routes. The set was validated on the scent axis
by testing, per route, the proportion of member genes differentially expressed between High- and
Low-scent lines (adjP < 0.1).

### 2.5 Radiation gene-set testing
Whether radiation preferentially perturbs scent genes was tested with a competitive gene-set
approach on the 40-vs-0 cGy effect (`phase3_robust.py`). For each gene a confound-adjusted dose
effect was computed from the 39-library log-normalised matrix as the mean, across the four
genotype×preservative strata, of the difference between the mean expression at 40 cGy and at 0 cGy
within each stratum. For a given gene set, the observed statistic was the difference between the mean
dose effect of set members and that of all other genes (background). Significance was assessed two
ways: (i) an analytic two-sided Mann–Whitney *U* test of set versus background dose effects; and
(ii) a **correlation-aware permutation test** in which the dose labels were shuffled within each
genotype×preservative stratum (5,000 permutations), the per-gene dose effect and the set-versus-
background statistic recomputed each time, and a two-sided empirical *p*-value obtained. Permuting
samples rather than genes preserves inter-gene correlation, which the analytic test ignores and which
otherwise inflates significance for co-expressed gene sets. Tests were run for the whole set, for
Tier-1 core genes, and for each biosynthetic route; permutation *p*-values across these seven tests
were corrected by the Benjamini–Hochberg procedure (reported as *q*). An earlier, unadjusted version
of the test (`phase3_geneset_test.py`) is retained for comparison.

### 2.6 Confound and interaction analysis
To assess whether any radiation response of scent genes depended on the experimental confounds, the
number of scent-set genes with a significant dose×genotype or dose×preservative interaction term
(adjP < 0.1) was counted from the iDEP interaction model.

### 2.7 Cross-species orthology and conservation
Because Ensembl Plants Compara does not include the Chiifu assembly used here (only
`brassica_rapa_ro18`), cross-species conservation was assessed with an *Arabidopsis*-anchored panel of
twelve canonical floral-scent genes spanning the four routes (`scent_orthology.py`). Orthologues of
each *A. thaliana* gene across all plant genomes were retrieved from the Ensembl Plants Compara REST
homology endpoint (`rest.ensembl.org`, `compara=plants`, orthologues, condensed format), and ortholog
counts were tabulated per species to summarise conservation and lineage-specific copy-number
expansion.

### 2.8 Model-species scent genes and full ortholog protocol
The classic model scent species (*Petunia hybrida*, *Antirrhinum majus*, *Rosa* spp.) are absent from
Ensembl Plants Compara (and *A. majus* from OrthoDB), precluding programmatic orthology. A curated
table of their functionally characterised scent genes was assembled from the primary literature and
cross-referenced to the *B. rapa* routes (`scent_reference_species.tsv`); protein accessions for the
landmark enzymes were retrieved from NCBI via E-utilities (`phase4c_link.py`;
`scent_query_accessions.tsv`). A reproducible reciprocal-orthology protocol (proteome sources,
OrthoFinder/DIAMOND commands and the intersection with the scent set) is provided for execution on
external compute (`phase4c_protocol.md`).

### 2.9 Software and data availability
Analyses used Python 3 with NumPy, SciPy and Matplotlib; figures were produced by `make_figures.py`
and tables by `make_tables.py`. Sequencing data trace to SRA (`SRR4417237`–`SRR4417244`) and NASA
OSDR (radiation accession [to be added]); model-species protein accessions are listed in
`scent_query_accessions.tsv`. All code and derived tables are archived at Zenodo ([DOI]).

## 3. Results
- **3.1 A joined two-experiment resource.** The scent and radiation transcriptomes share the
  Chiifu locus namespace; 31,756 genes join 1:1. → *Fig 1* (design + join schematic).
- **3.2 A validated, tiered floral-scent gene set.** 363 genes (108 Tier-1 core), all four routes;
  validated on the scent axis (DE High vs Low across every route). → *Fig 2*, *Table 1*.
- **3.3 Radiation does not significantly perturb the scent set at 40 cGy — an informative null.**
  Stratified (genotype×preservative) dose contrast; competitive test with a correlation-aware
  permutation null + BH FDR. Core enzymes unaffected (p≥0.37); the whole set and the
  ester/methyltransferase route give the only trends (analytic MWU p≈0.004 but permutation p≈0.06,
  **BH q≈0.16 — n.s.**), illustrating ~10× correlation inflation of the naive test. A power/
  sensitivity analysis shows the set-level test is well-powered — it would detect a coordinated
  shift of ≥~0.04 log2 (~3%) at 80% power — and the observed effects fall below that, so any
  coordinated scent-set response at 40 cGy is bounded to <~3–4% (not merely undetected). Effect
  sizes tiny; confounds clean (0 significant interactions). Report as a bounded null + a
  pre-specified trend + candidates (`Bra013161`, `Bra028224`, `Bra039555`). → *Fig 3*, *Fig S1*,
  *Table 2*.
- **3.4 The affected route is the evolutionarily labile scent node.** All scent routes are deeply
  conserved (orthologs to tomato/rice), but the volatile-tailoring methyltransferases
  (COMT/SABATH/BAHD) are lineage-specifically expanded (COMT1 29 Brassica copies; BSMT1 13) and are
  the scent-defining enzymes in Petunia (PhBSMT), snapdragon (AmBAMT) and rose (RhOMT). Rose's
  RhNUDX1 shows the same volatile via a non-orthologous enzyme (convergence caveat). → *Fig 4*,
  *Table 3*.

## 4. Discussion

We set out to test whether GCR-range radiation perturbs floral-scent biosynthesis in *B. rapa* and,
if so, where in the pathway. Having reconciled a scent and a radiation transcriptome and defined a
validated, tiered scent gene set, we found **no statistically significant coordinated effect** of
40 cGy on the set after correlation-aware, multiple-testing correction. Crucially, this is an
*informative* null rather than a failure to detect: the gene-set test is well powered, able to
resolve a coordinated shift as small as ~3–4% in the set's dose response at 80% power, and the
observed effects fall below that bound. The core volatile *synthase* enzymes are entirely
unresponsive; the only sub-threshold trend is a mild relative under-response of the
ester/methyltransferase tailoring route.

That trend is noteworthy less for its (non-significant) magnitude than for its location. Three
independent lines of evidence converge on the same part of the pathway. First, the radiation
gene-set test — where any signal it carries sits in the ester/methyltransferase route. Second, the
comparative genomics: while every scent route is deeply conserved across flowering plants, the
volatile-tailoring methyltransferases are lineage-specifically expanded — most strikingly the
Brassicaceae expansion of the O-methyltransferase/SABATH families (e.g. *COMT1*, with many *Brassica*
paralogues where Solanaceae and monocots retain one). Third, the characterised scent genes of the
model species reinforce it: methyl-benzoate and related volatiles are set by SABATH/BAHD tailoring
enzymes in *Petunia* (PhBSMT), snapdragon (AmBAMT) and rose (RhOMT) alike. Across rosids and
asterids, scent identity is fixed at this decorating step, and it is precisely there that our data
place whatever weak signal exists.

A coherent interpretation follows from the layered architecture of the pathway. The conserved,
dosage-sensitive core biosynthetic enzymes appear buffered against a low radiation dose, whereas the
regulatorily flexible, lineage-specifically expanded tailoring enzymes are the compartment in which
a subtle perturbation — and, more speculatively, scent plasticity in general — would first become
visible. This yields a concrete, falsifiable prediction: **if radiation alters floral scent, it
should act on the tailoring methyltransferases and acyltransferases, not the core synthases** — a
hypothesis that a higher-powered dataset can directly confirm or refute.

For space biology, the immediate implication is measured. At a modest 40 cGy we find no evidence
that radiation broadly reprograms scent biosynthesis, which is mildly reassuring for the prospect of
maintaining pollinator-attractive crops in shielded off-Earth growth systems — but the bound is only
as strong as the dose, the single time point, and the transcript-level readout allow, and the
framework's value is that it transfers directly to higher doses, chronic exposure, spaceflight
datasets and other species.

Beyond the biology, the study makes several transferable methodological points. The two datasets
initially appeared un-joinable because of a gene-identifier namespace mismatch that proved to be
mere letter-case — a cautionary example of how superficial identifier differences can obscure
directly comparable data. Where an annotation database lacks orthologue mappings for a non-model
species, a scent gene set can still be defined natively from GO terms and Pfam domains. The ~10-fold
gap we observed between naïve and correlation-aware p-values underscores that competitive gene-set
tests on co-expressed sets must account for inter-gene correlation. And the rose *RhNUDX1* case —
geraniol produced by a Nudix hydrolase rather than the usual terpene synthase — shows that
orthology-only comparisons will miss functionally convergent scent genes, so orthology must be paired
with pathway/function curation.

**Limitations.** The radiation data are a single low dose (40 cGy), one time point, and an in-house
re-analysis of collaborator-led data; the readout is transcriptional, not the emitted VOCs
themselves. The gene-set test is well powered for *coordinated* shifts (~3–4%) but has limited power
to detect individual large-effect genes (only 13 genome-wide differentially expressed genes at
40 cGy), so isolated responsive scent genes could be missed. The competitive-test assumptions are
addressed by the permutation null but not eliminated. Ensembl Plants Compara represents *B. rapa*
only by the R-o-18 assembly rather than our Chiifu data, and the classic scent-model species
(*Petunia*, snapdragon, rose) are absent from it entirely, so their comparison rests on curated
literature rather than programmatic orthology.

**Future directions.** The decisive next step is to re-run the identical, pre-specified gene-set and
power analysis on the collaborators' higher-powered radiation counts and, ideally, across a dose
series, to convert the present bound into either a confirmed tailoring-layer effect or a stronger
null. Pairing transcriptomics with VOC metabolomics would connect gene expression to emitted scent
and close the trait-level gap. Finally, executing the reciprocal-orthology protocol against the
*Petunia*, snapdragon and rose proteomes would place the *B. rapa* tailoring genes into explicit
cross-species orthogroups and test the conservation argument directly.

## 5. Data & code availability
- Code + derived tables: GitHub `dr-richard-barker/B_rappa_LLGCSS`, archived at Zenodo **[DOI]**.
- Scent raw reads: SRA `SRR4417237–244` (via OSDR/GeneLab). Radiation: [OSDR accession/DOI].
- Model-species scent proteins: NCBI/UniProt accessions in `annotation/scent_query_accessions.tsv`.

## Figures & Tables (all derivable from `annotation/`; rendered figures in `annotation/figures/`)

### Figure legends
**Figure 1. Two-experiment design and gene-identifier reconciliation.** Schematic of the two
*B. rapa* RNA-seq datasets combined here — the floral-scent experiment (8 libraries, High- vs
Low-scent lines, Chiifu `Bra######` identifiers) and the radiation dose-response experiment (39
libraries; dose 0/40 cGy × genotype WT/*anthocyaninless* × preservative DRS/RL, `BRA######`
identifiers). The two identifier sets denote the same Chiifu loci differing only in case;
case-normalisation joins 31,756 genes 1:1, enabling the tiered scent gene set and the downstream
radiation and cross-species analyses.

**Figure 2. The curated floral-scent gene set and its validation.** (**A**) The 363-gene set by
biosynthetic route and confidence tier (Tier 1, core volatile-forming enzyme families; Tier 2,
supporting/route-level), from GO-term and Pfam-domain curation. (**B**) Scent-axis validation: the
percentage of each route's genes differentially expressed between High- and Low-scent lines
(adjP < 0.1); the number of testable genes per route is shown above each bar. Every route contains
scent-responsive genes, confirming the set captures floral-scent biology in this species.

**Figure 3. Radiation dose effect on the scent gene set, by biosynthetic route.** Box plots of the
genotype×preservative-stratified radiation dose effect (40 vs 0 cGy, log2) for each scent route
against the genome-wide background (grey; dashed line, background median). Boxes show median and
IQR, whiskers 1.5×IQR, outliers omitted. No route is significant after a correlation-aware
permutation test (5,000 sample-label shuffles stratified by genotype and preservative) with
Benjamini–Hochberg correction; the ester/methyltransferase "tailoring" route gives the strongest
but still non-significant trend (permutation *p* ≈ 0.07, *q* ≈ 0.16). Core volatile-synthase
enzymes (Tier 1) are unaffected.

**Figure 4. Cross-species conservation of the floral-scent pathway.** Heatmap of Ensembl Plants
Compara ortholog counts for a panel of 12 Arabidopsis-anchored scent genes (rows) across
representative plant genomes (columns); colour encodes log(1 + ortholog count), cell numbers are
raw counts. All four routes are deeply conserved (orthologs present through Solanaceae and
monocots), whereas volatile-tailoring enzymes expand lineage-specifically — note the Brassicaceae
expansion of *COMT1* (rows for *Brassica* spp.) and the terpene-synthase expansions in grape and
tomato.

**Figure S1. Sensitivity of the radiation gene-set test.** Statistical power (α = 0.05, two-sided)
of the competitive gene-set test to detect a coordinated additive shift δ (log2) in the scent set's
40-vs-0 cGy dose effect, estimated by shifting the stratified permutation null. Curves for the whole
scent set, the ester/methyltransferase route and the Tier-1 core; filled markers indicate δ at 80%
power (≈0.04, 0.06, 0.04 log2 respectively); dashed vertical lines show the observed |effect| for
each set. The observed effects lie below the 80%-power thresholds, so a coordinated scent-set
response larger than ~3–4% would have been detected; the null is therefore informative, not merely
under-powered (at the set level). Source: `phase3_power.py`.

### Table legends
**Table 1. Composition of the curated *Brassica rapa* floral-scent gene set.** The 363 genes are
assigned by GO term and Pfam domain to biosynthetic routes (terpenoid; benzenoid/phenylpropanoid;
fatty-acid/green-leaf-volatile; apocarotenoid) and to an ester/methyltransferase tailoring category,
and to a confidence tier (Tier 1, core volatile-forming enzyme families — terpene synthases, SABATH
methyltransferases, lipoxygenases, carotenoid-cleavage dioxygenases; Tier 2, supporting/route-level
genes). Columns give per-route Tier-1 and Tier-2 counts, route totals, and example families.
Because a gene may carry evidence for more than one route, route totals sum to more than 363. Source:
`scent_geneset.tsv` (`curate_scent_geneset.py`).

**Table 2. Scent-associated candidate genes and their radiation dose response.** *B. rapa*
scent-set genes that are differentially expressed between High- and Low-scent lines
(Benjamini–Hochberg adjP < 0.1), ranked by the absolute genotype×preservative-stratified radiation
dose effect (40 vs 0 cGy, log2). Columns: gene (Chiifu `Bra` ID), biosynthetic route, confidence
tier, stratified dose effect, scent High/Low log2 fold-change and adjP, and the Pfam/GO enzyme
evidence. **Note:** the scent gene set as a whole shows *no* significant radiation response after
correlation-aware, multiple-testing correction (Figure 3); the genes listed are individual
candidates for targeted follow-up (e.g. in higher-powered data), not a validated radiation-responsive
set. adjP, BH-adjusted *p*-value; GLV, green-leaf volatile; MT, methyltransferase; SABATH,
salicylic-acid/benzoic-acid/theobromine methyltransferase family.

**Table 3. Characterised floral-scent biosynthesis genes in model scent species absent from Ensembl
Plants.** Landmark, functionally characterised scent genes from *Petunia hybrida*, *Antirrhinum
majus* (snapdragon) and *Rosa* spp. — none represented in Ensembl Plants Compara (snapdragon is also
absent from OrthoDB) — curated from the primary literature and cross-referenced to the *B. rapa*
biosynthetic routes. NCBI/UniProt accessions were retrieved programmatically for the landmark
enzymes (Methods; `scent_query_accessions.tsv`); the remaining citations are indicative and should be
verified against the primary literature before submission. The table illustrates that the
volatile-tailoring methyltransferase/acyltransferase step is the scent-defining node across lineages
(PhBSMT, AmBAMT, RhOMT) and that functionally convergent genes (rose RhNUDX1, a Nudix hydrolase
producing geraniol) can lack orthologs elsewhere in the pathway — a caveat for orthology-only
comparisons. Source: `scent_reference_species.tsv` + `scent_query_accessions.tsv`.

## Key references (verified this project)
Dudareva 2013 New Phytol 198:16; Verdonk 2005 Plant Cell (PhODO1, PMID 15805488); Nagegowda 2008
Plant J (AmNES/LIS, PMID 18363779); Magnard 2015 Science (RhNUDX1, PMID 26138978); Chen 2011 Plant J
(TPS). [Add: GeneLab/OSDR pipeline; DESeq2; iDEP; Ensembl Plants; Brassica triplication.]
