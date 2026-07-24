# Manuscript outline — *Floral scent under space radiation in Brassica rapa*

**Status:** draft scaffold (2026-07). Grounded in the analyses in `annotation/`. Placeholders
in **[brackets]**. Framing is deliberately modest: this is a **resource + analytical
framework + a preliminary, testable signal**, not a claimed discovery — the 40 cGy radiation
effect is small and needs the collaborators' higher-power counts to confirm.

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
- Space agriculture & the Lunar LEAF context; why plant reproduction/pollination matters off-Earth.
- Floral scent VOCs: terpenoid, benzenoid/phenylpropanoid, fatty-acid/GLV, apocarotenoid routes;
  scent's ecological role; its biochemical tractability (Dudareva 2013).
- GCR / low-dose ionizing radiation effects on plant transcriptomes (GeneLab/OSDR background).
- *B. rapa* Fast Plants as a model; the two datasets available.
- **Gap & aim:** no study links a radiation dose-response to a scent-gene framework; we build and
  test that link, and ask whether any effect falls on conserved vs labile parts of the pathway.

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
- **3.3 Radiation does not significantly perturb the scent set at 40 cGy.** Stratified
  (genotype×preservative) dose contrast; competitive test with a correlation-aware permutation
  null + BH FDR. Core enzymes unaffected (p≥0.37); the whole set and the ester/methyltransferase
  route give the only trends (analytic MWU p≈0.004 but permutation p≈0.06, **BH q≈0.16 — n.s.**),
  illustrating ~10× correlation inflation of the naive test. Effect sizes tiny; confounds clean
  (0 significant interactions). Report as a non-significant, pre-specified trend + candidates
  (`Bra013161`, `Bra028224`, `Bra039555`). → *Fig 3*, *Table 2*.
- **3.4 The affected route is the evolutionarily labile scent node.** All scent routes are deeply
  conserved (orthologs to tomato/rice), but the volatile-tailoring methyltransferases
  (COMT/SABATH/BAHD) are lineage-specifically expanded (COMT1 29 Brassica copies; BSMT1 13) and are
  the scent-defining enzymes in Petunia (PhBSMT), snapdragon (AmBAMT) and rose (RhOMT). Rose's
  RhNUDX1 shows the same volatile via a non-orthologous enzyme (convergence caveat). → *Fig 4*,
  *Table 3*.

## 4. Discussion
- Convergence of three independent analyses (radiation test, Brassica expansion, cross-lineage
  characterised genes) on the methyltransferase tailoring layer.
- Interpretation: core biosynthesis is buffered; the regulatorily flexible decorating step is where
  a subtle radiation effect (and, plausibly, scent plasticity generally) would first appear.
- **Space-biology implication:** a specific, testable prediction for higher-dose / spaceflight data.
- **Limitations (state plainly):** single low dose (40 cGy); in-house re-analysis; small,
  uncorrected/borderline effect; competitive-test assumptions; Chiifu vs R-o-18/compara mismatch;
  key scent models (Petunia/snapdragon/rose) absent from Ensembl.
- **Future work:** rerun on collaborators' counts + higher doses; run the OrthoFinder protocol;
  metabolite (VOC) measurement to connect transcript to emitted scent.

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
