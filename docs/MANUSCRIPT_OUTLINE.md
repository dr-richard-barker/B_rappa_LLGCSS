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

## 2. Materials & Methods
- **2.1 Datasets.** Scent: SRA `SRR4417237–244`, High/Low, processed via NASA OSDR/GeneLab.
  Radiation: [collaborator/OSDR accession], 0/40 cGy × WT/*anthocyaninless* × DRS/RL preservative,
  39 libraries; in-house DESeq2/iDEP DE (main + interaction terms).
- **2.2 Gene-ID reconciliation.** Scent `Bra######` (Chiifu/Ensembl `Brapa_1.0`) vs radiation
  `BRA######` shown to be the same loci differing only in case (100% join); `build_annotation.py`.
- **2.3 Scent gene set.** GO + Pfam from Ensembl Plants BioMart; tiered curation
  (Tier-1 core enzyme families vs Tier-2 supporting), 4 routes; `curate_scent_geneset.py`. Sources:
  Dudareva 2013; Chen 2011; D'Auria 2006; Effmert 2005.
- **2.4 Statistics.** Competitive gene-set test = Mann-Whitney U of the set's radiation dose
  log2FC vs genome background (signed and |log2FC|), per tier and route; sign test for direction;
  interaction terms for genotype/preservative confounds; `phase3_geneset_test.py`. State
  multiple-testing handling.
- **2.5 Cross-species conservation.** Ensembl Plants Compara REST homology anchored on 12
  Arabidopsis scent genes (`scent_orthology.py`); literature + NCBI-grounded genes for
  Petunia/snapdragon/rose absent from Ensembl (`scent_reference_species.tsv`,
  `scent_query_accessions.tsv`); OrthoFinder protocol for full linkage (`phase4c_protocol.md`).

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
