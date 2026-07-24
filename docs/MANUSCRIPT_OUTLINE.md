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
- **Results:** [scent set of 363 genes across 4 routes, validated]; core volatile enzymes are
  unaffected by 40 cGy, but the ester/methyltransferase "tailoring" route shows a weak,
  directional down-shift (MWU p≈0.01, uncorrected) — the same family that is lineage-specifically
  expanded and scent-defining across flowering plants.
- **Conclusions:** We provide a reusable framework and a specific, falsifiable prediction: if
  radiation alters scent, it acts on the tailoring methyltransferases, not the core synthases.

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
- **3.3 Radiation spares the core but weakly suppresses the tailoring route.** Tier-1 core: no
  shift (p=0.40); broad set: down-shift vs background (signed MWU p=0.032); concentrated in the
  ester/methyltransferase route (p=0.012); magnitude unchanged; confounds clean (0 significant
  interactions). Dual-hit candidates `Bra013161`, `Bra028224`, `Bra039555`. → *Fig 3*, *Table 2*.
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

## Figures & Tables (all derivable from `annotation/`)
- **Fig 1** two-experiment design + `Bra`↔`BRA` join schematic.
- **Fig 2** scent gene set: routes × tiers; scent-axis validation (DE High/Low per route).
- **Fig 3** radiation gene-set test: per-route MWU (forest/violin of dose log2FC, scent vs background).
- **Fig 4** cross-species conservation heatmap (`scent_orthology_matrix.tsv`) + expansion callouts.
- **Table 1** the 363-gene scent set summary (from `scent_geneset.tsv`).
- **Table 2** radiation×scent candidate genes (dual hits).
- **Table 3** characterised scent genes across model species (`scent_reference_species.tsv`).

## Key references (verified this project)
Dudareva 2013 New Phytol 198:16; Verdonk 2005 Plant Cell (PhODO1, PMID 15805488); Nagegowda 2008
Plant J (AmNES/LIS, PMID 18363779); Magnard 2015 Science (RhNUDX1, PMID 26138978); Chen 2011 Plant J
(TPS). [Add: GeneLab/OSDR pipeline; DESeq2; iDEP; Ensembl Plants; Brassica triplication.]
