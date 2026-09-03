# Project Plan — *Brassica rapa*: Floral Scent Under Radiation

**Status:** DRAFT for review by R. Barker · **Drafted:** 2026-07-23
**Repo:** `B_rappa_LLGCSS` · **License:** CC0-1.0

This file is a working plan, not a results document. It (1) states honestly what
is and is not in the repo today, (2) answers — as far as the current data allows —
the question *"does GCR radiation influence floral scent in B. rapa?"*, and
(3) lays out the path to turn this repo into the scent-under-radiation follow-up
manuscript, including the comparative transcriptomics / comparative genomics
(Ensembl orthology) component. It closes with a repo tidy-up and Zenodo checklist.

---

## 1. What is actually in this repo right now

Two **separate** RNA-seq experiments are interleaved here. Keeping them straight is
the single most important thing for the new manuscript.

| | **Scent dataset** | **Radiation / GCR dataset** |
|---|---|---|
| Biology | High- vs Low-scent *B. rapa* lines | Galactic-cosmic-ray dose response |
| Samples | 8 (`SRR4417237`–`SRR4417244`) | 39 `*_GLbulkRNAseq_sub` libraries |
| Factor(s) | `Group` = High / Low | `condition` 0 / 40 cGy · `genotype` WT / antho_less · `preservative` DRS / RL |
| Gene IDs | `Bra000001` (Chiifu / Ensembl Plants) | `BRA028087` (FPsc / RefSeq annotation) |
| Where | `Metadata/`, `NewTest/`, notebooks `01`/`02` | `counts and factors/`, `Radiation model/`, `treatment_genotype + preservative/`, `Exporatory_Model_V1/` |
| Status | Processed → counts → DE table (High vs Low) | Processed → counts → iDEP DESeq2 DEGs (4 contrasts + interaction terms) |

**Provenance / why it's messy.** All the notebooks are adaptations of the NASA
GeneLab **GL4U OSD-104** mouse-spaceflight training module, repurposed first for the
scent data (the practice run) and later extended toward the radiation data. Work in
Gitpod and with Jules/Gemini left behind duplicated notebooks, half-finished edits,
and — importantly — leftover mouse/OSD-104 text still sitting in the markdown and
metadata references.

**Known problems to be aware of before trusting anything here:**

- **`README.md` is aspirational, not real.** It describes a clean `notebooks/`,
  `data/`, `reference/`, `results/`, `scripts/` layout plus `METHODS.md`,
  `DATA-DICTIONARY.md`, `QUICKSTART.md`, `CITATION.cff` — **none of those files or
  folders exist.** It must be rewritten to match reality before any deposit.
- **Gene-ID namespaces don't match.** Scent = `Bra000001`; radiation = `BRA028087`.
  These are different annotations of *B. rapa* and **cannot be joined directly** —
  this is the core obstacle to answering the scent×radiation question (see §2).
- **Radiation DEGs are ~95% unannotated.** 40,540 of 42,278 rows are bare `BRA…`
  locus IDs; the only "named" entries are rRNA / ERCC spike-ins. No scent gene is
  visible by symbol in the radiation results.
- ~~**The ortholog map is partial.**~~ **RESOLVED.** `NewTest/B.rapa_to_Ara_mart_export.txt`
  reached ~3,500 genes, and `annotation/brapa_to_arabidopsis_orthologs.tsv` was a 100-byte
  file holding a BioMart error (Ensembl Plants exposes no Arabidopsis-homolog attribute
  for *B. rapa*, so that query cannot be made to work). Rebuilt by DIAMOND reciprocal best
  hit against the TAIR10 proteome: **18,770 one-to-one pairs**, agreeing with Ensembl's own
  one-way call on 96.6% of genes where both make one. See `annotation/build_annotation.py`.
- ~~**The KEGG mapping is built on a wrong-organism file.**~~ **RESOLVED.** The code `bra`
  is *Bradyrhizobium* sp. ORS 278; *B. rapa* is `brp`. Correcting it is necessary but not
  sufficient: KEGG keys `brp` on NCBI GeneIDs, which join **0** of this repo's `Bra######`
  counts, and only 495 of 44,411 brp genes carry a symbol. The annotation is therefore
  routed through Arabidopsis orthologs — KEGG `ath` is keyed on AGI codes and richly
  symbolled — giving **10,922 pathway links over 4,400 genes** that join to the counts.
  Correct-organism brp reference files are kept in `annotation/kegg/`.
- **Git bloat.** `.git` is ~241 MB: genome FASTAs, STAR indexes, and 12–16 MB DEG
  CSVs were committed. No `.gitignore`.
- **`brapa_gsea.R` is orphaned from its inputs.** It expects `Brapa_metadata.csv`,
  `Brapa_analysis/dds.rds`, and `…/differential_expression.csv`, none of which are in
  the repo — so it will not run as-is.

---

## 2. The scientific question — can we answer it yet?

> **Q: Does GCR radiation potentially influence floral scent in *B. rapa*?**

**Short answer: not from the repo as it stands today — but the data needed to test
it is here, and the analysis is tractable.** Three things currently block a direct
answer:

1. **Different experiments, different genomes.** The scent signal (High vs Low) lives
   in `Bra…` IDs; the radiation signal (40 vs 0 cGy) lives in `BRA…` IDs. Nothing
   currently links a scent gene to its radiation-response value.
2. **No scent gene set is defined.** Neither dataset has been analysed against a
   curated floral-volatile biosynthesis gene set (terpene synthases, benzenoid /
   phenylpropanoid, LOX / green-leaf-volatile, apocarotenoid pathways).
3. **Annotation is too sparse to bridge on.** A quick probe using the existing
   ortholog map surfaced essentially no annotated scent genes — because the map only
   covers ~15% of genes, not because scent biology is absent. This is an annotation
   gap, not a biological result.

**So the honest, defensible statement for a manuscript today is:** *this repository
establishes the two halves needed to address the question — a radiation dose-response
transcriptome and a scent-associated transcriptome in the same species — but the link
between them has not yet been tested.* The plan in §4 is exactly the work that turns
"we have both halves" into a real, testable answer, with a clear null hypothesis
(radiation does not preferentially perturb scent-biosynthesis genes) to reject or fail
to reject.

---

## 3. Assets worth keeping (the reusable / novel pieces)

Everything below is genuinely reusable; the rest is stock GeneLab scaffolding or cruft.

- **`brapa_gsea.R`** (531 lines) — GeneLab-style DGE→PCA→volcano→heatmap→**KEGG GSEA +
  pathview** pipeline adapted for a non-model plant. The custom `TERM2GENE` /
  `TERM2NAME` GSEA and `pathview(species="bra")` calls are the novel part. Needs its
  inputs reconnected (see §1).
- **`brapa_sbgnview.R`** — SBGNview pathway-painting workflow (glucosinolate example;
  keyword-swappable to scent pathways).
- **`create_mapping.R` / `gene_id_mapping.R`** — Ensembl Plants biomaRt queries for
  *B. rapa* symbol↔Ensembl↔KEGG. The right approach for rebuilding the ortholog map.
- **`r_error_handling.R`** — reusable `check_required_packages` / `safe_library` /
  `check_file_exists` guards; worth keeping as a shared utility.
- **`01-RNAseq_processing_Brapa_FPsc_v1_3.ipynb`** — the clean, self-consistent
  fastq→counts pipeline (FastQC→TrimGalore→STAR→RSeQC→RSEM). The canonical processing notebook.
- **`02-RNAseq_analysis_Brapa.ipynb`** (== `_executed`) — the fully ported *B. rapa*
  DGE + KEGG-GSEA + pathview notebook. The canonical analysis notebook.
- **`Radiation model/` interaction-term outputs** (`Radi_model_with_interactions_…csv`)
  — dose×preservative and dose×genotype interaction effects; not reproduced elsewhere.
- **`NewTest/` scent DE table + `B.rapa_to_Ara_mart_export.txt`** — the only scent-side
  DE results and the (partial) orthology bridge; the seed for §4's annotation work.
- **`counts and factors/LEAF_Metadata.csv`** — clean canonical radiation metadata.
- **`environment.yml`** — pinned, reproducible GeneLab RNAseq conda env.

**Redundant / droppable:** `01-…test-checkpoint.ipynb` (stale, mixed mouse/Brapa
refs), `02-…UPDATE_needed.ipynb` (un-ported mouse predecessor), `02-…Brapa.ipynb.bak`,
`02-…Brapa_executed.ipynb` (near-duplicate of `Brapa.ipynb`), numbered scratch PNG/TIFF
duplicates in `Exporatory_Model_V1/`, and the redundant factor-file versions in
`counts and factors/`. `treatment_genotype + preservative/` overlaps heavily with
`Radiation model/` — consolidate to one canonical genotype×dose model.

---

## 4. Plan — from messy repo to scent-under-radiation manuscript

Each phase is independent enough to stop and review after. Phase 0 is the original
tidy-up ask; phases 1–4 are the new research direction; phase 5 is deposit + writing.

### Phase 0 — Consolidate & de-duplicate the repo *(low risk, do first)*
- Adopt a real directory layout and move files into it: `notebooks/`, `scripts/`,
  `data/scent/`, `data/radiation/`, `reference/`, `results/`, `docs/`.
- Keep exactly one canonical copy of each notebook/model (list in §3); archive the
  rest under `docs/legacy/` rather than deleting outright.
- Add a `.gitignore` for genome FASTAs, STAR/RSEM indexes, and large derived CSVs.
- **Rewrite `README.md` to describe what actually exists** (replaces the current
  aspirational one). Delete or fill the templated `SECURITY.md`.
- Strip leftover mouse/OSD-104 markdown from the kept notebooks.
- Decide git-history strategy: the 241 MB history holds large binaries. Options —
  (a) leave as-is, (b) `git gc`/BFG to purge blobs from history (rewrites history),
  or (c) start a clean orphan branch for the deposit. **Needs your call** (see §5).

### Phase 1 — Reconcile gene IDs & build a real annotation layer — **DONE**
- ~~Build a full `Bra…` ↔ Arabidopsis (`AT…`) ortholog table via Ensembl Plants biomaRt~~
  — biomaRt cannot do this: there is no `athaliana_eg_homolog_ensembl_gene` attribute for
  *B. rapa*. Built instead by DIAMOND reciprocal best hit against TAIR10 (18,770 pairs).
- ~~Re-derive KEGG annotations from the correct organism~~ — done; `brp`, not `bra`, and
  routed via Arabidopsis so the ids join to the counts.
- Run `python3 annotation/build_annotation.py` to regenerate everything below from source.

| File | What it is | Joins to counts |
|---|---|---|
| `annotation/brapa_to_arabidopsis_orthologs.tsv` | 18,770 RBH pairs, `Bra… → AT…` | 18,770 genes |
| `annotation/brapa_gene_to_kegg_via_arabidopsis.tsv` | pathway ↔ Bra gene ↔ AGI ↔ symbol | 4,400 genes |
| `Brapa_analysis/Metadata/brapa_kegg_link.tsv` | TERM2GENE for GSEA, Bra-keyed | 4,400 genes |
| `Brapa_analysis/Metadata/brapa_kegg_pathway_names.tsv` | TERM2NAME, 162 pathways | — |
| `brapa_symbol_to_kegg_id_map.csv` | `kegg_id, symbol, brapa_gene, arabidopsis_agi` | 18,770 genes |
| `annotation/kegg/brp_*.tsv` | correct-organism *B. rapa* reference (NCBI GeneID) | 0 — reference only |

- Still open: GO terms are not yet folded into a single `gene_annotation.tsv`.

### Phase 2 — Define the floral-scent gene set
- Curate an Arabidopsis/Brassica floral-volatile biosynthesis gene set: terpenoid
  (TPS family, MEP/MVA), benzenoid/phenylpropanoid (PAL, BSMT, BEAT/BEBT), LOX/
  green-leaf-volatile (LOX, HPL, ADH), apocarotenoid (CCD1). Cite sources.
- Map it into `Bra…` and `BRA…` space via the Phase-1 table.
- Sanity check: confirm these genes move in the **scent** High-vs-Low DE table (they
  should — this validates the gene set in this species before testing radiation).

### Phase 3 — Test the core question: radiation × scent genes
- Re-annotate the **radiation** DEGs (WT 40 vs 0 cGy is the clean radiation contrast;
  also the antho_less dose contrast and the interaction terms).
- Test the null hypothesis: are scent-biosynthesis genes over-represented / directionally
  shifted among radiation-responsive genes? (GSEA / hypergeometric test with the §2
  gene set; report effect sizes and FDR, not just yes/no.)
- Cross-tabulate: genes that are BOTH scent-associated (High vs Low) AND radiation-
  responsive (40 vs 0 cGy) — the headline candidate list for the manuscript.
- Watch the confounds explicitly: `preservative` (DRS/RL) and `genotype` (antho_less
  is a pigment mutant) are baked into the radiation design and must stay in the model.

### Phase 4 — Comparative transcriptomics / genomics (the "scent in space" narrative)
- Extend orthology beyond Arabidopsis to other **scent-producing flowering plants**
  via Ensembl Plants orthology (e.g. *Petunia*, *Antirrhinum*/snapdragon, rose where
  available, other Brassicaceae) — one-to-one and one-to-many ortholog groups for the
  scent gene set.
- Ask: is the radiation-perturbed scent module conserved across scent producers? Which
  scent pathways are lineage-specific vs shared? This is the comparative-genomics spine
  of the "floral scent may be altered in space" story.
- Optional: overlay onto public spaceflight/irradiation plant datasets (GeneLab OSDR)
  for an external replication angle.

### Phase 5 — Deposit & manuscript framing
- Finalise Zenodo deposit (see §6) with a DOI, and align README/CITATION with it.
- Draft the manuscript around: (i) two-experiment design, (ii) scent gene set,
  (iii) radiation×scent test result, (iv) cross-species conservation.

---

## 5. Decisions I need from you before running Phase 0+

1. **Git history / size** — leave the 241 MB history, purge large blobs (rewrites
   history — anyone who cloned must re-clone), or start a clean deposit branch?
2. **Delete vs archive** — OK to move redundant notebooks/models to `docs/legacy/`
   rather than delete? Or delete outright?
3. **Scent dataset provenance** — do you have the citation/accession for the
   `SRR4417237–244` High/Low-scent study? Needed for DATA docs and to confirm the
   High/Low design.
4. **Genome of record** — for the reconciled analysis, standardise on one *B. rapa*
   assembly (Chiifu `Bra…` vs FPsc/RefSeq `BRA…`)? This drives the Phase-1 mapping.
5. **Radiation data sharing** — the radiation experiment is your collaborators'; are
   its counts/DEGs OK to include in a public Zenodo deposit, or should the deposit be
   scent-focused with radiation referenced by their DOI?
6. **Scope of this repo** — keep it as one combined repo, or split "scent analysis"
   (deposit-ready) from "radiation model" (collaborator-owned)?

---

## 6. Zenodo readiness checklist

- [ ] `README.md` rewritten to match real contents (currently fabricated)
- [ ] `CITATION.cff` created (does not exist yet) with your ORCID `0000-0001-5681-9857`
- [ ] `.zenodo.json` with authors, keywords, funding (Lunar LEAF), related identifiers
- [ ] `LICENSE` present ✅ (CC0-1.0 already in repo)
- [ ] Large binaries handled (`.gitignore` + decide what data ships in the deposit)
- [ ] `environment.yml` verified to build ✅ (present; test once)
- [ ] Data provenance documented (SRA accessions for scent; DOI for radiation)
- [ ] Leftover mouse/OSD-104 text removed from notebooks
- [ ] One canonical notebook + script per step; legacy archived
- [ ] `DATA-DICTIONARY.md` / `METHODS.md` written (referenced by README but missing)
