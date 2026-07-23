# *Brassica rapa* — Radiation-Induced Change in Floral Scent (LLGCSS)

RNA-seq analysis code and data for a *Brassica rapa* (Wisconsin Fast Plant) project
investigating whether **galactic cosmic ray (GCR) radiation alters floral scent**.
Part of the **Lunar LEAF** (Life sciences & Engineering for lunar Agriculture and
Flight) effort.

> **Status:** actively being consolidated and developed into a follow-up manuscript.
> See **[`README_PLAN.md`](README_PLAN.md)** for the full research + tidy-up plan,
> the current state of the analysis, and open decisions. Read that first if you are
> picking this work up.

---

## What this repository actually contains

This repo brings together **two separate RNA-seq experiments** in the same species.
Keeping them distinct matters — they use different samples *and* different genome
annotations.

### 1. Scent dataset (the practice / seed dataset)
- 8 samples, **High- vs Low-scent** *B. rapa* lines.
- SRA accessions `SRR4417237`–`SRR4417244`; reads were processed through NASA
  **OSDR/GeneLab**-style tooling to obtain counts.
- Gene IDs: Chiifu / Ensembl Plants (`Bra000001` style).
- Files: `Metadata/Brapp_Scent_metadata_factors.csv`, `NewTest/` (DE table +
  B. rapa→Arabidopsis ortholog export), notebooks `01`/`02`.

### 2. Radiation / GCR dataset
- 39 GeneLab-style libraries. Factors: `condition` = 0 / 40 cGy, `genotype` = WT /
  *anthocyaninless*, `preservative` = DRS / RL.
- Analysed with iDEP-generated DESeq2 models (main effects + interaction terms).
- Gene IDs: FPsc / RefSeq (`BRA028087` style).
- Files: `counts and factors/`, `Radiation model/`, `treatment_genotype + preservative/`,
  `Exporatory_Model_V1/`.
- **The counts/results here are our own exploratory analysis.** The primary radiation
  study is being carried forward separately by collaborators; treat OSDR as the source
  of record (accession/DOI to be linked).

> **Important:** the two datasets' gene-ID namespaces (`Bra…` vs `BRA…`) do **not**
> join directly. Reconciling them is the first analytical step toward the scent×radiation
> question — see `README_PLAN.md` §4.

---

## Repository layout

```
01-RNAseq_processing_Brapa_FPsc_v1_3.ipynb   fastq → counts (FastQC→TrimGalore→STAR→RSeQC→RSEM)
02-RNAseq_analysis_Brapa.ipynb               counts → DESeq2 DGE → PCA/volcano/heatmap → KEGG GSEA + pathview
brapa_gsea.R                                 standalone DGE + KEGG GSEA + pathview pipeline (non-model plant)
brapa_sbgnview.R                             SBGNview pathway painting
create_mapping.R / gene_id_mapping.R         Ensembl Plants biomaRt: B. rapa symbol ↔ Ensembl ↔ KEGG
r_error_handling.R                           reusable package/file/dir guards, sourced by the R scripts
download_file_urls.py                        fetch files from an OSDR/GLDS accession by pattern
process_genes.py                             build a symbol↔KEGG map (see caveat below)
environment.yml                              pinned conda env (GeneLab RNAseq)
Metadata/, counts and factors/               experiment metadata + count/factor matrices
Brapper_fastq/                               scent-dataset processing tree (QC, STAR/RSEM outputs)
Radiation model/, treatment_genotype + preservative/, Exporatory_Model_V1/   radiation DESeq2 models
NewTest/                                     scent DE results + B. rapa→Arabidopsis ortholog export
Brapa_analysis/                              KEGG annotation files (gene→pathway, pathway→name)
FPsc_genome/, Brapa_FPsc_v1_3/               Wisconsin Fast Plant (FPsc v1.3) reference
images*/                                     figures from exploratory analyses
docs/                                        LEGACY.md manifest + on-disk legacy archive
README_PLAN.md                               research + consolidation plan (read this)
```

### Reference genomes (kept on disk, not in git)
Two *B. rapa* assemblies are present for the comparative-genomics work:
- **Wisconsin Fast Plant — FPsc v1.3** (Phytozome `BrapaFPsc_277_v1.3`), in `FPsc_genome/`.
- **NCBI RefSeq `GCF_000309985.2` (CAAS Brap v3.01, Chiifu-based)**, in `Brapper_fastq/Reference/`.

Large reference, index, and derived-result files are **git-ignored** and provided via
the data deposit rather than committed (they are what previously bloated the history).

---

## Reusable / notable pieces

- **`brapa_gsea.R`** — a KEGG GSEA + `pathview` pipeline for a *non-model* plant
  (custom `TERM2GENE`/`TERM2NAME`), the most reusable code here. *Currently orphaned
  from its inputs (`Brapa_metadata.csv`, `Brapa_analysis/dds.rds`) — see plan.*
- **`create_mapping.R` / `gene_id_mapping.R`** — the right approach for rebuilding a
  full ortholog/ID map via Ensembl Plants biomaRt.
- **`01`/`02` notebooks** — the canonical processing and analysis pipeline (adapted from
  the NASA GeneLab GL4U training module).

## Known caveats (before you trust an output)

- `process_genes.py` and the shipped `brapa_symbol_to_kegg_id_map.csv` were derived from
  a file that turned out to contain the **wrong organism** (*Bradyrhizobium*, not
  *Brassica*); that source file is archived in `docs/legacy/`. Rebuild KEGG/ortholog
  maps before relying on pathway results.
- The notebooks still carry some leftover mouse/OSD-104 narrative text from the GeneLab
  template; a cleanup pass is planned.
- The in-repo ortholog map (`NewTest/B.rapa_to_Ara_mart_export.txt`) covers only ~15% of
  genes and needs a full rebuild.

---

## Quick start

```bash
conda env create -f environment.yml
conda activate gl4u_rnaseq_2024
jupyter lab
# open 01-RNAseq_processing_Brapa_FPsc_v1_3.ipynb, then 02-RNAseq_analysis_Brapa.ipynb
```

## Data availability
- **Scent reads:** NCBI SRA `SRR4417237`–`SRR4417244`; processed via NASA OSDR/GeneLab.
- **Radiation data:** collaborator-led; cite OSDR accession/DOI (link to be added).

## Citation
See [`CITATION.cff`](CITATION.cff). A Zenodo DOI will be minted on deposit.

## License
Released under **CC0 1.0 Universal** (public domain). See [`LICENSE`](LICENSE).

## Author
**Dr. Richard Barker** — ORCID [0000-0001-5681-9857](https://orcid.org/0000-0001-5681-9857)
