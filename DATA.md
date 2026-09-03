# Data files & accessions

This repo tracks **code and small metadata** in git. Large data — reference genomes,
alignment/quantification indexes, and bulky derived result tables — is **kept on disk
but git-ignored** (see `.gitignore`) and is intended to travel with the Zenodo data
deposit rather than in the git history. This is what previously inflated `.git` to
~240 MB.

## Raw sequencing data

| Dataset | Accessions | Source | Notes |
|---|---|---|---|
| Floral scent | `SRR4417237`–`SRR4417244` (8 runs) | NCBI SRA → processed via NASA OSDR/GeneLab | High- vs Low-scent lines |
| Radiation / GCR | OSDR accession — **TODO add** | NASA OSDR (collaborator-led) | 0/40 cGy × WT/antho_less × DRS/RL preservative; counts here are our own exploratory analysis |

## Reference genomes (git-ignored, on disk)

| Assembly | Location | Use |
|---|---|---|
| Wisconsin Fast Plant FPsc v1.3 (`BrapaFPsc_277_v1.3`, Phytozome) | `FPsc_genome/`, `Brapa_FPsc_v1_3/` | scent-dataset alignment; `Bra…`-adjacent |
| NCBI RefSeq `GCF_000309985.2` (CAAS Brap v3.01, Chiifu) | `Brapper_fastq/Reference/` | alternate annotation; `BRA…` IDs |

## Git-ignored file categories (provided via deposit)

- Reference sequence/annotation: `*.fa.gz`, `*.fna`, `*.gff3.gz`, `*.gtf.zip`
- STAR / RSEM indexes and intermediates: `Brapper_fastq/STAR_index/`,
  `Brapper_fastq/RSEM_output/`, `*.tab`, `*.results`, `*.theta`, `*.stat`
- Bulky HTML reports (iDEP DEG reports, MultiQC, pathway workflows): `*.html`
- Large derived tables: `*Results_LFC_Pval_DESeq2.csv`,
  `differential_expression*GLbulkRNAseq*.csv`, raw count matrices

## Key small tables kept in git

- `Metadata/Brapp_Scent_metadata_factors.csv` — scent sample sheet (High/Low)
- `counts and factors/LEAF_Metadata.csv` — radiation sample sheet (clean, canonical)
- `brapa_symbol_to_kegg_id_map.csv` — **rebuilt** (was derived from a wrong-organism
  source; now routed via Arabidopsis orthologs, carries `brapa_gene` as a join key)
- `annotation/brapa_to_arabidopsis_orthologs.tsv` — **rebuilt**, 18,770 RBH pairs
  (was a 100-byte BioMart error message)
- `annotation/build_annotation.py` — regenerates all of the above from KEGG + the RBH map
- `Brapa_analysis/Metadata/brapa_kegg_*.tsv` — KEGG gene→pathway / pathway→name
