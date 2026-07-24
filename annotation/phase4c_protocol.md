# Phase 4c — reproducible ortholog protocol (run on your compute)

Goal: replace the literature curation of Phase 4b with a **programmatic** ortholog
assignment linking the model scent species (Petunia, snapdragon, rose) to the B. rapa
scent set. This sandbox has no alignment toolchain, so the linkage itself is documented
here to run where OrthoFinder/DIAMOND are available. The grounded inputs
(`scent_query_proteins.faa`, `scent_query_accessions.tsv`) are already produced by
`phase4c_link.py`.

## 1. Proteomes to collect (grounded sources)

| Species | Source | Accession / location |
|---|---|---|
| *Arabidopsis thaliana* | Araport11 (TAIR) or Ensembl Plants | `Araport11_pep` / `arabidopsis_thaliana` pep FASTA |
| *Brassica rapa* (Chiifu) | Ensembl Plants (matches our `Bra######` IDs) | FTP `…/plants/fasta/brassica_rapa/pep/` or BioMart `brapa_eg_gene` attribute `peptide` |
| *Petunia axillaris* | **Sol Genomics** (Bombarely 2016; NCBI copies are unannotated) | `https://solgenomics.net/ftp/genomes/Petunia_axillaris/` (reachable) |
| *Antirrhinum majus* | **Li et al. 2019, Nat. Plants** genome portal (absent from NCBI datasets & Ensembl) | A. majus cv. JI7 proteome from the dedicated genome DB |
| *Rosa chinensis* 'Old Blush' | **NCBI RefSeq** (Raymond 2018; annotated) | `GCF_002994745.2` (RchiOBHm-V2) |

```bash
# example: Rosa proteome via NCBI datasets CLI (robust, no fragile FTP path)
datasets download genome accession GCF_002994745.2 --include protein
# B. rapa Chiifu peptides via Ensembl Plants BioMart (attribute 'peptide' on brapa_eg_gene)
```

Put every proteome (one FASTA per species) + `scent_query_proteins.faa` into `proteomes/`.

## 2. Run OrthoFinder (or a DIAMOND reciprocal-best-hit)

```bash
conda create -n ortho -c bioconda orthofinder diamond mcl -y && conda activate ortho
orthofinder -f proteomes/ -t 8            # -> Orthogroups.tsv, Orthologues/
# lightweight alternative — DIAMOND RBH of the 7 query proteins vs B. rapa Chiifu:
diamond makedb --in brassica_rapa.pep.fa -d brapa
diamond blastp -q scent_query_proteins.faa -d brapa -k 5 --very-sensitive -o q2brapa.tsv
diamond makedb --in scent_query_proteins.faa -d q
diamond blastp -q brassica_rapa.pep.fa -d q -k 5 --very-sensitive -o brapa2q.tsv
# reciprocal best hits = candidate orthologs
```

## 3. Intersect with the scent set

```bash
# keep orthogroups that contain BOTH a query accession (scent_query_accessions.tsv)
# AND a B. rapa gene in annotation/scent_geneset.tsv -> the model-species↔B. rapa scent links.
```

Expected: the ester/benzenoid methyltransferase queries (PhBSMT `A4ZDG8`, AmBAMT
`Q9FYZ9`, RhOMT `AAM23005`) should recover a B. rapa SABATH/O-MT orthogroup that includes
`Bra029041` — turning the Phase-4/4b "same tailoring family" argument into a hard
orthogroup membership call. RhNUDX1 (`M4I1C6`, a Nudix hydrolase) is expected to have **no**
terpene-synthase ortholog — the built-in positive control for pathway convergence.
