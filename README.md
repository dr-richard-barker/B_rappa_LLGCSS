# B. rappa RNA-seq Analysis: Lunar LEAF LLGCSS Project

## Overview

This repository contains the complete RNA-sequencing analysis pipeline for *Brassica rappa* (Wisconsin Fast Plant) exposed to simulated galactic cosmic ray radiation. This work supports the Lunar LEAF (Life Sciences & Engineering for Lunar Environments and Applications for Flight) program and provides baseline molecular responses for plant adaptation studies relevant to future lunar agriculture.

## Project Description

**Purpose:** Investigate transcriptomic responses of *B. rappa* to cosmic ray simulation to establish baseline data for the Lunar LEAF project's genetic expression studies.

**Organism:** *Brassica rappa* (Wisconsin Fast Plant)  
**Reference Genome Version:** BrapaFPsc v1.3 (GCF_000309985.2)  
**Analysis Method:** GeneLab standard RNAseq pipeline  
**License:** CC0 1.0 Universal (Public Domain)

## Contents

### Notebooks (Analysis Pipeline)
- **`notebooks/01-preprocessing-RNAseq.ipynb`** - Raw data quality control, trimming, alignment, and quantification (fastq → counts)
- **`notebooks/02-differential-expression.ipynb`** - Normalization, differential expression analysis, visualization, and pathway enrichment

### Data Organization
- **`data/metadata/`** - Experimental design and sample information
- **`data/processed/`** - Intermediate and final analysis outputs
- **`reference/`** - Reference genome and annotation files
- **`results/`** - Final figures and summary tables

### Code & Scripts
- **`scripts/`** - Reusable R and Python functions for analysis
- **`environment.yml`** - Complete conda environment specification for reproducibility
- **`.gitpod.yml`** - Configuration for cloud-based development

### Documentation
- **`METHODS.md`** - Detailed methods and pipeline specifications
- **`DATA-DICTIONARY.md`** - Description of all data files and formats
- **`QUICKSTART.md`** - Quick start guide for running the analysis
- **`CITATION.cff`** - Citation metadata

## Quick Start

### Option 1: Local Installation

```bash
# Clone repository
git clone https://github.com/dr-richard-barker/B_rappa_LLGCSS.git
cd B_rappa_LLGCSS

# Create conda environment
conda env create -f environment.yml
conda activate gl4u_rnaseq_2024

# Start Jupyter Lab
jupyter lab
```

Then open `notebooks/01-preprocessing-RNAseq.ipynb` to begin.


## System Requirements

- **Python:** 3.10
- **R:** 4.4.1
- **Memory:** 16 GB recommended
- **Disk Space:** ~100 GB for raw sequencing data (not included in repository)
- **Operating System:** Linux/macOS preferred; Windows with WSL2 supported

## Installation

### Prerequisites
- Conda or Manconda installed ([Installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html))
- Git

### Setup Steps

```bash
# Create environment with all dependencies
conda env create -f environment.yml

# Activate environment
conda activate gl4u_rnaseq_2024

# (Optional) Install additional development tools
pip install nbclassic
```

**Dependencies Include:**
- **Bioinformatics tools:** STAR, RSEM, Trim Galore!, FastQC, MultiQC, samtools
- **R packages:** DESeq2, tidyverse, ggplot2, clusterProfiler, pathview
- **Python packages:** pandas, numpy, scipy, jupyter

See `environment.yml` for complete version specifications.

## Usage

### Step 1: Prepare Raw Data
```bash
# Download raw FASTQ files from NCBI SRA (see DATA.md for links)
# Place in: data/raw/fastq/
```

### Step 2: Run Preprocessing Pipeline
Open and execute `notebooks/01-preprocessing-RNAseq.ipynb`:
- Quality control with FastQC/MultiQC
- Read trimming with Trim Galore!
- Genome alignment with STAR
- Read quantification with RSEM

**Outputs:** `data/processed/counts/`

### Step 3: Run Differential Expression Analysis
Open and execute `notebooks/02-differential-expression.ipynb`:
- Count normalization with DESeq2
- Statistical testing for differential expression
- Visualization (PCA, heatmaps, volcano plots)
- Gene set enrichment analysis

**Outputs:** `results/`

## Data Availability

### Raw Sequencing Data
- **Source:** NCBI Sequence Read Archive (SRA)
- **Bioproject ID:** [SRA_PROJECT_ID]
- **Download with:**
  ```bash
  fastq-dump --gzip [SRA_RUN_ID]
  ```
- See `DATA.md` for complete list of sample accessions

### Reference Files
- **Genome:** *B. rappa* v1.3 (GCF_000309985.2) from NCBI Datasets
- **Included:** `reference/` directory contains reference sequences
- **Alternative sources:** 
  - NCBI: https://www.ncbi.nlm.nih.gov/datasets
  - Phytozome: https://phytozome-next.jgi.doe.gov/

## Reproducibility

This repository is designed for full reproducibility:

1. **Exact versions:** All software versions pinned in `environment.yml`
2. **Containerization:** Gitpod configuration for cloud reproducibility
3. **Documentation:** Step-by-step methods in `METHODS.md`
4. **Code availability:** All analysis scripts included and annotated

To reproduce all results:
```bash
conda env create -f environment.yml
conda activate gl4u_rnaseq_2024
jupyter execute notebooks/01-preprocessing-RNAseq.ipynb
jupyter execute notebooks/02-differential-expression.ipynb
```

## Methods Summary

**RNAseq Pipeline Overview:**
1. **Quality Assessment:** FastQC evaluates raw read quality metrics
2. **Preprocessing:** Trim Galore! removes adapters and low-quality bases
3. **Alignment:** STAR maps reads to *B. rappa* v1.3 reference genome
4. **Quantification:** RSEM generates transcript abundance estimates
5. **Normalization:** DESeq2 applies variance-stabilizing transformation
6. **Statistical Testing:** DESeq2 performs negative binomial GLM differential expression testing
7. **Visualization & Enrichment:** ggplot2 and clusterProfiler for downstream analysis

For detailed methods, see `METHODS.md`.

## Key Results

[Summary of main findings - to be added after analysis completion]

- Number of differentially expressed genes: [N]
- Treatment effect: [description]
- Key pathways identified: [pathways]
- [Other key results]

See `results/` directory for complete outputs and figures.

## Authors & Contact

**Project Lead:** Dr. Richard Barker  
**Affiliation:** [Your Institution]  
**Contact:** [Your Email]  
**ORCID:** [Your ORCID if available]

## How to Cite

If you use this code or data, please cite:

```bibtex
@software{barker2024brappa,
  title={B. rappa RNA-seq Analysis: Lunar LEAF LLGCSS Project},
  author={Barker, Richard},
  year={2024},
  url={https://github.com/dr-richard-barker/B_rappa_LLGCSS},
  doi={[Zenodo DOI]},
  note={Available on Zenodo}
}
```

Or use the CITATION.cff file: `cite this repository` button on GitHub.

## License

This work is released under the **Creative Commons Zero v1.0 Universal (CC0 1.0)** license, effectively placing it in the public domain. You are free to:
- Use for any purpose
- Modify and adapt
- Distribute commercially and non-commercially
- No attribution required (though appreciated!)

See `LICENSE` file for full legal text.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b improve/feature`)
3. Make your changes
4. Submit a pull request with description

For bug reports or feature requests, please open an issue.

## Funding & Acknowledgments

This research was supported by the **Lunar LEAF (Life Sciences & Engineering for Lunar Environments and Applications for Flight)** program.

We acknowledge:
- [NASA/funding agency]
- GeneLab consortium for RNAseq pipeline standards
- [Any collaborators or data providers]

## References

- GeneLab RNAseq Processing Pipeline: [https://github.com/nasa/GeneLab_Data_Processing](https://github.com/nasa/GeneLab_Data_Processing)
- FastQC: [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- STAR Aligner: [https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)
- RSEM: [https://github.com/deweylab/RSEM](https://github.com/deweylab/RSEM)
- DESeq2: [https://bioconductor.org/packages/DESeq2/](https://bioconductor.org/packages/DESeq2/)

## Troubleshooting

### Issue: Environment creation fails
**Solution:** Try updating conda: `conda update conda`

### Issue: STAR index building runs out of memory
**Solution:** Reduce `--limitGenomeGenerateRAM` parameter in notebooks

### Issue: Gitpod times out
**Solution:** Use local installation for long-running analyses; Gitpod sessions have 30-minute inactivity limit

For more issues, see GitHub Issues or contact the authors.

---

**Last Updated:** July 2024  
**Repository Status:** Active Development  
**Maintained by:** Dr. Richard Barker
