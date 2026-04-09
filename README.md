# Integrative multi-omics analysis identifies putative causal genes underlying ADHD risk

## Overview
This repository provides a reproducible multi-omics analysis pipeline integrating GWAS, brain eQTL, Bayesian colocalization, and SMR analyses to prioritize putative causal genes associated with attention-deficit/hyperactivity disorder (ADHD).

This framework systematically integrates genetic association signals with transcriptomic regulatory evidence to refine GWAS loci into biologically interpretable candidate genes.

---

## Project structure

ADHD_multiomics_prioritization/

- README.md
- data/                    # Input data (not included)
  - Step1/                # GWAS summary statistics
  - Step2/                # LD reference panel
  - Step3/                # GTEx brain eQTL
  - Step5/                # BrainMeta eQTL and auxiliary files

- scripts/                # Analysis scripts
  - 01_GWAS_QC.R
  - 02_LD_clumping.R
  - 03_eQTL_annotation.R
  - 04_coloc_analysis.R
  - Step5/
    - 05_01_SMR_prepare_GWAS.R
    - 05_02_run_SMR_BrainMeta.bat
    - 05_03_merge_SMR_results.R
    - 05_04_SMR_summary_and_filter.R
    - 05_05_SMR_plot_pipeline.R
    - 05_06_run_SMR_8genes.bat
    - 05_07_regional_association_plot_publication.R
  - 06_integration_and_visualization.R

- results/                # Output results
  - Step1 to Step6/

- tools/                  # External tools (optional, not required in repo)
  - plink/
  - smr/

---

## Data sources

The following publicly available datasets were used:

- ADHD GWAS summary statistics  
  Psychiatric Genomics Consortium (PGC)  
  Download: https://pgc.unc.edu/for-researchers/download-results/  
  Demontis et al., Nature Genetics, 2023

- Brain eQTL data  
  GTEx v10 (brain tissues)  
  Download: https://www.gtexportal.org/home/downloads/adult-gtex  
  GTEx Consortium, Science, 2020

- BrainMeta cis-eQTL summary data  
  SMR resource  
  Download: https://yanglab.westlake.edu.cn/software/smr/#DataResource  
  Qi et al., Nature Genetics, 2022

- LD reference panel  
  1000 Genomes Project (Phase 3, European ancestry)  
  Data portal: https://www.internationalgenome.org/data  
  1000 Genomes Project Consortium, Nature, 2015

Due to data size and licensing restrictions, these datasets are not included in this repository.

---

## Software requirements

- R (¡Ý 4.0)

Required R packages:
- data.table
- coloc
- ggplot2
- patchwork
- ggrepel

External tools:
- PLINK (v1.9 or later)
- SMR software (v1.3.1)

---

## Workflow

### Step 1. GWAS quality control
Script:  
scripts/01_GWAS_QC.R

- Standardizes GWAS summary statistics
- Converts OR to log effect size (¦Â = ln(OR))
- Performs quality control filtering

---

### Step 2. LD clumping and SNP expansion
Script:  
scripts/02_LD_clumping.R

- Identifies independent lead SNPs
- Expands proxy SNPs using LD reference

---

### Step 3. Brain eQTL annotation
Script:  
scripts/03_eQTL_annotation.R

- Maps GWAS SNPs to brain eQTL signals
- Extracts top eQTL signals per gene

---

### Step 4. Bayesian colocalization analysis
Script:  
scripts/04_coloc_analysis.R

- Performs coloc analysis
- Estimates posterior probability (PP.H4)
- Classifies colocalization evidence

---

### Step 5. SMR analysis
Directory:  
scripts/Step5/

Includes:

- 05_01_SMR_prepare_GWAS.R
- 05_02_run_SMR_BrainMeta.bat
- 05_03_merge_SMR_results.R
- 05_04_SMR_summary_and_filter.R
- 05_05_SMR_plot_pipeline.R
- 05_06_run_SMR_8genes.bat
- 05_07_regional_association_plot_publication.R

Main steps:
- Preparation of GWAS for SMR
- Chromosome-wise SMR analysis
- Result merging
- HEIDI filtering
- Regional association visualization

---

### Step 6. Multi-evidence integration and visualization
Script:  
scripts/06_integration_and_visualization.R

- Integrates:
  - GTEx brain eQTL evidence
  - Colocalization results
  - SMR + HEIDI results

- Assigns gene-level evidence scores
- Prioritizes candidate genes
- Generates summary tables and figures

---

## Key outputs

### Tables (results/Step6/tables/)
- Step6_final_candidate_genes.tsv
- Step6_final_candidate_snps.tsv
- Step6_highconfidence_genes_summary.tsv
- Step6_multi_evidence_snps.tsv

### Figures (results/Step6/figures/)
- Step6_combined_integrated_evidence.pdf (main figure)
- Step6_gene_evidence_levels.pdf
- Step6_gene_evidence_matrix.pdf
- Step6_snp_evidence_distribution.pdf

### SMR plots (results/Step5/)
- Regional association plots for representative genes (e.g., LSM6)

---

## Reproducibility

All scripts use relative paths based on the repository structure.

For scripts in `scripts/`:
`project_root <- ".."`

For scripts in `scripts/Step5/`, relative paths are defined consistently within the corresponding R scripts or batch files.

To run the pipeline:

1. Place required input data into the `data/` directory
2. Install required tools (PLINK and SMR)
3. Run scripts sequentially from Step1 to Step6

---

## Notes

- Large datasets (GWAS, GTEx, BrainMeta, LD reference) are not included
- Relative paths assume that the repository directory structure is preserved
- Some SMR steps require Windows batch execution

---

## Citation

If you use this pipeline, please cite:

Zhao Q, Li F, et al.  
Integrative multi-omics analysis identifies putative causal genes underlying ADHD risk.

---

## Contact

Zhao Qianlong  
The Second Hospital & Clinical Medical School, Lanzhou University  

Email: zhaoqianl2025@gmail.com