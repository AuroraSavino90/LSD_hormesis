# LSD, Hormesis, and Stress: Cross-Dataset Transcriptomic Analysis

This repository contains the analysis code used to compare transcriptional responses induced by LSD with those observed under hormetic stress paradigms (e.g. fasting, exercise, preconditioning) and chronic stress across multiple publicly available transcriptomic datasets.

The central goal of the project is to test whether LSD-induced gene-expression programs are more closely aligned with adaptive (hormetic) stress responses rather than maladaptive chronic stress signatures, using direction-aware, pathway-level, and multivariate approaches.

---

## Analysis workflow

The scripts are intended to be run in roughly the following order:

**`01_LSD_PFC.R`**
Quality control, differential expression, and pathway enrichment analysis for LSD-treated samples.

**`02_External_datasets.R`**
Loads and standardizes external transcriptomic datasets representing hormetic interventions (fasting, exercise, ischemic preconditioning) and chronic stress paradigms.

**`03_Overlap.R`**
Basic overlap analysis between LSD differentially expressed genes (DEGs) and external datasets.

**`04_Fisher_ormesivsstress.R`**
Fisher exact tests comparing enrichment of LSD signatures in hormesis versus chronic stress datasets.

**`05_Overlap_alldir.R`**
Direction-resolved overlap analyses distinguishing same-direction regulation (up–up, down–down) and opposite-direction regulation (up–down).

**`06_Tests_for_sign.R`**
Focuses on directional enrichment of LSD transcriptional signatures within hormetic datasets.

**`07_Tests_for_sign_stress.R`**
Applies the same direction-resolved analytical framework used for hormesis to chronic stress datasets.

**`08_alldatasets_PCA.R`**
Unsupervised PCA on pathway enrichment matrices across all datasets.

**`09_NES_sPLSDA.R`**
Sparse PLS-DA trained on pathway enrichment values to discriminate hormesis versus chronic stress datasets, with LSD datasets projected into the learned space.

**`10_LODO.R`**
Leave-one-dataset-out (LODO) analysis to assess robustness and dataset dependence of the sPLS-DA results.

---

## Reproducibility

Exact package versions used in this analysis are recorded in:

- `package_versions.txt`


