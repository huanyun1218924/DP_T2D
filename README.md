# Dietary Patterns, Circulating Metabolome, and Risk of Type 2 Diabetes


## 📌 Overview

This repository contains the simulated sample data and main analysis code for a integrative study investigating the associations between diet, circulating metabolome in relation to incident type 2 diabetes. The analysis utilizes diet, metabolomics, genetics, metagenomics, and clinical phenotypes.


## 📂 Repository Structure
```text
├── data/
│   ├── model_training_sample.RData
│   ├── signature_sample.RData
│   ├── t2d_sample.RData
│   ├── micro_sample.RData
│
├── code/
│   ├── Figure 2.R
│   ├── Figure 3.R
│   ├── Figure 4.R
│   ├── Figure 5.Rmd
│   └── Figure 6.R
│
└── README.md

Note: Individual-level data are not included in this repository because of data-use and participant confidentiality restrictions.
```

## 🔬 Simulated Sample Data
- `model_training_sample.Rdata`: used for elastic net model training and initial testing for developing dietary metabolomic signatures, with sample size and variable structure same with raw dataset
  
- `signature_sample.RData`: contain signature information, test sample for signature calculation, and metabolite annotation file

- `t2d_sample`: used for prospective analysis between dietary pattern scores and metabolomic signatures and type 2 diabetes

- `micro_sample`: used for microbiome analysis, with sample size and variable structure same with raw dataset


## 🧬 Main Analysis Code
- `Figure 2.R`: Develop metabolomic signatures for dietary patterns using elastic net regression; Caculation of metabolomic signature in each study; Pearson 
   correlation between dietary pattern scores and corresponding metabolomic signatures; Count of metabolite subclass in each metabolomic signature.

- `Figure 3.R`: Feature of dietary metabolomic signatures, including its coefficient from elastic net model, its association with corresponding metabolomic signatures using linear regression, and the association between composite metabolites and incident type 2 diabetes.

- `Figure 4.R`: Association between dietary metabolomic signatures and incident type 2 diabetes using either cox or logistic regression model; Mediation analysis among dietary pattern scores, composite metabolites, and incident type 2 diabetes.

- `Figure 5.Rmd`: Genetic heritability estimation using LDSC; TWAS analysis using FUSION; tissue-specific colocalization analysis; genetic correlation analysis; and two-sample Mendelian randomization analysis. Genetic analyses were primarily conducted using a Linux server.
 
- `Figure 6.R`: Identify species associated metabolomic signatures and composite metabolites; Mediation analysis among dietary pattern scores, species, and corresponding metabolomic signatures.


## 💻 Environment & Dependencies
- R version: `v4.3.3 (2024-02-29)`
- Platform: `x86_64-pc-linux-gnu (64-bit)`
- Running under: `Rocky Linux 9.7 (Blue Onyx)`
- BLAS:   `/app/R-4.3.3@i86-rhel9.0/lib64/R/lib/libRblas.so`
- LAPACK: `v3.9`

The following core packages are required (full list available in the /functions folder):
- **Data processing**: `tidyr_1.3.2`,    `plyr_1.8.9`,  `dplyr_1.2.0`,  `data.table_1.18.2.1`
- **Model training**:  `cvTools_0.3.3`,  `glmnet_4.1-10`
- **Data analysis**:   `survival_3.8-3`, `metafor_4.8-0`,  `coxme_2.2-22`,  `lmerTest_3.2-1`,  `Maaslin2_1.16.0`
- **Visualization**:   `ggplot2_4.0.2`, `ComplexHeatmap_2.25.3`,  `RColorBrewer_1.1-3`,  `circlize_0.4.15`,  `gridExtra_2.3`,  `venn_1.12`,  `cowplot_1.2.0`   


## 📧 Contact
For questions regarding the analysis or code, please contact:

Huan Yun, Harvard T.H. Chan School of Public Health, huanyun@hsph.harvard.edu
