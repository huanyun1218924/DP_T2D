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

## 🔬 Sample data
1. model_training_sample.Rdata: used for elastic net model training and initial testing

2. signature_sample.RData: contain signature information, test sample for signature calculation, and annotation file

3. t2d_sample: used for prospective analysis between dietary pattern scores and metabolomic signatures and type 2 diabetes

4. micro_sample: used for microbiome analysis


## 🧬 Main analysis code
1. Figure 2.R: Model training to develop metabolomic signatures for dietary patterns; Caculation of metabolomic signature; Pearson 
   correlation between dietary pattern scores and metabolomic signatures; Count of metabolite subclass in each metabolomic signature.
   
2. Figure 3.R: Feature of dietary metabolomic signatures.
  
3. Figure 4.R: Association between dietary metabolomic signatures and incident type 2 diabetes and mediation analysis.
  
4. Figure 5.Rmd: Two-sample MR analysis.
 
5. Figure 6.R: Identify species asscoaited metabolomic signatures and mediation analysis.


## 💻 Environment & Dependencies
R version 4.3.3 (2024-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Rocky Linux 9.7 (Blue Onyx)

Matrix products: default
BLAS:   /app/R-4.3.3@i86-rhel9.0/lib64/R/lib/libRblas.so 
LAPACK: FlexiBLAS OPENBLAS-OPENMP;  LAPACK version 3.9.0

The following core packages are required (full list available in the /functions folder):
- Data processing: tidyr_1.3.2    plyr_1.8.9  dplyr_1.2.0  data.table_1.18.2.1
- Model training:  cvTools_0.3.3  glmnet_4.1-10
- Data analysis:   survival_3.8-3 metafor_4.8-0  coxme_2.2-22  lmerTest_3.2-1  Maaslin2_1.16.0
- Visualization:   ComplexHeatmap_2.25.3  RColorBrewer_1.1-3  circlize_0.4.15  gridExtra_2.3  ComplexHeatmap_2.25.3   circlize_0.4.15  venn_1.12  cowplot_1.2.0    


## 📧 Contact
For questions regarding the analysis or code, please contact:

Huan Yun, Harvard T.H. Chan School of Public Health
huanyun@hsph.harvard.edu
