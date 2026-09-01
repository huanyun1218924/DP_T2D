**Dietary Patterns, Circulating Metabolome, and Risk of Type 2 Diabetes**


Overview:
This repository contains the simulated sample data and main analysis code for a integrative study investigating the associations between diet, circulating metabolome in relation to incident type 2 diabetes. The analysis utilizes diet, metabolomics, genetics, metagenomics, and clinical phenotypes.


Sample data:
1. model_training_sample.Rdata: used for elastic net model training and initial testing

2. signature_sample.RData: contain signature information, test sample for signature calculation, and annotation file

3. t2d_sample: used for prospective analysis between dietary pattern scores and metabolomic signatures and type 2 diabetes

4. micro_sample: used for microbiome analysis


Code files in this repository covers the main analyses:
1. Figure 2.R: Model training to develop metabolomic signatures for dietary patterns; Caculation of metabolomic signature; Pearson 
   correlation between dietary pattern scores and metabolomic signatures; Count of metabolite subclass in each metabolomic signature.
   
2. Figure 3.R: Feature of dietary metabolomic signatures.
  
3. Figure 4.R: Association between dietary metabolomic signatures and incident type 2 diabetes and mediation analysis.
  
4. Figure 5.Rmd: Two-sample MR analysis.
 
5. Figure 6.R: Identify species asscoaited metabolomic signatures and mediation analysis.


Environment & Dependencies
R Packages
The following core packages are required (full list available in the /functions folder):

Statistics: lmerTest (v3.1.3), ppcor (v1.1), vegan (v2.6.8), bnlearn (v5.1), mediation (v4.5.0)
Machine Learning: randomForest (v4.7.1.2), iml
Bioinformatics: clusterProfiler (v4.14.4), org.Hs.eg.db (v3.20.0)
Visualization: ggplot2 (v3.5.1), ggpubr (v0.6.0), igraph (v2.1.2)
Data Handling: Hmisc (v5.2.2), compositions (v2.0.8)


Contact
[Huan Yun/Department of Epidemiology, Harvard T.H. Chan School of Public Health]
[huanyun@hsph.harvard.edu]
