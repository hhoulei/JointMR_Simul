****R Scripts for JointMR Analysis and Data Formats****  

This repository contains R scripts for JointMR simulations and real-data applications, along with the required input data format specifications.    

Article: JointMR: A joint likelihood-based approach for causal effect estimation using multiple GWAS summary databases in Mendelian randomization. Sijia Wu, Lei Hou, Hongkai Li, Fuzhong Xue  


Table 1. Code files and corresponding illustration
| Code file | Purpose | Corresponding output |
| :--- | :--- | :--- |
| primary-JointMR.R | Primary simulation under the NOME setting without horizontal pleiotropy. This script evaluates different numbers of exposure and outcome GWAS databases. | Numerical results for the primary simulations |
| primary-plot-tau-compare.R | Plotting script comparing fixed-effect and random-effect settings in the primary simulations. | Figures 2 and 3, Figures S1 |
| primary-plot-datanum-compare.R | Plotting script for the database-number comparison in the primary simulations. | Figures S2-S4 |
| all-snp-compare.R | Sensitivity simulation in which all methods use the same valid SNP set. | Numerical results for the harmonized-SNP sensitivity analysis (Figure S5). |
| all-snp-compare-plot.R | Plotting script for the harmonized-SNP sensitivity analysis. | Figure S5 |
| mvmeta-JointMR-compare.R | Simulation comparing JointMR with estimate-level multivariate meta-analysis. | Numerical results for the multivariate meta-analysis comparison (Figure S6). |
| mvmeta-JointMR-compare-plot.R | Plotting script for the multivariate meta-analysis comparison. | Figure S6 |
| relaxF-JointMR.R | Sensitivity simulation using a relaxed instrument-strength rule. SNPs are retained if they satisfy the F-statistic criterion in at least one exposure database. | Numerical results for the relaxed-F sensitivity analysis (Figure S7-S9). |
| relaxF-JointMR-plot.R | Plotting script for the relaxed-F sensitivity analysis. | Figure S7-S9 |
| ple-plugin-JointMR.R | Simulation under directional horizontal pleiotropy using the MR-Egger intercept-plugin correction. | Numerical results for the pleiotropy-adjusted JointMR simulations (Figure 4 and Figure S10). |
| ple-plugin-JointMR-plot.R | Plotting script for the horizontal-pleiotropy simulation. | Figure 4 and Figure S10 |
| relaxed-nome-ple-JointMR.R | Simulation under the relaxed-NOME setting, allowing uncertainty in SNP-exposure effects and additional covariance components induced by exposure-exposure and exposure-outcome correlations. | Main and supplementary plots for the relaxed-NOME analyses (Figure 5-6 and Figure S11-13). |
| relaxed-nome-ple-JointMR-plot.R | Plotting script for the relaxed-NOME simulations. | Figure 5-6 and Figure S11-13 |
| application-final.Rapplication-run.R | Real-data application to blood lipid traits and type 2 diabetes. | Application estimates, heterogeneity diagnostics, pleiotropy diagnostics and relaxed-NOME estimates. |

 Table 2. Required input data format for the real-data application
| Field | Description | Use in the workflow |
| :--- | :--- | :--- |
| SNP | Variant identifier, usually rsID. | Matching variants across exposure and outcome datasets. |
| Beta | Estimated SNP-trait association. | Construction of Wald-ratio estimates and MR estimators. |
| Se | Standard error of the SNP-trait association. | F-statistic calculation, variance estimation and weighting. |
| P | SNP-trait association P value. | Genome-wide significance filtering and instrument selection. |
| Effect allele | Effect allele corresponding to beta. | Allele harmonization across datasets. |
| Other allele | Non-effect allele. | Allele harmonization across datasets. |
| Samplesize | GWAS sample size. | Covariance calculation in the relaxed-NOME and application analyses. |
| Dataset name | Identifier of each exposure or outcome GWAS dataset. | Matching sample sizes, correlation matrices and output annotations. |
