# Obesity and Cancer Risk – Meta-analyses and Figures

Code to run meta-analysis for the paper "Adiposity and cancer: systematic review and meta-analysis" based on the SupplementaryData.xlsx file, which is publicly available alongside the manuscript

This repository contains R scripts used to conduct meta-analyses and generate figures for a research project examining the association between adiposity measures and cancer risk, using observational epidemiology and Mendelian randomization approaches.

The code supports primary analyses, sensitivity analyses, and manuscript figures.

## Repository Contents

### Meta-analysis scripts
- `primary_obs_meta_analysis.R`  
  Primary observational meta-analyses across cancer sites

- `compare_waist_bmi.R`  
  Comparative analyses of BMI versus waist-based adiposity measures

- `MR_metaanalysis.R`  
  Mendelian randomization meta-analyses

### Figure scripts
- `Fig1.R` – `Fig7.R`  
  Scripts used to generate manuscript figures and selected supplementary figures

Each figure script is intended to be run independently once the required input data are available.

## Software Requirements

- R (version 4.4.1)
