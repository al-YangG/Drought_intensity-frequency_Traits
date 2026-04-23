# DIF-Traits: Drought intensity and frequency effects on plant traits

Data and R scripts from a factorial drought experiment investigating the effects of drought intensity and frequency on plant functional traits and fitness.

This repository contains a complete and reproducible workflow, including model selection, statistical analyses, visualization, and multivariate approaches.

## Repository structure

- code/ — R scripts for all analyses  
- data/ — raw data and metadata  
- data_processed/ — intermediate datasets  
- results/ — figures and tables  

## Code overview

All scripts are numbered to reflect the analysis workflow:

1. Model selection (AIC/BIC comparison)  
2. ANOVA and variance partitioning  
3. Variance heatmaps  
4. Trait responses (standardized)  
5. Main effects (significant effects only)  
6. Interactive effects (EMMs and contrasts)  
7. Fitness traits (sprouting, survival, reproduction)  
8. Fitness plots  
9. Trait correlations and network analysis  
10. PCA and RDA analyses  
11. PCA visualization  

## Data

- DIF-Traits_data.csv — main dataset (traits, treatments, fitness)  
- Metadata.pdf — description of variables and experimental design  

## Results

### Figures
Located in results/figures/, including:

- Trait responses and main effects  
- Interaction plots  
- Fitness responses  
- Variance partitioning heatmaps  
- Correlation matrix and trait networks  
- PCA (all traits and by category)

### Tables
Located in results/tables/, including:

- Model selection and ANOVA results  
- Fitness model outputs
- Trait correlations and network metrics
- PCA scores and trait loadings  
- RDA results  

## Requirements

Main R packages used:

tidyverse, dplyr, ggplot2, glmmTMB, emmeans, car, DHARMa, performance, vegan, factoextra, Hmisc, ggraph, tidygraph, ggpubr, ggrepel  

## Workflow

Run scripts sequentially from 01 to 11.  
Each script reads from data/ and writes outputs to results/.

## License

See LICENSE file for details.

## Author
Yang [yangg@natur.cuni.cz]
April 2026
