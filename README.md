# BootStratMI
### **Boot**strap timepoint-**Strat**ified **M**utual **I**nformation

A statistical pipeline using R to evaluate differential network analysis among time series proteomics assays. Particularly suited for MaxQuant output files.

## Overview
This code aims to uncover differential changes in protein association between treatment groups, specifically under a time-series experimental design. 
This is accomplished by using minet to calculate the mutual information between proteins of MaxQuant output data. For statistical significance, results are bootstrapped 5000 times by random sampling with replacement from observations **from the same timepoint**. P-values are then calculated empirically by how often bootstrapped results of one treatment group fall above or below the observed value in the alternative treatment group. This is done reciprocally, so bootstrapped treatment samples are tested against the observed control values, and bootstrapped control is tested against the observed treatment values. Storey q-values are implemented to correct for multiple hypothesis testing. 

## Installation and Setup

### 1. Install R
This pipeline was developed and tested on **R 4.5.2**. Download R from
[CRAN](https://cran.r-project.org/). (RStudio is optional but recommended.)

### 2. Install the required packages
Most dependencies are on CRAN; `minet` and `qvalue` come from Bioconductor.

```r
# CRAN packages
install.packages(c(
  "tibble", "dplyr", "tidyr", "purrr",
  "ggplot2", "stringr", "doSNOW", "foreach"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("minet", "qvalue"))
```

`parallel` ships with base R, so no separate installation is needed.

### 3. Get the code
```bash
git clone https://github.com/bsaiz123/BootStratMI.git
cd BootStratMI
```

### 4. Load the functions
There is no package to install — source the two function files at the start
of your R session:

```r
source("MQ-DataProcessingFunctions.R")
source("MutualInformationFunctions.R")
```

You're now ready to run the pipeline (see **Usage**).

## Prerequisite R Packages
The following R packages required for this workflow and their tested versions are:
### CRAN
- tibble (v 3.3.0)
- dplyr (v 1.1.4)
- tidyr (v 1.3.1)
- purrr (v 1.2.0)
- ggplot2 (v 4.0.1)
- stringr (v 1.6.0)
- doSNOW (v 1.0.20)
- foreach (v 1.5.2)
### Base R
- parallel (v 4.5.2)
### Bioconductor
- minet (v 3.66.0)
- qvalue (v 3.21)

## Input Data Format
- This work is intended to be used with the proteinGroups.txt output file generated from MaxQuant LFQ analysis

## Usage
1. Load the pipeline functions
    ```r
    source("MQ-DataProcessingFunctions.R")
    source("MutualInformationFunctions.R")
    ```
2. Import the MaxQuant proteinGroups output
    ```r
    proteinGroups <- read.delim("proteinGroups.txt", stringsAsFactors = FALSE)
    ```
3. Define Treatment groups (string patterns that appear in LFQ columns)
    - For example, if groups are DD and WW the code would read:
     ```r
     treatments <- c("DD","WW")
     ```
4. Data Preparation
    - Here the MaxQuant Data will undergo standard filtration steps and transformation, and will also be filtered by a user-specified valid-values percentage, evaluated per group (group completeness). Final results will impute values for remaining NAs
      ```r
      preparedData <- RawDataPreparation(
        df                = proteinGroups, # Original MQ Data
        treatments        = treatments, # Previously defined treatment groups
        dataColPattern    = "^LFQ.intensity", # String pattern to select LFQ intensity columns
        majorityProteinCol = "Majority.protein.IDs", # Column used for Protein Group ID
        NA.Threshold      = 0.7 # filter for 70% valid values after Log2 transformation
      )
      ```
5. Metadata Guide Construction
    - This creates a dataframe which guides the downstream process
      ```r
      metadata <- ExtractMetadata(
        df         = proteinGroups, # Original MQ Data
        treatments = treatments, # Previously defined treatment groups
        timepoints = c("T0", "T1", "T2") # Timepoint string used in MQ Column
      )
      ```
6. Run the bootstrapping analysis
     ```r
    results <- BootMI(
     data            = preparedData,
     metadata        = metadata,
     treatmentColumn = "treatment",
     n_bootstrap     = 5000,
     treatments      = treatments,
     estimator       = "spearman",
     timeStrat       = TRUE,   # stratify resampling within timepoints
     parallel        = TRUE,   # FALSE to run single-threaded
     workers         = 0       # Cores for parallel processing; 0 = detectCores() - 1
    )
     ```
7. Categorize protein pairs based on statistical significance (q-value) and effect size (observed change in MI)
     ```r
      results <- CategorizeResults(results, qThreshold = 0.05, effectThreshold = 0.01)
     ```
### Categorization Output Columns
- **Protein 1**: The first protein group of the pair used in MI analysis
- **Protein 2**: The second protein group of the pair used in MI analysis
- **OBS_MI1**: the observed mutual information from the first experimental group
- **OBS_MI2**: the observed mutual information from the second experimental group
- **OBS_Dif**: the observed difference in mutual information between experimental groups
- **pval_1vs2**: Empirically derived statistical significance from the first group's observed value against the second group's bootstrapped distribution
- **pval_2vs1**: Empirically derived statistical significance from the second group's observed value against the first group's bootstrapped distribution
- **Qval_1vs2**: Multiple testing corrected q-values derived from pval_1vs2
- **Qval_2vs1**: Multiple testing corrected q-values derived from pval_2vs1
- **category**: The observed change based on statistical significance (Significant Up, Significant Down, Non-significant)

## Repository Structure
- **MQ-DataProcessingFunctions.R**
  - This R file contains the library imports as well as the functions involved in data formatting and preparation. Used to set up for downstream MI analysis.
-  **MutualInformationFunctions.R**
  - This R file contains the functions used to carry out the statistical analysis of treatment groups.
- **LICENSE**
  - GPL-3.0 license text.
- **README.md**
  - (This file)

## License
This project is licensed under GPL-3.0. See LICENSE.

## Contact
Email regarding questions and comments at brandonsaiz.mail@gmail.com
