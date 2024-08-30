# IdentifiHR                                                                
<img src="https://github.com/user-attachments/assets/4f46868c-df3b-4525-8fa7-e049dd74508c" width="100" height="100" align="right">
IdentifiHR is a predictive machine learning model of homologous recombination (HR) status in high-grade serous ovarian carcinoma (HGSC) that uses only gene expression.


Please see the [wiki](https://github.com/DavidsonGroup/IdentifiHR/wiki) for installation and detailled usage documentation, in addititon to worked examples.

For more information about model training and testing, please read our [preprint](https://www.biorxiv.org/content/10.1101/2024.08.15.608185v1).

## Installation: 

To install the package, please use the following code in R:

```
library(devtools)   
devtools::install_github("DavidsonGroup/IdentifiHR", build_vignettes = TRUE)
library(IdentifiHR)
```
## How to predict HR status in HGSC quickly:

IdentifiHR can predict HR status in a single HGSC sample or across many samples. It requires only a matrix or data frame of raw gene expression counts, collected through bulk or pseudobulked single-cell RNA sequencing, with samples given as columns and genes given by rows (in rownames). To use the trained IdentifiHR model, counts must first be transformed and scaled by "processcounts()", after which "predictHr()" can be used to predict the discrete HR status, being HR deficient ("HRD") or proficient ("HRP"), in addition to the probability that a sample is HRD. For basic use, please use the following code in R:

```
processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
predictions <- predictHr(processedCounts)
```
Note: IdentifiHR uses the expression of 2604 genes for normalisation, and further, uses only 209 of these genes when predicting HR status in HGSC. All genes are needed to ensure optimal model accuracy, however, if a sample is "missing" some of the required genes, "processCounts()" will warn users and we recommended that the missing genes are investigated, to determine their use and contribution to IdentifiHR. The package includes functions to allow this missingness to be investigated, including "interrogateMissingness()" and "plotMissingness()".

## Package overview:

![identifiHRPackageOverview](https://github.com/user-attachments/assets/1ac74610-004a-44c1-824b-ed4af67f5994)
