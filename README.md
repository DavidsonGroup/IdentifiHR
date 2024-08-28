# IdentifiHR                                                                
<img src="https://github.com/user-attachments/assets/4f46868c-df3b-4525-8fa7-e049dd74508c" width="100" height="100" align="right">
IdentifiHR is a predictive machine learning model of homologous recombination (HR) status in high-grade serous ovarian carcinoma (HGSC) that uses only gene expression.


Please see the [wiki](https://github.com/DavidsonGroup/IdentifiHR/wiki) for installation and usage documentation, in addititon to worked examples.

For more information about model training and testing, please read our [preprint](https://www.biorxiv.org/content/10.1101/2024.08.15.608185v1).

## Installation: 

To install the package, please use the following code in R:

```
# install.packages("devtools")
devtools::install_github("DavidsonGroup/IdentifiHR")
```
## How to predict HR status in HGSC quickly:

IdentifiHR can predict HR status in a single HGSC sample or across many samples. It requires only a matrix or data frame of raw gene expression counts, collected through bulk or pseudobulked single-cell RNA sequencing, with samples given as columns and genes given by rows (in rownames). To use the trained IdentifiHR model, counts must first be transformed and scaled by "processcounts()", after which "predictHr()" can be used to predict the discrete HR status, being HR deficient ("HRD") or proficient ("HRP"), in addition to the probability that a sample is HRD. For basic use, please use the following code in R:

```
processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
predictions <- predictHr(processedCounts)
```
