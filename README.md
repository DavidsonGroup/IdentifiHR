# IdentifiHR
IdentifiHR is a predictive machine learning model of homologous recombination (HR) status in high-grade serous ovarian carcinoma (HGSC) that uses only gene expression.

## Background
Approximately half of all HGSCs have a therapeutically targetable defect in the HR DNA repair mechanism. HGSC is the most commonly HR deficient (HRD) cancer type, largely due to the frequency of germline and somatic mutations in the HR-genes, BRCA1/2. While there are genomic methods and some transcriptomic signatures, developed for alternate cancers, to identify HRD patients, there are no gene expression-based tools to predict HR repair status in HGSC specifically. We have built the first HGSC-specific model to predict HR repair status using gene expression.

## Installation 

To install the package, use the following script in R:

```
# install.packages("devtools")
devtools::install_github("DavidsonGroup/IdentifiHR")
```

## How does it work?

The IdentifiHR R package has several functions to support use, and can be used to predict HR status in a single sample, or across several samples. The model requires only a matrix of raw gene expression counts, with genes annotated with ENSEMBL identifiers.

The processCounts() function subsets the input matrix to only the genes required for predicition. It then transforms counts with log2 counts-per-million (CPM) to normalise for library size differences, and scales genes using a z-score, whereby the mean and standard deviation are taken from our training dataset.

Processed counts can then be used by the predictHR() function to estimate the probability that a sample is HRD. PredictHR() uses the trained IdentiiHR model to infer HR status from the expression of only 209 genes.

The output of IdentifiHR is a data frame containing both a discrete prediction of HR status, being HRD OR HR proficient (HRP), in addiiton to the probability that a sample is HRD.

