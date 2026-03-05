# Internal data objects stored in R/sysdata.rda
# Author: Ashley L Weir <weir.a@wehi.edu.au>

# modelGeneId
# A data frame with 2604 rows and 4 columns:
#   ensembl_id     - character vector of ensembl identifiers
#   betaCoef       - numeric vector of weighted genes in the trained IdentifiHR model
#   hgnc_symbol    - character vector of hgnc symbols
#   entrezgene_id  - numeric vector of entrez identifiers
# Source: generated from top table of DE analysis.

# modelMeanGenesIdentifiHR
# A data frame with 2604 rows and 1 column, with ensembl identifiers as rownames:
#   mean - numeric vector of mean log2 counts-per-million values for genes
#          required by the IdentifiHR model
# Source: created from training cohort for z score scaling.

# modelSDGenesIdentifiHR
# A data frame with 2604 rows and 1 column, with ensembl identifiers as rownames:
#   sd - numeric vector of standard deviations of log2 counts-per-million values
#        for genes required by the IdentifiHR model
# Source: created from training cohort for z score scaling.

# modelIdentifiHR
# An object of class "glmnet".
# Trained IdentifiHR elastic net penalised logistic regression model that
# predicts HR status using the weighted expression of 209 genes.
# Source: trained IdentifiHR model for predicting HR status.
