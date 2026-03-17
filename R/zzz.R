# zzz.R

# Import 'is' function from methods
#' @importFrom methods is
#' @importFrom glue glue
NULL

# Declare global variables for NSE
utils::globalVariables(c(
  "ensembl_id", "entrezgene_id", "hgnc_symbol",
  "inputStatus", "modelGeneId", 
  "modelMeanGenesIdentifiHR", "modelSDGenesIdentifiHR",
  "normWeightedGene"
))