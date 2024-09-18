# inst/scripts/rawCountsSymProcess
# code to change rawCounts gene identifiers from ensembl to hgnc symbols

rawCountsSymProcess <- function(rawCounts, modelGeneId) {
  
  # load dplyr as required
  library(dplyr)
  
  # add ensembl_id column to rawCounts for left_join
  rawCounts$ensembl_id <- rownames(rawCounts)
  
  # left_join with modelGeneId on ensembl_id
  rawCountsEnt <- left_join(rawCounts, modelGeneId, by = "ensembl_id")
  
  # remove rows where hgnc_symbol is NA
  rawCountsSym <- rawCountsSym |>
    dplyr::filter(!is.na(hgnc_symbol))
  
  # reset row names of counts to be hgnc_symbol
  rownames(rawCountsSym) <- rawCountsSym$hgnc_symbol
  
  # remove other gene identifier columns
  rawCountsSym <- rawCountsSym |>
    dplyr::select(-ensembl_id, -hgnc_symbol, -entrezgene_id, -betaCoef)
  
  # return rawCounts with hgnc symbols
  return(rawCountsSym)
  
}

# example:
# rawCountsSym <- rawCountsSymProcess(rawCounts, modelGeneId)
  