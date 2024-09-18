# inst/scripts/rawCountsEntProcess
# code to change rawCounts gene identifiers from ensembl to entrez identifiers

rawCountsEntProcess <- function(rawCounts, modelGeneId) {
  
  # load dplyr as required
  library(dplyr)
  
  # add ensembl_id column to rawCounts for left_join
  rawCounts$ensembl_id <- rownames(rawCounts)
  
  # left_join with modelGeneId on ensembl_id
  rawCountsEnt <- left_join(rawCounts, modelGeneId, by = "ensembl_id")
  
  # remove rows where entrezgene_id is NA
  rawCountsEnt <- rawCountsEnt |>
    dplyr::filter(!is.na(entrezgene_id))
  
  # reset row names of counts to be enterezgene_id
  rownames(rawCountsEnt) <- rawCountsEnt$entrezgene_id
  
  # remove other gene identifier columns
  rawCountsEnt <- rawCountsEnt |>
    dplyr::select(-ensembl_id, -hgnc_symbol, -entrezgene_id, -betaCoef)
  
  # return rawCounts with entrez identifiers
  return(rawCountsEnt)
  
}

# example:
# rawCountsEnt <- rawCountsEntProcess(rawCounts, modelGeneId)
