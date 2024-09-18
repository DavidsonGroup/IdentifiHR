# inst/scripts/rawCountsSingSampProcess
# code to change rawCounts gene identifiers from ensembl to entrez identifiers

rawCountsSingSampProcess <- function(rawCounts) {

  rawCounts <- rawCounts |>
    dplyr::select(Sample1) 
  # return rawCounts with entrez identifiers
  return(rawCounts)
  
}

# example:
# rawCountsSingSamp <- rawCountsSingSampProcess(rawCounts)