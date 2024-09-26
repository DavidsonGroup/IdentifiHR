# inst/scripts/rawCountsSummarisedExperiment
# code to change rawCounts matrix to a SummarisedExperiment object

rawCountsSummarisedExperiment <- function(rawCounts) {
  
  library(SummarizedExperiment)
  
  # convert rawCounts data frame to a SummarisedExperiment object
  rawCountsSe <- SummarizedExperiment(assays = list(counts = rawCounts))
  
  # return rawCounts as a SummarisedExperiment object
  return(rawCountsSe)

}

# example:
# rawCountsSe <- rawCountsSummarisedExperiment(rawCounts)