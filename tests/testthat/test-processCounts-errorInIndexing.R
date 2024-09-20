test_that("error is called if models genes do not match input after indexing", {
  
  data("rawCounts")
  data("modelGeneId")
  data("modelMeanGenesIdentifiHR")
  geneId <- modelGeneId$ensembl_id
  rownames(rawCounts) <- gsub("\\..*","", rownames(rawCounts))
  rawCounts <- rawCounts[rownames(rawCounts) %in% geneId, ]
  completeCounts <- rawCounts[geneId, ]
  completeCounts <- completeCounts[-1, ]
  
  geneIdx <- function(modelMeanGenesIdentifiHR, completeCounts) { 
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(completeCounts)) == TRUE) {
      
      # Transform to log2 counts-per-million
      countsCpm <- edgeR::cpm(completeCounts, log = TRUE) # Rows give genes, columns give samples
      message("Normalising sequencing depth.")
      # Z score scale counts using vectors saved from model training
      message("Z-score scaling genes using mean and sd from training cohort.")
      countsCpmZ <- (countsCpm - meanGene) / sdGene
      message("Counts successfully processed for IdentifiHR.")
      return(countsCpmZ)
      
    } else if (identical(rownames(modelMeanGenesIdentifiHR), rownames(completeCounts)) == FALSE) {
      
      stop("Error in gene indexing. Please check input matrix and/or notify package author.")
      
    }
    
  }
  
  expect_error(geneIdx(modelMeanGenesIdentifiHR, completeCounts),
               regexp = "Error in gene indexing. Please check input matrix and/or notify package author.") 
  
  expect_error(is.matrix(geneIdx(modelMeanGenesIdentifiHR, completeCounts)))
  
})
