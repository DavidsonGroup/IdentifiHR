test_that("plotProb object is a list with two elements", {
  
  data(rawCounts)
  # Removing a model gene gene to give an example of the plotMissingness() function.
  rawCountsMissing <- rawCounts[rownames(rawCounts) != "ENSG00000160959", ]
  processedCounts <- processCounts(y = rawCountsMissing, geneIds = "ENSEMBL")
  missingGenes <- interrogateMissingness(y = rawCountsMissing, geneIds = "ENSEMBL")
  plottedMissing <- plotMissingness(missingGenes)
  
  expect_length(plottedMissing, 2)
  
  })
