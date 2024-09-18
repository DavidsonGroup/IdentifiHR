test_that("plotProb object is a list with two elements", {
  
  data(rawCounts)
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  missingGenes <- interrogateMissingness(y = rawCounts, geneIds = "ENSEMBL")
  plottedMissing <- plotMissingness(missingGenes)
  
  expect_length(plottedMissing, 2)
  
  })
