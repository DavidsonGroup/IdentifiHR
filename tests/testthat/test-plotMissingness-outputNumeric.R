test_that("plotMissingness object is a list with numeric elements", {
  
  data(rawCounts)
  # Removing a model gene gene to give an example of the plotMissingness() function.
  rawCountsMissing <- rawCounts[rownames(rawCounts) != "ENSG00000160959", ]
  processedCounts <- processCounts(y = rawCountsMissing, geneIds = "ENSEMBL")
  missingGenes <- interrogateMissingness(y = rawCountsMissing, geneIds = "ENSEMBL")
  plottedMissing <- plotMissingness(missingGenes)
  
  expect_true(is(plottedMissing, "list"))
  expect_true(is.numeric(plottedMissing[[1]][[1]]))
  expect_true(is.numeric(plottedMissing[[1]][[2]]))
  expect_true(is.numeric(plottedMissing[[1]][[3]]))
  expect_true(is.numeric(plottedMissing[[1]][[4]]))
  expect_true(is.numeric(plottedMissing[[2]][[1]]))
  expect_true(is.numeric(plottedMissing[[2]][[2]]))
  
})
