test_that("plotMissingness object is a list with numeric elements", {
  
  data(rawCounts)
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  missingGenes <- interrogateMissingness(y = rawCounts, geneIds = "ENSEMBL")
  plottedMissing <- plotMissingness(missingGenes)
  
  expect_true(is(plottedMissing, "list"))
  expect_true(is.numeric(plottedMissing[[1]][[1]]))
  expect_true(is.numeric(plottedMissing[[1]][[2]]))
  expect_true(is.numeric(plottedMissing[[1]][[3]]))
  expect_true(is.numeric(plottedMissing[[1]][[4]]))
  expect_true(is.numeric(plottedMissing[[2]][[1]]))
  expect_true(is.numeric(plottedMissing[[2]][[2]]))
  
})
