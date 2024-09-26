test_that("SummarisedExperiment object are suitable inputs", {
  
  data("rawCounts")

  source(system.file("scripts", "rawCountsSummarisedExperiment.R", package = "IdentifiHR"))
  rawCountsSe <- rawCountsSummarisedExperiment(rawCounts)
  processedCounts <- processCounts(y = rawCountsSe, geneIds = "ENSEMBL")
  expect_true(is.matrix(processedCounts))
  expect_true(is.numeric(processedCounts))
  
})
