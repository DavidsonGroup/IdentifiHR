testthat::test_that("processed counts are present for all samples, with hgnc symbols as input", {
  
  data("rawCounts")
  data("modelGeneId")
  source(system.file("scripts", "rawCountsSymProcess.R", package = "IdentifiHR"))
  rawCountsSym <- rawCountsSymProcess(rawCounts, modelGeneId)
  processedCounts <- processCounts(y = rawCountsSym, geneIds = "HGNC")
  expect_equal(ncol(rawCounts), ncol(processedCounts))
  
})