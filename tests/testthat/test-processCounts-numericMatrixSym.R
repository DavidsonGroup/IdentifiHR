testthat::test_that("processed counts are numeric and not NA, with hgnc symbols as input", {
  
  data("rawCounts")
  modelGeneId <- IdentifiHR:::modelGeneId
  source(system.file("scripts", "rawCountsSymProcess.R", package = "IdentifiHR"))
  rawCountsSym <- rawCountsSymProcess(rawCounts, modelGeneId)
  processedCounts <- processCounts(y = rawCountsSym, geneIds = "HGNC")
  expect_true(is.matrix(processedCounts))
  expect_true(is.numeric(processedCounts))
  
})