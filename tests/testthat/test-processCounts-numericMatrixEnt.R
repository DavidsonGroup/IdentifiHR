testthat::test_that("processed counts are numeric and not NA, with entrez identifiers as input", {
  
  data("rawCounts")
  modelGeneId <- IdentifiHR:::modelGeneId
  source(system.file("scripts", "rawCountsEntProcess.R", package = "IdentifiHR"))
  rawCountsEnt <- rawCountsEntProcess(rawCounts, modelGeneId)
  
  expect_warning(processedCounts <- processCounts(y = rawCountsEnt, geneIds = "ENTREZ"))
  expect_true(is.matrix(processedCounts))
  expect_true(is.numeric(processedCounts))
  
})
