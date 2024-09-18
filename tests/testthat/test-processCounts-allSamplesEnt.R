testthat::test_that("processed counts are present for all samples, with entrez identifiers as input", {
  
  data("rawCounts")
  data("modelGeneId")
  source(system.file("scripts", "rawCountsEntProcess.R", package = "IdentifiHR"))
  rawCountsEnt <- rawCountsEntProcess(rawCounts, modelGeneId)
  processedCounts <- processCounts(y = rawCountsEnt, geneIds = "ENTREZ")
  expect_equal(ncol(rawCounts), ncol(processedCounts))
  
})
