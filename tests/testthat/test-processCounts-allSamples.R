testthat::test_that("processed counts are present for all samples", {
  
  data("rawCounts")
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  expect_equal(ncol(rawCounts), ncol(processedCounts))
  
})