
testthat::test_that("processed counts are numeric and not NA", {
  
  data("rawCounts")
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  expect_true(is.matrix(processedCounts))
  expect_true(is.numeric(processedCounts))
  
})


