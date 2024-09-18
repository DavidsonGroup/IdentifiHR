test_that("predictHr produces a data frame with predictions for the coreect number of samples", {
  
  data("rawCounts")
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  predicitons <- predictHr(processedCounts)
  
  expect_equal(nrow(predicitons), ncol(rawCounts))
  expect_equal(nrow(predicitons), ncol(processedCounts))
  
})
