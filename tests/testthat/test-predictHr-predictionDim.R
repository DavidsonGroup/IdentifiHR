test_that("predictHr produces a data frame with predictions for the correct number of samples", {
  
  data("rawCounts")
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  predicitons <- predictHr(processedCounts)
  
  expect_equal(nrow(predicitons), ncol(rawCounts))
  expect_equal(nrow(predicitons), ncol(processedCounts))
  
})
