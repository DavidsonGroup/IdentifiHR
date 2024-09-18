test_that("predictHr produces a data frame with the correct columns", {
 
  data("rawCounts")
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  predicitons <- predictHr(processedCounts)
  colNames <- c("Sample", "hrPrediction", "predictionProbability")
  
  expect_true(is.data.frame(predicitons))
  expect_equal(ncol(predicitons), 3)
  expect_identical(colnames(predicitons), colNames)
  
})
