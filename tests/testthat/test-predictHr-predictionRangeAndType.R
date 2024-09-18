test_that("predictHR gives thge correct prediciton labels and outputs a probability within the correct range", {
 
  data("rawCounts")
  hrLevels <- c("HRD", "HRP")
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  predicitons <- predictHr(processedCounts)
  
  expect_in(predicitons$hrPrediction, hrLevels)
  expect_true(all(predicitons$predictionProbability <= 1 & predicitons$predictionProbability >= 0))
  
})
