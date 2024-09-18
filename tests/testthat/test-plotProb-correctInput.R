test_that("check input is correct for plotting probability", {
  
  data(rawCounts)
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  predictions <- predictHr(processedCounts)
  predictions$cohort <- "input"
  joinCohort <- rbind(predictions, tcgaTestPred)
  colNamesPlot <- c("Sample", "hrPrediction", "predictionProbability", "cohort")
  
  expect_identical(colNamesPlot, colnames(joinCohort))


})
