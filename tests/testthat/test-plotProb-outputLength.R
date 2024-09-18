test_that("plotProb object is a list with two elements", {
  
  data(rawCounts)
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  predictions <- predictHr(processedCounts)
  plottedProb <- plotProb(predictions)
  
  expect_length(plottedProb, 2)
  
})
