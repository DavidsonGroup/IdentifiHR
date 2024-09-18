test_that("plotProb object is a list with numeric elements", {
 
  data(rawCounts)
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  predictions <- predictHr(processedCounts)
  plottedProb <- plotProb(predictions)
  
  expect_true(is(plottedProb, "list"))
  expect_true(is.numeric(plottedProb[[1]][[1]]))
  expect_true(is.numeric(plottedProb[[1]][[2]]))
  expect_true(is.numeric(plottedProb[[1]][[3]]))
  expect_true(is.numeric(plottedProb[[1]][[4]]))
  expect_true(is.numeric(plottedProb[[2]][[1]]))
  expect_true(is.numeric(plottedProb[[2]][[2]]))
  
})
