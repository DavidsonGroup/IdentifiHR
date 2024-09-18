test_that("call error message if the input of processCounts is not a data frame or a matrix", {
  
  data("rawCounts")
  rawCounts <- as.list(rawCounts)
  expect_error(processCounts(rawCounts, geneIds = "ENSEMBL"))
  
})