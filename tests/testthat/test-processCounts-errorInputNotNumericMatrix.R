test_that("call error message if the input of processCounts is not numeric matrix", {
  
  rawCounts <- as.data.frame(LETTERS)
  rawCounts <- as.matrix(rawCounts)
  expect_error(processCounts(rawCounts, geneIds = "ENSEMBL"))
  
})