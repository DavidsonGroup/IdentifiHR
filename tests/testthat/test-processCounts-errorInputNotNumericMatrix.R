test_that("call error message if the input of processCounts is not numeric matrix", {
  
  data("rawCounts")
  rawCounts <- as.character(rawCounts)
  expect_error(processCounts(rawCounts, geneIds = "ENSEMBL"))
  
})