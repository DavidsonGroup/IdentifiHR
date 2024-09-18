test_that("call error message if the input of ProcessCounts is not numeric", {

  data("rawCounts")
  rawCounts[2,3] <- "string"
  expect_error(processCounts(rawCounts, geneIds = "ENSEMBL"))
  
})
