test_that("when all required genes are present in input", {
  
  data("rawCounts")
  
  expect_error(interrogateMissingness(y = rawCounts, geneIds = "ENSEMBL"))
  
})
