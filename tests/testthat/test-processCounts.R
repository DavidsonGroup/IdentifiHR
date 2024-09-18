testthat::test_that("processed counts contain all required genes", {
  
  data("rawCounts")
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  
  expect_identical(rownames(processedCounts), rownames(modelMeanGenesIdentifiHR))
  expect_identical(rownames(processedCounts), rownames(modelSDGenesIdentifiHR))
  expect_identical(rownames(processedCounts), modelGeneId$ensembl_id)
  
})

testthat::test_that("processed counts are present for all samples", {
  
  data("rawCounts")
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  expect_equal(ncol(rawCounts), ncol(processedCounts))
  
})

testthat::test_that("processed counts gives a matrix", {
  
  data("rawCounts")
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  expect_s4_class(object = processedCounts, class = "matrix")
  
})

testthat::test_that("processed counts are numeric and not NA", {
  
})


