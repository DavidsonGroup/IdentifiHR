testthat::test_that("processed counts contain all required genes", {
  
  data("rawCounts")
  modelGeneId <- IdentifiHR:::modelGeneId
  modelMeanGenesIdentifiHR <- IdentifiHR:::modelMeanGenesIdentifiHR
  modelSDGenesIdentifiHR <- IdentifiHR:::modelSDGenesIdentifiHR
  processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
  
  expect_identical(rownames(processedCounts), rownames(modelMeanGenesIdentifiHR))
  expect_identical(rownames(processedCounts), rownames(modelSDGenesIdentifiHR))
  expect_identical(rownames(processedCounts), modelGeneId$ensembl_id)
  
})