testthat::test_that("processed counts contain all required genes, when only a single sample is present", {
  
  data("rawCounts")
  source(system.file("scripts", "rawCountsSingSampProcess.R", package = "IdentifiHR"))
  rawCountsSingSamp <- rawCountsSingSampProcess(rawCounts)
  processedCounts <- processCounts(y = rawCountsSingSamp, geneIds = "ENSEMBL")
  
  expect_identical(rownames(processedCounts), rownames(modelMeanGenesIdentifiHR))
  expect_identical(rownames(processedCounts), rownames(modelSDGenesIdentifiHR))
  expect_identical(rownames(processedCounts), modelGeneId$ensembl_id)
  
})