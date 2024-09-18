testthat::test_that("processed counts contain all required genes, with hgnc symbols as input, for only a single sample", {
  
  data("rawCounts")
  data("modelGeneId")
  source(system.file("scripts", "rawCountsSingSampProcess.R", package = "IdentifiHR"))
  rawCountsSingSamp <- rawCountsSingSampProcess(rawCounts)
  source(system.file("scripts", "rawCountsSymProcess.R", package = "IdentifiHR"))
  rawCountsSym <- rawCountsSymProcess(rawCountsSingSamp, modelGeneId)
  processedCounts <- processCounts(y = rawCountsSym, geneIds = "HGNC")
  
  expect_identical(rownames(processedCounts), rownames(modelMeanGenesIdentifiHR))
  expect_identical(rownames(processedCounts), rownames(modelSDGenesIdentifiHR))
  expect_identical(rownames(processedCounts), modelGeneId$ensembl_id)
  
})
