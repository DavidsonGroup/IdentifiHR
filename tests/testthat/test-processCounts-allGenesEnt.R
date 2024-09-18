testthat::test_that("processed counts contain all required genes, with entrez identifiers as input", {
  
  data("rawCounts")
  data("modelGeneId")
  source(system.file("scripts", "rawCountsEntProcess.R", package = "IdentifiHR"))
  rawCountsEnt <- rawCountsEntProcess(rawCounts, modelGeneId)
  processedCounts <- processCounts(y = rawCountsEnt, geneIds = "ENTREZ")
  
  expect_identical(rownames(processedCounts), rownames(modelMeanGenesIdentifiHR))
  expect_identical(rownames(processedCounts), rownames(modelSDGenesIdentifiHR))
  expect_identical(rownames(processedCounts), modelGeneId$ensembl_id)
  
})
