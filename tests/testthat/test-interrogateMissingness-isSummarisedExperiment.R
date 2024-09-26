test_that("SummarisedExperiment object are suitable inputs", {
  
  data("rawCounts")
  
  source(system.file("scripts", "rawCountsSummarisedExperiment.R", package = "IdentifiHR"))
  
  rawCountsMissing <- rawCounts[rownames(rawCounts) != "ENSG00000160959", ]
  rawCountsSe <- rawCountsSummarisedExperiment(rawCountsMissing)
  
  colNamesMiss <- c("ensembl_id", "hgnc_symbol", "entrezgene_id", "betaCoef", "inputStatus",  "normWeightedGene")
  inputStatusVec <- c("present", "missing")
  missingGenes <- interrogateMissingness(y = rawCountsSe, geneIds = "ENSEMBL")
  
  expect_true(is.data.frame(missingGenes))
  expect_equal(ncol(missingGenes), 6)
  expect_identical(colnames(missingGenes), colNamesMiss)
  expect_in(missingGenes$inputStatus, inputStatusVec)
  
})