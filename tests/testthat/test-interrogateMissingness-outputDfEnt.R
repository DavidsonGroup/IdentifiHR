test_that("interrogateMissingness object is a data frame with the correct columns, with entrez identifiers as input", {
  
  data("rawCounts")
  data("modelGeneId")
  source(system.file("scripts", "rawCountsEntProcess.R", package = "IdentifiHR"))
  rawCountsEnt <- rawCountsEntProcess(rawCounts, modelGeneId)
  colNamesMiss <- c("ensembl_id", "hgnc_symbol", "entrezgene_id", "betaCoef", "inputStatus",  "normWeightedGene")
  inputStatusVec <- c("present", "missing")
  missingGenes <- interrogateMissingness(y = rawCountsEnt, geneIds = "ENTREZ")
  
  expect_true(is.data.frame(missingGenes))
  expect_equal(ncol(missingGenes), 6)
  expect_identical(colnames(missingGenes), colNamesMiss)
  expect_in(missingGenes$inputStatus, inputStatusVec)
  
})

