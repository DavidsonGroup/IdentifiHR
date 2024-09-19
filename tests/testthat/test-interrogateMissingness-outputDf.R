test_that("interrogateMissingness object is a data frame with the correct columns", {
  
  data(rawCounts)
  # Removing a model gene gene to give an example of the interrogateMissingness() function.
  rawCountsMissing <- rawCounts[rownames(rawCounts) != "ENSG00000160959", ]
  colNamesMiss <- c("ensembl_id", "hgnc_symbol", "entrezgene_id", "betaCoef", "inputStatus",  "normWeightedGene")
  inputStatusVec <- c("present", "missing")
  missingGenes <- interrogateMissingness(y = rawCountsMissing, geneIds = "ENSEMBL")
  
  expect_true(is.data.frame(missingGenes))
  expect_equal(ncol(missingGenes), 6)
  expect_identical(colnames(missingGenes), colNamesMiss)
  expect_in(missingGenes$inputStatus, inputStatusVec)
  
})
