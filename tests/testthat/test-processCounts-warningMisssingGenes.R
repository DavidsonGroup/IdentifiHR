test_that("a warning message is generated when genes are missing", {
  
  data("rawCounts")
  # Removing a model gene gene to give a warning message.
  rawCountsMissing <- rawCounts[rownames(rawCounts) != "ENSG00000160959", ]
  data("modelGeneId")
  geneId <- modelGeneId$ensembl_id

  expect_warning(processCounts(y = rawCountsMissing, geneIds = "ENSEMBL"))
  
})
