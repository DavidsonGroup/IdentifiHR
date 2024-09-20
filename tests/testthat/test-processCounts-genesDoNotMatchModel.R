test_that("when all required genes are present in input, they do not match the model genes", {
  
  data("rawCounts")
  # Removing a model gene gene to give an example of non-matching.
  rawCountsMissing <- rawCounts[rownames(rawCounts) != "ENSG00000160959", ]
  y <- rawCountsMissing
  # Removing a model gene gene to give a warning message.
  data("modelGeneId")
  geneId <- modelGeneId$ensembl_id
  # Strip any potential ensembl version numbers
  rownames(y) <- gsub("\\..*","", rownames(y))
  
  # Subset to differentially expressed genes identified in model training (using ensembl IDs)
  y <- y[rownames(y) %in% geneId, ]
  
  expect_true(identical(sort(rownames(y)), sort(geneId)) == FALSE)
  
})
