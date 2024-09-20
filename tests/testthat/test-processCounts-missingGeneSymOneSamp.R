test_that("for a single sample represented by hgnc symbols, missing genes are identified and reported", {
  
  data("rawCounts")
  data("modelGeneId")
  source(system.file("scripts", "rawCountsSingSampProcess.R", package = "IdentifiHR"))
  rawCountsSingSamp <- rawCountsSingSampProcess(rawCounts)
  source(system.file("scripts", "rawCountsSymProcess.R", package = "IdentifiHR"))
  rawCountsSym <- rawCountsSymProcess(rawCountsSingSamp, modelGeneId)
  
  expect_warning(processCounts(y = rawCountsSym, geneIds = "HGNC"))
  
})
