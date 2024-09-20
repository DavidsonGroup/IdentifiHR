test_that("for a single sample represented by entrez ids, missing genes are identified and reported", {
  
  data("rawCounts")
  data("modelGeneId")
  source(system.file("scripts", "rawCountsSingSampProcess.R", package = "IdentifiHR"))
  rawCountsSingSamp <- rawCountsSingSampProcess(rawCounts)
  source(system.file("scripts", "rawCountsEntProcess.R", package = "IdentifiHR"))
  rawCountsEnt <- rawCountsEntProcess(rawCountsSingSamp, modelGeneId)
  
  expect_warning(processCounts(y = rawCountsEnt, geneIds = "ENTREZ"))
  
})

