test_that("for a single sample represented by ensembl ids, missing genes are identified and reported", {

  data("rawCounts")
  source(system.file("scripts", "rawCountsSingSampProcess.R", package = "IdentifiHR"))
  rawCountsSingSamp <- rawCountsSingSampProcess(rawCounts)
  rawCountsSingSamp <- rawCountsSingSamp |>
    rownames_to_column(var = "ensembl_id") |>
    dplyr::filter(ensembl_id != "ENSG00000160959") |>
    column_to_rownames(var = "ensembl_id")

  expect_warning(processCounts(y = rawCountsSingSamp, geneIds = "ENSEMBL"))
  
  })

