#' IdentifiHR: Example raw counts matrix for 10 samples
#'
#' Example raw counts matrix for 10 samples, though note that the gene "ENSG00000107175" is missing from this data set and required for model predictions. This allows for the "interrogateMissingness()" and "plotMissingness()" functions to be explored.
#'
#' @docType data
#' @format A gene expression counts table with 52124 rows, containing ensembl gene identifiers as rownames, and 10 columns, representing distinct samples.
#' 
#' @source Created as an example dataset for package testing and use.
#' 
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
"rawCounts"