#' IdentifiHR: means of 2604 DE genes
#'
#' Mean expression values for model genes, to be used in z score scaling.
#'
#' @format A data frame with 2604 rows and 1 column, with ensembl identifiers given as rownames:
#' \describe{
#'   \item{mean}{A numeric vector of mean log2 counts-per-million values for genes required by the identifiHR model}
#' }
#' @source Created from training cohort for z score scaling.
#' 
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
"modelMeanGenesIdentifiHR"