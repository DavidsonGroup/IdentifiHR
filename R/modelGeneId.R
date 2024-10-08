#' IdentifiHR: gene annotations and beta coefficients
#'
#' Identifiers for model genes
#'
#' @docType data
#' @format A data frame with 2604 rows and 3 columns:
#' \describe{
#'   \item{ensembl_id}{A character vector of ensembl identifiers}
#'   \item{betaCoef}{A numeric vector of weighted genes in the trained IdentifiHR model}
#'   \item{hgnc_symbol}{A character vector of hgnc symbols}
#'   \item{entrezgene_id}{A numeric vector of entrez identifiers}

#' }
#' 
#' @source Created for subsetting genes down to those required by identifiHR. Generated from top table of DE analysis.
#' 
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
"modelGeneId"