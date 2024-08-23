#' IdentifiHR
#'
#' The "interrogateMissingness" function allows you to examine which genes, if any, required for the IdentifiHR model are missing in your dataset. 
#'
#' @param y A numeric matrix of raw gene expression counts, with genes presented in rows and samples presented in columns.
#' @param geneIds How are genes annotated? Specify either "ENSEMBL", "HGNCSYMBOL" or "ENTREZID" ("ENSEMBL" is preferred).
#' @return A numeric matrix of counts for model genes, that has been log2-counts-per-million transformed and z-score scaled across genes.
#' @import dplyr
#' @importFrom edgeR cpm
#' @export
#'
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
#'
#' @examples
#' 
#' 