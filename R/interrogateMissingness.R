#' Interrogate which model genes are missing from input data
#'
#' The "interrogateMissingness" function allows you to examine which genes, if any, required for the IdentifiHR model are missing in your dataset. There are 2604 genes required for normalising library size differences and 209 gene required to predict HR status by IdentifiHR. 
#'
#' @param y A numeric matrix of raw gene expression counts, with genes presented in rows and samples presented in columns.
#' @param geneIds How are genes annotated? Specify either "ENSEMBL", "HGNCSYMBOL" or "ENTREZID" ("ENSEMBL" is preferred).
#' @return A plot of missing genes and their beta coefficient, or "weight" in the IdentifiHR model
#' @import dplyr
#' @importFrom edgeR cpm
#' @export
#'
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
#'
#' @examples
#' # to add

interrogateMissingness <- function(y,
                                   geneIds) {
  if (!is.matrix(y)) {
    
    stop("Input is not a matrix. The input must be a numeric matrix, with genes presented in rownames and samples presented in columns.")
    
  }
  
  if (!is.numeric(y)) {
    
    stop("Matrix is not numeric. The input must be a numeric matrix, with genes presented in rownames and samples presented in columns.")
    
  }
  
  if (geneIds == "ENSEMBL") {
    
    
    
  }
  
}