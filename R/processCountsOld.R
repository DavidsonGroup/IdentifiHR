#' IdentifiHR
#'
#' The "processCounts" function prepares gene expression counts for input into the IdentifiHR classifier. 
#'
#' @param y A numeric matrix of gene expression counts, with genes presented in rows and samples presented in columns.
#' @param logCpm Logical. Have the raw gene expression counts been transformed into log2-counts-per-million? Default is TRUE, being that the counts matrix has been transformed. If FALSE, a log2-counts-per-million transformation will be applied.
#' @param scaled Logical. Has z-score scaling been performed? Default is TRUE, being that the counts matrix has been scaled. If FALSE, z-score scaling of each gene will be performed.
#' @param geneIds How are genes annotated? Specify either "ENSEMBL", "SYMBOL" or "ENTREZID".
#' @return A numeric matrix of counts for model genes, that has been log2-counts-per-million transformed and z-score scaled across genes.
#' @import dplyr
#' @importFrom edgeR cpm
#' @export
#'
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}

processCountsOld <- function(y,
                          logCpm,
                          scaled,
                          geneIds) {
  
  if(geneIds == "ENSEMBL"){
    if(logCpm == "TRUE") {
      if(scaled == "TRUE") {
        print("Counts have already been logCpm transformed")
        print("Predicting HR status using scaled gene expression counts")
      }
    }
    
    if(logCpm == "TRUE") {
      if(scaled == "FALSE") {
        print("Counts have already been logCpm transformed")
        print("Z-score scaling gene expression counts and then predicting HR status using scaled gene expression counts")
        # Z-score scale each row of the counts matrix
        y <- scale(y,
                   center = TRUE,
                   scale = TRUE)
      }
    }
    
    if(logCpm == "FALSE") {
      if(scaled == "FALSE") {
        print("Counts are being logCpm transformed")
        y <- edgeR::cpm(y, log = TRUE)
        print("Z-score scaling gene expression counts and then predicting HR status using scaled gene expression counts")
        # Z-score scale each row of the counts matrix
        y <- scale(y,
                   center = TRUE,
                   scale = TRUE)
      }
    }
    
    if(logCpm == "FALSE") {
      if(scaled == "TRUE") {
        print("Cannot transform data if it has already been scaled. Please use un-scaled counts data as an input.")
      }
    }
    
    # Strip any potential ensembl version numbers
    rownames(y) <- gsub("\\..*","", rownames(y))
    # Subset to differentially expressed genes identified in model training (using ensembl IDs)
    y <- y[rownames(y) %in% trainingCohortDeEns, ]
    print("Subsetting to only model genes.")
    return(y)
  }
  
  if(geneIds == "SYMBOL"){
    if(logCpm == "TRUE") {
      if(scaled == "TRUE") {
        print("Counts have already been logCpm transformed")
        print("Predicting HR status using scaled gene expression counts")
      }
    }
    
    if(logCpm == "TRUE") {
      if(scaled == "FALSE") {
        print("Counts have already been logCpm transformed")
        print("Z-score scaling gene expression counts and then predicting HR status using scaled gene expression counts")
        # Z-score scale each row of the counts matrix
        y <- scale(y,
                   center = TRUE,
                   scale = TRUE)
      }
    }
    
    if(logCpm == "FALSE") {
      if(scaled == "FALSE") {
        print("Counts are being logCpm transformed")
        y <- edgeR::cpm(y, log = TRUE)
        print("Z-score scaling gene expression counts and then predicting HR status using scaled gene expression counts")
        # Z-score scale each row of the counts matrix
        y <- scale(y,
                   center = TRUE,
                   scale = TRUE)
      }
    }
    
    if(logCpm == "FALSE") {
      if(scaled == "TRUE") {
        print("Cannot transform data if it has already been scaled. Please use un-scaled counts data as an input.")
      }
    }
    
    # Subset to differentially expressed genes identified in model training (using gene symbols)
    y <- y[rownames(y) %in% trainingCohortDeSym, ]
    return(y)
  }
  
  if(geneIds == "ENTREZ"){
    if(logCpm == "TRUE") {
      if(scaled == "TRUE") {
        print("Counts have already been logCpm transformed")
        print("Predicting HR status using scaled gene expression counts")
      }
    }
    
    if(logCpm == "TRUE") {
      if(scaled == "FALSE") {
        print("Counts have already been logCpm transformed")
        print("Z-score scaling gene expression counts and then predicting HR status using scaled gene expression counts")
        # Z-score scale each row of the counts matrix
        y <- scale(y,
                   center = TRUE,
                   scale = TRUE)
      }
    }
    
    if(logCpm == "FALSE") {
      if(scaled == "FALSE") {
        print("Counts are being logCpm transformed")
        y <- edgeR::cpm(y, log = TRUE)
        print("Z-score scaling gene expression counts and then predicting HR status using scaled gene expression counts")
        # Z-score scale each row of the counts matrix
        y <- scale(y,
                   center = TRUE,
                   scale = TRUE)
      }
    }
    
    if(logCpm == "FALSE") {
      if(scaled == "TRUE") {
        print("Cannot transform data if it has already been scaled. Please use un-scaled counts data as an input.")
      }
    }
    
    # Subset to differentially expressed genes identified in model training (using entrez IDs)
    y <- y[rownames(y) %in% trainingCohortDeEnt, ]
    return(y)
  }
  
}
