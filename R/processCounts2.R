#' IdentifiHR
#'
#' The "processCounts" function prepares gene expression counts for input into the IdentifiHR classifier. 
#'
#' @param y A numeric matrix of raw gene expression counts, with genes presented in rownames and samples presented in columns.
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
processCounts <- function(y,
                          geneIds) {
  
  if (!is.matrix(y)) {
    stop("Input is not a matrix. The input must be a numeric matrix, with genes presented in rownames and samples presented in columns.")
  }
  
  if (!is.numeric(y)) {
    stop("Matrix is not numeric. The input must be a numeric matrix, with genes presented in rownames and samples presented in columns.")
  }
  
  if(geneIds == "ENSEMBL") {
    # Extract vector of ensembl identifiers required
    geneId <- modelGeneId$ENSEMBL
  
    # Strip any potential ensembl version numbers
    rownames(y) <- gsub("\\..*","", rownames(y))
    
    # Subset to differentially expressed genes identified in model training (using ensembl IDs)
    y <- y[rownames(y) %in% geneId, ]
    
    print("Subsetting to only model genes.")
    
    # Check all genes required are present
    countsT <- t(y)
    
    # If not all genes are present:
    if(identical(sort(colnames(countsT)), sort(geneId)) == FALSE){
      suppressWarnings(cat(paste0(((table(colnames(countsT) == geneId)/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present. Missing genes will have their counts set to zero.")))
      
      countsTSub <- as.data.frame(countsT)
      # impute those that are missing and count how many are in the final 209 model genes
      
      emptyDf <- data.frame(matrix(ncol = 2604, nrow = 0))
      colnames(emptyDf) <- geneId
      table(colnames(emptyDf) == geneId)
      # TRUE 
      # 2604
      # Coerce to data.table
      setDT(emptyDf)
      
      # rbind vector (set as a list)
      modelFillNa <- rbind(emptyDf, countsTSub, fill = TRUE)
      
      # Coerce back to a data.frame if you wish
      setDF(modelFillNa)
      rownames(modelFillNa) <- rownames(countsTSub)
      table(colnames(modelFillNa) == geneId)
      # TRUE 
      # 2604  
      
      naRetain <- modelFillNa %>%
        summarise_if(~ any(is.na(.)), ~sum(is.na(.))) %>% 
        t()
      naEns <- modelGeneId[modelGeneId$ENSEMBL %in% rownames(naRetain), ]
      # n = 8 genes with weights in the model
      
      # impute missing genes
      modelFillImp <- modelFillNa %>% 
        replace(is.na(.), 0)
      
      countsT <- modelFillImp
    }
    
    # If all gene are present:
    if(identical(sort(colnames(countsT)), sort(geneId)) == TRUE){
      cat(paste0(((table(colnames(countsT) == geneId)/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present."))
      countsT <- countsT[ , geneId] # order genes
    }
    
    # Transform to log2 counts-per-million
    countsCpm <- edgeR::cpm(countsT, log = TRUE)
    print("Normalising sequencing depth.")
    
    # Z score scale counts using vectors saved from model training
    countsCpmZ <- (countsCpm - meanModelGenesIdentifiHR) / sdModelGenesIdentifiHR
    countsCpmZ <- t(countsCpmZ)
    
    return(countsCpmZ)
  }
  
}

y <- tcgaOvTest %>% 
  column_to_rownames(var = "ENSEMBL") %>% 
  as.matrix()

processCounts(y = y, geneIds = "ENSEMBL")
