#' IdentifiHR
#'
#' The "processCounts" function prepares gene expression counts for input into the IdentifiHR classifier. 
#'
#' @param y A numeric matrix of raw gene expression counts, with genes presented in rownames and samples presented in columns.
#' @param geneIds How are genes annotated? Specify either "ENSEMBL", "HGNC" or "ENTREZ" ("ENSEMBL" is preferred).
#' @return A numeric matrix of counts for model genes, that has been log2-counts-per-million transformed and z-score scaled across genes.
#' @import dplyr
#' @importFrom data.table setDT
#' @importFrom data.table setDF
#' @importFrom edgeR cpm
#' @export
#'
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
#' @examples
#' # to add

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
    
    message("Subsetting to only model genes.")
    
    # Check all genes required are present
    countsT <- t(y)
    
    # If not all genes are present:
    if (identical(sort(colnames(countsT)), sort(geneId)) == FALSE) {
      
      message(suppressWarnings(cat(paste0(((table(colnames(countsT) == geneId)/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present. Missing genes will have their counts set to zero."))))
      
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
      
      # impute missing genes
      modelFillImp <- modelFillNa %>% 
        replace(is.na(.), 0)
      
      countsT <- t(modelFillImp[ , geneId]) # order genes
      
    } else if (identical(sort(colnames(countsT)), sort(geneId)) == TRUE) {
      
      countsT <- countsT[ , geneId] # order genes
      message(cat(paste0(((table(colnames(countsT) == geneId)/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present.")))
      counts <- t(countsT)
      
    }
    
    # Transform to log2 counts-per-million
    countsCpm <- edgeR::cpm(countsT, log = TRUE)
    message("Normalising sequencing depth.")
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(modelSDGenesIdentifiHR)) == TRUE) {
      
      if (identical(rownames(modelMeanGenesIdentifiHR), colnames(countsCpm)) == TRUE) {
        
        message("Z-score scaling genes using mean and sd from training cohort.")
        # Z score scale counts using vectors saved from model training
        countsCpmZ <- sweep(countsCpm, 2, modelMeanGenesIdentifiHR$mean, FUN = "-")
        countsCpmZ <- sweep(countsCpmZ, 2, modelSDGenesIdentifiHR$sd, FUN = "/")
        countsCpmZ <- t(countsCpmZ) 
        message("Counts successfully processed for IdentifiHR.")
        return(countsCpmZ)
        
      }
      
    } else if (identical(rownames(modelMeanGenesIdentifiHR), colnames(countsCpm)) == FALSE) {
        
        stop("Error in gene indexing. Please check input matrix and/or notify package author.")
        
      }
      
    }
    
  if(geneIds == "HGNC") {
    
    # Replace hgnc with ensembl
    y <- y |>
      as.data.frame() |>
      rownames_to_column(var = "hgnc_symbol") |>
      left_join(geneId, by = "hgnc_symbol") |>
      dplyr::filter(!is.na(ENSEMBL)) |>
      column_to_rownames(var = "ENSEMBL") |>
      dplyr::select(-c(hgnc_symbol, entrezgene_id))
    
    # Extract vector of ensembl identifiers required
    geneId <- modelGeneId$ENSEMBL
    
    # Strip any potential ensembl version numbers
    rownames(y) <- gsub("\\..*","", rownames(y))
    
    # Subset to differentially expressed genes identified in model training (using ensembl IDs)
    y <- y[rownames(y) %in% geneId, ]
    
    message("Subsetting to only model genes.")
    
    # Check all genes required are present
    countsT <- t(y)
    
    # If not all genes are present:
    if (identical(sort(colnames(countsT)), sort(geneId)) == FALSE) {
      
      message(suppressWarnings(cat(paste0(((table(colnames(countsT) == geneId)/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present. Missing genes will have their counts set to zero."))))
      
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
      
      # impute missing genes
      modelFillImp <- modelFillNa %>% 
        replace(is.na(.), 0)
      
      countsT <- modelFillImp
      countsT <- countsT[ , geneId] # order genes
      
    } else if (identical(sort(colnames(countsT)), sort(geneId)) == TRUE) {
      
      countsT <- countsT[ , geneId] # order genes
      message(cat(paste0(((table(colnames(countsT) == geneId)/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present.")))
      
    }
    
    # Transform to log2 counts-per-million
    countsCpm <- edgeR::cpm(countsT, log = TRUE)
    message("Normalising sequencing depth.")
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(modelSDGenesIdentifiHR)) == TRUE) {
      
      if (identical(rownames(modelMeanGenesIdentifiHR), colnames(countsCpm)) == TRUE) {
        
        message("Z-score scaling genes using mean and sd from training cohort.")
        # Z score scale counts using vectors saved from model training
        countsCpmZ <- sweep(countsCpm, 2, modelMeanGenesIdentifiHR$mean, FUN = "-")
        countsCpmZ <- sweep(countsCpmZ, 2, modelSDGenesIdentifiHR$sd, FUN = "/")
        countsCpmZ <- t(countsCpmZ) 
        message("Counts successfully processed for IdentifiHR.")
        return(countsCpmZ)
        
      }
      
    }
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(modelSDGenesIdentifiHR)) == TRUE) {
      
      if (identical(rownames(modelMeanGenesIdentifiHR), colnames(countsCpm)) == FALSE) {
        
        stop("Error in gene indexing. Please check input matrix and/or notify package author.")
        
      }
      
    }
    
  }
  
  if(geneIds == "ENTREZ") {
    
    # Replace entrez id with ensembl
    y <- y |>
      as.data.frame() |>
      rownames_to_column(var = "entrezgene_id") |>
      left_join(geneId, by = "entrezgene_id") |>
      dplyr::filter(!is.na(ENSEMBL)) |>
      column_to_rownames(var = "ENSEMBL") |>
      dplyr::select(-c(hgnc_symbol, entrezgene_id))
    
    # Extract vector of ensembl identifiers required
    geneId <- modelGeneId$ENSEMBL
    
    # Strip any potential ensembl version numbers
    rownames(y) <- gsub("\\..*","", rownames(y))
    
    # Subset to differentially expressed genes identified in model training (using ensembl IDs)
    y <- y[rownames(y) %in% geneId, ]
    
    message("Subsetting to only model genes.")
    
    # Check all genes required are present
    countsT <- t(y)
    
    # If not all genes are present:
    if (identical(sort(colnames(countsT)), sort(geneId)) == FALSE) {
      
      message(suppressWarnings(cat(paste0(((table(colnames(countsT) == geneId)/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present. Missing genes will have their counts set to zero."))))
      
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
      
      # impute missing genes
      modelFillImp <- modelFillNa %>% 
        replace(is.na(.), 0)
      
      countsT <- modelFillImp
      countsT <- countsT[ , geneId] # order genes
      
    } else if (identical(sort(colnames(countsT)), sort(geneId)) == TRUE) {
      
      countsT <- countsT[ , geneId] # order genes
      message(cat(paste0(((table(colnames(countsT) == geneId)/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present.")))
      
    }
    
    # Transform to log2 counts-per-million
    countsCpm <- edgeR::cpm(countsT, log = TRUE)
    message("Normalising sequencing depth.")
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(modelSDGenesIdentifiHR)) == TRUE) {
      
      if (identical(rownames(modelMeanGenesIdentifiHR), colnames(countsCpm)) == TRUE) {
        
        message("Z-score scaling genes using mean and sd from training cohort.")
        # Z score scale counts using vectors saved from model training
        countsCpmZ <- sweep(countsCpm, 2, modelMeanGenesIdentifiHR$mean, FUN = "-")
        countsCpmZ <- sweep(countsCpmZ, 2, modelSDGenesIdentifiHR$sd, FUN = "/")
        countsCpmZ <- t(countsCpmZ) 
        message("Counts successfully processed for IdentifiHR.")
        return(countsCpmZ)
        
      }
      
    }
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(modelSDGenesIdentifiHR)) == TRUE) {
      
      if (identical(rownames(modelMeanGenesIdentifiHR), colnames(countsCpm)) == FALSE) {
        
        stop("Error in gene indexing. Please check input matrix and/or notify package author.")
        
      }
      
    }
    
  }
  
}
