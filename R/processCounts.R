#' Process raw gene expression counts for prediction
#'
#' The "processCounts" function prepares gene expression counts for input into the IdentifiHR classifier. 
#'
#' @param y A numeric matrix or dataframe of raw gene expression counts, with genes presented as rownames and samples presented in columns.
#' @param geneIds How are genes annotated? Specify either "ENSEMBL", "HGNC" or "ENTREZ" ("ENSEMBL" is preferred).
#' @return A numeric matrix of counts for model genes, that has been log2-counts-per-million transformed and z-score scaled across genes.
#' @importFrom edgeR cpm
#' @importFrom tibble rownames_to_column 
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr select 
#' @export
#'
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
#' @examples
#' # to add

processCounts <- function(y,
                          geneIds) {
  
  if (!is.matrix(y) & !is.data.frame(y)) {
    
    stop("Input is not a matrix or data frame. The input must be a numeric matrix or data frame, with genes presented in rownames and samples presented in columns.")
    
  }
  
  if (is.matrix(y) & !is.numeric(y)) {
    
    stop("Input is not numeric. The input must be a numeric matrix or data frame, with genes presented in rownames and samples presented in columns.")
    
  }
  
  if (is.data.frame(y) & !is.numeric(as.matrix(y))) {
    
    stop("Input is not numeric. The input must be a numeric matrix or data frame, with genes presented in rownames and samples presented in columns.")
    
  }
  
  if(ncol(y) > 1 & geneIds == "ENSEMBL") {
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(modelSDGenesIdentifiHR)) == TRUE) {
      
      # Extract vector of ensembl identifiers required and associated means and standard deviations
      geneId <- modelGeneId$ENSEMBL
      meanGene <- modelMeanGenesIdentifiHR$mean
      sdGene <- modelSDGenesIdentifiHR$sd
      
    }
    
    # Strip any potential ensembl version numbers
    rownames(y) <- gsub("\\..*","", rownames(y))
    
    # Subset to differentially expressed genes identified in model training (using ensembl IDs)
    y <- y[rownames(y) %in% geneId, ]
    
    message("Subsetting to only model genes.")
    
    # If not all genes are present:
    if (identical(sort(rownames(y)), sort(geneId)) == FALSE) {
      
      message(suppressWarnings(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present."))))
      
      countsTSub <- t(y)
      countsTSub <- as.data.frame(countsTSub)
      countsTSub <- countsTSub[ , intersect(colnames(countsTSub), geneId), drop = FALSE]  # Keep existing columns in order
      
      # Report which genes are missing in warning message
      missingGenes <- setdiff(geneId, colnames(countsTSub))
      warning(paste0("The following genes are missing from the input: ", 
                     paste0(missingGenes, collapse = ", "),
                     ". Missing genes will have their counts set to zero though the model's accuracy will be reduced; we recommend also running the interrogateMissingness() function."))
      
      # Fill missing gene columns with zero values
      
      for (gene in geneId) {
        
        if (!gene %in% colnames(countsTSub)) {
          countsTSub[[gene]] <- 0
          
        }
        
      }
      
      countsTSub <- countsTSub[ , geneId] # Order genes
      completeCounts <- t(countsTSub)
      
    } else if (identical(sort(rownames(y)), sort(geneId)) == TRUE) { # if all genes are present
      
      message(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present.")))
      
      completeCounts <- y[geneId, ] # order genes
      
    }
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(completeCounts)) == TRUE) {
      
      # Transform to log2 counts-per-million
      countsCpm <- edgeR::cpm(completeCounts, log = TRUE) # Rows give genes, columns give samples
      message("Normalising sequencing depth.")
      
      # Z score scale counts using vectors saved from model training
      message("Z-score scaling genes using mean and sd from training cohort.")
      countsCpmZ <- (countsCpm - meanGene) / sdGene
      message("Counts successfully processed for IdentifiHR.")
      return(countsCpmZ)
      
    } else if (identical(rownames(modelMeanGenesIdentifiHR), rownames(completeCounts)) == FALSE) {
      
      stop("Error in gene indexing. Please check input matrix and/or notify package author.")
      
    }
    
  }
  
  if(ncol(y) == 1 & geneIds == "ENSEMBL") {
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(modelSDGenesIdentifiHR)) == TRUE) {
      
      # Extract vector of ensembl identifiers required and associated means and standard deviations
      geneId <- modelGeneId$ENSEMBL
      meanGene <- modelMeanGenesIdentifiHR$mean
      sdGene <- modelSDGenesIdentifiHR$sd
      
    }
    
    # Get sample identifier
    colnameSamp <- colnames(y)
    
    # Add additional ENSEMBL column to avoid R converting input into a vector
    y$ENSEMBL <- gsub("\\..*","", rownames(y))
    
    # Strip any potential ensembl version numbers
    rownames(y) <- gsub("\\..*","", rownames(y))
    
    # Subset to differentially expressed genes identified in model training (using ensembl IDs)
    y <- y[y$ENSEMBL %in% geneId, ]
    
    message("Subsetting to only model genes.")
    
    # If not all genes are present:
    if (identical(sort(rownames(y)), sort(geneId)) == FALSE) {
      
      message(suppressWarnings(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present."))))
      
      # Report which genes are missing in warning message
      missingGenes <- setdiff(geneId, rownames(y))
      warning(paste0("The following genes are missing from the input: ", 
                     paste0(missingGenes, collapse = ", "),
                     ". Missing genes will have their counts set to zero though the model's accuracy will be reduced; we recommend also running the interrogateMissingness() function."))
      
      if (length(missingGenes) > 0) {
        # Create a data frame with missing rows filled with zeros
        missingGeneDf <- data.frame(matrix(0, nrow = length(missingGenes), ncol = ncol(y)))
        rownames(missingGeneDf) <- missingGenes
        colnames(missingGeneDf) <- colnames(y)
        
        # Combine the existing data frame with the missing rows
        countsTSub <- rbind(y, missingGeneDf)
      }
      
      completeCounts <- countsTSub[geneId, , drop = FALSE] # Order genes
      rownameGeneIdx <- rownames(completeCounts)
      
    } else if (identical(sort(rownames(y)), sort(geneId)) == TRUE) { # if all genes are present
      
      message(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present.")))
      
      completeCounts <- y[geneId, ] # order genes
      rownameGeneIdx <- rownames(completeCounts)
      
    }
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownameGeneIdx) == TRUE) {
      
      completeCounts <- completeCounts[ , colnames(completeCounts) != "ENSEMBL"]
      
      # Transform to log2 counts-per-million
      countsCpm <- edgeR::cpm(completeCounts, log = TRUE) # Rows give genes, columns give samples
      message("Normalising sequencing depth.")
      
      # Z score scale counts using vectors saved from model training
      message("Z-score scaling genes using mean and sd from training cohort.")
      countsCpmZ <- (countsCpm - meanGene) / sdGene
      
      # Make into a data frame and re-assign rownames and colname
      countsCpmZ <- as.data.frame(countsCpmZ)
      rownames(countsCpmZ) <- rownameGeneIdx
      colnames(countsCpmZ) <- colnameSamp
      message("Counts successfully processed for IdentifiHR.")
      return(countsCpmZ)
      
    } else if (identical(rownames(modelMeanGenesIdentifiHR), rownameGeneIdx) == FALSE) {
      
      stop("Error in gene indexing. Please check input matrix and/or notify package author.")
      
    }
    
  }
  
  if(ncol(y) > 1 & geneIds == "HGNC") {
    
    # Replace hgnc with ensembl
    y <- y |>
      as.data.frame() |>
      rownames_to_column(var = "hgnc_symbol") |>
      left_join(geneId, by = "hgnc_symbol") |>
      dplyr::filter(!is.na(ENSEMBL)) |>
      column_to_rownames(var = "ENSEMBL") |>
      dplyr::select(-c(hgnc_symbol, entrezgene_id))
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(modelSDGenesIdentifiHR)) == TRUE) {
      
      # Extract vector of ensembl identifiers required and associated means and standard deviations
      geneId <- modelGeneId$ENSEMBL
      meanGene <- modelMeanGenesIdentifiHR$mean
      sdGene <- modelSDGenesIdentifiHR$sd
      
    }
    
    # Strip any potential ensembl version numbers
    rownames(y) <- gsub("\\..*","", rownames(y))
    
    # Subset to differentially expressed genes identified in model training (using ensembl IDs)
    y <- y[rownames(y) %in% geneId, ]
    
    message("Subsetting to only model genes.")
    
    # If not all genes are present:
    if (identical(sort(rownames(y)), sort(geneId)) == FALSE) {
      
      message(suppressWarnings(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present."))))
      
      countsTSub <- t(y)
      countsTSub <- as.data.frame(countsTSub)
      countsTSub <- countsTSub[ , intersect(colnames(countsTSub), geneId), drop = FALSE]  # Keep existing columns in order
      
      # Report which genes are missing in warning message
      missingGenes <- setdiff(geneId, colnames(countsTSub))
      warning(paste0("The following genes are missing from the input: ", 
                     paste0(missingGenes, collapse = ", "),
                     ". Missing genes will have their counts set to zero though the model's accuracy will be reduced; we recommend also running the interrogateMissingness() function."))
      
      # Fill missing gene columns with zero values
      
      for (gene in geneId) {
        
        if (!gene %in% colnames(countsTSub)) {
          countsTSub[[gene]] <- 0
          
        }
        
      }
      
      countsTSub <- countsTSub[ , geneId] # Order genes
      completeCounts <- t(countsTSub)
      
    } else if (identical(sort(rownames(y)), sort(geneId)) == TRUE) { # if all genes are present
      
      message(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present.")))
      
      completeCounts <- y[geneId, ] # order genes
      
    }
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(completeCounts)) == TRUE) {
      
      # Transform to log2 counts-per-million
      countsCpm <- edgeR::cpm(completeCounts, log = TRUE) # Rows give genes, columns give samples
      message("Normalising sequencing depth.")
      
      # Z score scale counts using vectors saved from model training
      message("Z-score scaling genes using mean and sd from training cohort.")
      countsCpmZ <- (countsCpm - meanGene) / sdGene
      message("Counts successfully processed for IdentifiHR.")
      return(countsCpmZ)
      
    } else if (identical(rownames(modelMeanGenesIdentifiHR), rownames(completeCounts)) == FALSE) {
      
      stop("Error in gene indexing. Please check input matrix and/or notify package author.")
      
    }
    
  }
  
  if(ncol(y) == 1 & geneIds == "HGNC") {
    
    # Replace hgnc with ensembl
    y <- y |>
      as.data.frame() |>
      rownames_to_column(var = "hgnc_symbol") |>
      left_join(geneId, by = "hgnc_symbol") |>
      dplyr::filter(!is.na(ENSEMBL)) |>
      column_to_rownames(var = "ENSEMBL") |>
      dplyr::select(-c(hgnc_symbol, entrezgene_id))
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(modelSDGenesIdentifiHR)) == TRUE) {
      
      # Extract vector of ensembl identifiers required and associated means and standard deviations
      geneId <- modelGeneId$ENSEMBL
      meanGene <- modelMeanGenesIdentifiHR$mean
      sdGene <- modelSDGenesIdentifiHR$sd
      
    }
    
    # Get sample identifier
    colnameSamp <- colnames(y)
    
    # Add additional ENSEMBL column to avoid R converting input into a vector
    y$ENSEMBL <- gsub("\\..*","", rownames(y))
    
    # Strip any potential ensembl version numbers
    rownames(y) <- gsub("\\..*","", rownames(y))
    
    # Subset to differentially expressed genes identified in model training (using ensembl IDs)
    y <- y[y$ENSEMBL %in% geneId, ]
    
    message("Subsetting to only model genes.")
    
    # If not all genes are present:
    if (identical(sort(rownames(y)), sort(geneId)) == FALSE) {
      
      message(suppressWarnings(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present."))))
      
      # Report which genes are missing in warning message
      missingGenes <- setdiff(geneId, rownames(y))
      warning(paste0("The following genes are missing from the input: ", 
                     paste0(missingGenes, collapse = ", "),
                     ". Missing genes will have their counts set to zero though the model's accuracy will be reduced; we recommend also running the interrogateMissingness() function."))
      
      if (length(missingGenes) > 0) {
        # Create a data frame with missing rows filled with zeros
        missingGeneDf <- data.frame(matrix(0, nrow = length(missingGenes), ncol = ncol(y)))
        rownames(missingGeneDf) <- missingGenes
        colnames(missingGeneDf) <- colnames(y)
        
        # Combine the existing data frame with the missing rows
        countsTSub <- rbind(y, missingGeneDf)
      }
      
      completeCounts <- countsTSub[geneId, , drop = FALSE] # Order genes
      rownameGeneIdx <- rownames(completeCounts)
      
    } else if (identical(sort(rownames(y)), sort(geneId)) == TRUE) { # if all genes are present
      
      message(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present.")))
      
      completeCounts <- y[geneId, ] # order genes
      rownameGeneIdx <- rownames(completeCounts)
      
    }
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownameGeneIdx) == TRUE) {
      
      completeCounts <- completeCounts[ , colnames(completeCounts) != "ENSEMBL"]
      
      # Transform to log2 counts-per-million
      countsCpm <- edgeR::cpm(completeCounts, log = TRUE) # Rows give genes, columns give samples
      message("Normalising sequencing depth.")
      
      # Z score scale counts using vectors saved from model training
      message("Z-score scaling genes using mean and sd from training cohort.")
      countsCpmZ <- (countsCpm - meanGene) / sdGene
      
      # Make into a data frame and re-assign rownames and colname
      countsCpmZ <- as.data.frame(countsCpmZ)
      rownames(countsCpmZ) <- rownameGeneIdx
      colnames(countsCpmZ) <- colnameSamp
      message("Counts successfully processed for IdentifiHR.")
      return(countsCpmZ)
      
    } else if (identical(rownames(modelMeanGenesIdentifiHR), rownameGeneIdx) == FALSE) {
      
      stop("Error in gene indexing. Please check input matrix and/or notify package author.")
      
    }
    
  }
  
  if(ncol(y) > 1 & geneIds == "ENTREZ") {
    
    # Replace entrez with ensembl
    y <- y |>
      as.data.frame() |>
      rownames_to_column(var = "entrezgene_id") |>
      left_join(geneId, by = "entrezgene_id") |>
      dplyr::filter(!is.na(ENSEMBL)) |>
      column_to_rownames(var = "ENSEMBL") |>
      dplyr::select(-c(hgnc_symbol, entrezgene_id))
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(modelSDGenesIdentifiHR)) == TRUE) {
      
      # Extract vector of ensembl identifiers required and associated means and standard deviations
      geneId <- modelGeneId$ENSEMBL
      meanGene <- modelMeanGenesIdentifiHR$mean
      sdGene <- modelSDGenesIdentifiHR$sd
      
    }
    
    # Strip any potential ensembl version numbers
    rownames(y) <- gsub("\\..*","", rownames(y))
    
    # Subset to differentially expressed genes identified in model training (using ensembl IDs)
    y <- y[rownames(y) %in% geneId, ]
    
    message("Subsetting to only model genes.")
    
    # If not all genes are present:
    if (identical(sort(rownames(y)), sort(geneId)) == FALSE) {
      
      message(suppressWarnings(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present."))))
      
      countsTSub <- t(y)
      countsTSub <- as.data.frame(countsTSub)
      countsTSub <- countsTSub[ , intersect(colnames(countsTSub), geneId), drop = FALSE]  # Keep existing columns in order
      
      # Report which genes are missing in warning message
      missingGenes <- setdiff(geneId, colnames(countsTSub))
      warning(paste0("The following genes are missing from the input: ", 
                     paste0(missingGenes, collapse = ", "),
                     ". Missing genes will have their counts set to zero though the model's accuracy will be reduced; we recommend also running the interrogateMissingness() function."))
      
      # Fill missing gene columns with zero values
      
      for (gene in geneId) {
        
        if (!gene %in% colnames(countsTSub)) {
          countsTSub[[gene]] <- 0
          
        }
        
      }
      
      countsTSub <- countsTSub[ , geneId] # Order genes
      completeCounts <- t(countsTSub)
      
    } else if (identical(sort(rownames(y)), sort(geneId)) == TRUE) { # if all genes are present
      
      message(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present.")))
      
      completeCounts <- y[geneId, ] # order genes
      
    }
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(completeCounts)) == TRUE) {
      
      # Transform to log2 counts-per-million
      countsCpm <- edgeR::cpm(completeCounts, log = TRUE) # Rows give genes, columns give samples
      message("Normalising sequencing depth.")
      
      # Z score scale counts using vectors saved from model training
      message("Z-score scaling genes using mean and sd from training cohort.")
      countsCpmZ <- (countsCpm - meanGene) / sdGene
      message("Counts successfully processed for IdentifiHR.")
      return(countsCpmZ)
      
    } else if (identical(rownames(modelMeanGenesIdentifiHR), rownames(completeCounts)) == FALSE) {
      
      stop("Error in gene indexing. Please check input matrix and/or notify package author.")
      
    }
    
  }
  
  if(ncol(y) == 1 & geneIds == "HGNC") {
    
    # Replace entrez with ensembl
    y <- y |>
      as.data.frame() |>
      rownames_to_column(var = "entrezgene_id") |>
      left_join(geneId, by = "entrezgene_id") |>
      dplyr::filter(!is.na(ENSEMBL)) |>
      column_to_rownames(var = "ENSEMBL") |>
      dplyr::select(-c(hgnc_symbol, entrezgene_id))
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownames(modelSDGenesIdentifiHR)) == TRUE) {
      
      # Extract vector of ensembl identifiers required and associated means and standard deviations
      geneId <- modelGeneId$ENSEMBL
      meanGene <- modelMeanGenesIdentifiHR$mean
      sdGene <- modelSDGenesIdentifiHR$sd
      
    }
    
    # Get sample identifier
    colnameSamp <- colnames(y)
    
    # Add additional ENSEMBL column to avoid R converting input into a vector
    y$ENSEMBL <- gsub("\\..*","", rownames(y))
    
    # Strip any potential ensembl version numbers
    rownames(y) <- gsub("\\..*","", rownames(y))
    
    # Subset to differentially expressed genes identified in model training (using ensembl IDs)
    y <- y[y$ENSEMBL %in% geneId, ]
    
    message("Subsetting to only model genes.")
    
    # If not all genes are present:
    if (identical(sort(rownames(y)), sort(geneId)) == FALSE) {
      
      message(suppressWarnings(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present."))))
      
      # Report which genes are missing in warning message
      missingGenes <- setdiff(geneId, rownames(y))
      warning(paste0("The following genes are missing from the input: ", 
                     paste0(missingGenes, collapse = ", "),
                     ". Missing genes will have their counts set to zero though the model's accuracy will be reduced; we recommend also running the interrogateMissingness() function."))
      
      if (length(missingGenes) > 0) {
        # Create a data frame with missing rows filled with zeros
        missingGeneDf <- data.frame(matrix(0, nrow = length(missingGenes), ncol = ncol(y)))
        rownames(missingGeneDf) <- missingGenes
        colnames(missingGeneDf) <- colnames(y)
        
        # Combine the existing data frame with the missing rows
        countsTSub <- rbind(y, missingGeneDf)
      }
      
      completeCounts <- countsTSub[geneId, , drop = FALSE] # Order genes
      rownameGeneIdx <- rownames(completeCounts)
      
    } else if (identical(sort(rownames(y)), sort(geneId)) == TRUE) { # if all genes are present
      
      message(cat(paste0(((length(intersect(rownames(y), geneId))/2604)[-2])*100, "%", sep = " ", "of the 2604 genes required for IdentifiHR are present.")))
      
      completeCounts <- y[geneId, ] # order genes
      rownameGeneIdx <- rownames(completeCounts)
      
    }
    
    if (identical(rownames(modelMeanGenesIdentifiHR), rownameGeneIdx) == TRUE) {
      
      completeCounts <- completeCounts[ , colnames(completeCounts) != "ENSEMBL"]
      
      # Transform to log2 counts-per-million
      countsCpm <- edgeR::cpm(completeCounts, log = TRUE) # Rows give genes, columns give samples
      message("Normalising sequencing depth.")
      
      # Z score scale counts using vectors saved from model training
      message("Z-score scaling genes using mean and sd from training cohort.")
      countsCpmZ <- (countsCpm - meanGene) / sdGene
      
      # Make into a data frame and re-assign rownames and colname
      countsCpmZ <- as.data.frame(countsCpmZ)
      rownames(countsCpmZ) <- rownameGeneIdx
      colnames(countsCpmZ) <- colnameSamp
      message("Counts successfully processed for IdentifiHR.")
      return(countsCpmZ)
      
    } else if (identical(rownames(modelMeanGenesIdentifiHR), rownameGeneIdx) == FALSE) {
      
      stop("Error in gene indexing. Please check input matrix and/or notify package author.")
      
    }
    
  }
  
}

