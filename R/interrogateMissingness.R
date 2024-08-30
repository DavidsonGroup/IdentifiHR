#' Interrogate which model genes are missing from input data
#'
#' The "interrogateMissingness()" function allows you to examine which genes, if any, required for the IdentifiHR model are missing in your dataset. There are 2604 genes required for normalising library size differences and 209 gene required to predict HR status by IdentifiHR. 
#'
#' @param y A numeric matrix or data frame of raw gene expression counts, with genes presented as rownames and samples presented in columns.
#' @param geneIds How are genes annotated? Specify either "ENSEMBL", "HGNC" or "ENTREZ" ("ENSEMBL" is preferred).
#' @return A data frame containing the following columns:
#' \itemize{
#'   \item ensembl_id - A character vector of ensembl identifiers for genes needed in the model.
#'   \item hgnc_symbol - A character vector of hgnc symbols for genes needed in the model.
#'   \item entrezgene_id - A numeric vector of entrez identifiers for genes needed in the model.
#'   \item betaCoef - A numeric vector of beta coefficients for genes that are weighted towards IdentifiHR predictions. NAs indicate that a gene is only used in normalisation.
#'   \item inputStatus - A character vector indicating if a gene was "present"  or "missing in the input data. 
#'   \item normWeightedGene - A character vector indicating if a gene contributes to only normalisation, being "normGene", or contributes to both normalisation and model predictions, being "normWeightedGene".
#' }
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom dplyr case_when
#' @importFrom dplyr filter
#' @export
#'
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
#' @examples
#' data(rawCounts)
#' missingGenes <- interrogateMissingness(y = rawCounts, geneIds = "ENSEMBL")
#' head(missingGenes)

interrogateMissingness <- function(y,
                                   geneIds) {
  
  if (geneIds == "ENSEMBL") {
    
    # Make input a data frame
    y <- as.data.frame(y)
    
    # Create a column of ensembl identifiers for joining
    y$ensembl_id <- rownames(y)
    
    # Make a column to note which genes were present in both the input and the model
    y$inputStatus <- "present"
    
    # Join with modelGeneId, and select only columns that are needed to interrogate missingness
    joinGenes <- dplyr::left_join(modelGeneId, y, by = "ensembl_id") |>
      dplyr::select(ensembl_id, hgnc_symbol, entrezgene_id, betaCoef, inputStatus)
    joinGenes$normWeightedGene <- dplyr::case_when(!is.na(joinGenes$betaCoef) ~ "normWeightedGene",
                                            is.na(joinGenes$betaCoef) ~ "normGene")
    joinGenes$inputStatus <- dplyr::case_when(!is.na(joinGenes$inputStatus) ~ "present",
                                              is.na(joinGenes$inputStatus) ~ "missing")
    
    if(all(joinGenes$inputStatus == "present")) {
      
      stop("None of the required genes are missing from the input.")
      
    }
    
    missingNormGenes <- joinGenes |> 
      dplyr::filter(normWeightedGene == "normGene") |>
      dplyr::filter(is.na(inputStatus))
    message(paste0("The following genes used in normalising library size are missing from the input: ", 
                   paste0(missingNormGenes$ensembl_id, collapse = ", "),
                   ". Their counts will have been set to zero prior to processing and predicting HR status. The  model's accuracy may have been slightly reduced; we recommend also running the plotMissingness() function."))
    
    missingNormWeightedGenes <- joinGenes |> 
      dplyr::filter(normWeightedGene == "normWeightedGene") |>
      dplyr::filter(is.na(inputStatus))  
    message(paste0("The following genes used in normalising library size, that are also weighted for IdentifiHR predictions, are missing from the input: ", 
                   paste0(missingNormWeightedGenes$ensembl_id, collapse = ", "),
                   ". Their counts will have been set to zero prior to processing and predicting HR status. The  model's accuracy will have been reduced; we recommend also running the plotMissingness() function."))
    
    return(joinGenes)
    
  }

  if (geneIds == "HGNC") {
    
    # Make input a data frame
    y <- as.data.frame(y)
    
    # Create a column of hgnc identifiers for joining
    y$hgnc_symbol <- rownames(y)
    y$hgnc_symbol <- as.character(y$hgnc_symbol)
    
    # Make a column to note which genes were present in both the input and the model
    y$inputStatus <- "present"
    
    # Join with modelGeneId, and select only columns that are needed to interrogate missingness
    joinGenes <- dplyr::left_join(modelGeneId, y, by = "hgnc_symbol") |>
      dplyr::filter(!is.na(ensembl_id)) |>
      dplyr::select(ensembl_id, hgnc_symbol, entrezgene_id, betaCoef, inputStatus)
    joinGenes$normWeightedGene <- dplyr::case_when(!is.na(joinGenes$betaCoef) ~ "normWeightedGene",
                                                   is.na(joinGenes$betaCoef) ~ "normGene")
    joinGenes$inputStatus <- dplyr::case_when(!is.na(joinGenes$inputStatus) ~ "present",
                                              is.na(joinGenes$inputStatus) ~ "missing")
    
    if(all(joinGenes$inputStatus == "present")) {
      
      stop("None of the required genes are missing from the input.")
      
    }
    
    missingNormGenes <- joinGenes |> 
      dplyr::filter(normWeightedGene == "normGene") |>
      dplyr::filter(is.na(inputStatus))
    message(paste0("The following genes used in normalising library size are missing from the input: ", 
                   paste0(missingNormGenes$ensembl_id, collapse = ", "),
                   ". Their counts will have been set to zero prior to processing and predicting HR status. The  model's accuracy may have been slightly reduced; we recommend also running the plotMissingness() function."))
    
    missingNormWeightedGenes <- joinGenes |> 
      dplyr::filter(normWeightedGene == "normWeightedGene") |>
      dplyr::filter(is.na(inputStatus))  
    message(paste0("The following genes used in normalising library size, that are also weighted for IdentifiHR predictions, are missing from the input: ", 
                   paste0(missingNormWeightedGenes$ensembl_id, collapse = ", "),
                   ". Their counts will have been set to zero prior to processing and predicting HR status. The  model's accuracy will have been reduced; we recommend also running the plotMissingness() function."))
    
    return(joinGenes)
    
  }
  
  if (geneIds == "ENTREZ") {
    
    # Make input a data frame
    y <- as.data.frame(y)
    
    # Create a column of hgnc identifiers for joining
    y$entrezgene_id <- rownames(y)
    y$entrezgene_id <- as.numeric(y$entrezgene_id)
    
    # Make a column to note which genes were present in both the input and the model
    y$inputStatus <- "present"
    
    # Join with modelGeneId, and select only columns that are needed to interrogate missingness
    joinGenes <- dplyr::left_join(modelGeneId, y, by = "entrezgene_id") |>
      dplyr::filter(!is.na(ensembl_id)) |>
      dplyr::select(ensembl_id, hgnc_symbol, entrezgene_id, betaCoef, inputStatus)
    joinGenes$normWeightedGene <- dplyr::case_when(!is.na(joinGenes$betaCoef) ~ "normWeightedGene",
                                                   is.na(joinGenes$betaCoef) ~ "normGene")
    joinGenes$inputStatus <- dplyr::case_when(!is.na(joinGenes$inputStatus) ~ "present",
                                              is.na(joinGenes$inputStatus) ~ "missing")
    
    if(all(joinGenes$inputStatus == "present")) {
      
      stop("None of the required genes are missing from the input.")
      
    }
    
    missingNormGenes <- joinGenes |> 
      dplyr::filter(normWeightedGene == "normGene") |>
      dplyr::filter(is.na(inputStatus))
    message(paste0("The following genes used in normalising library size are missing from the input: ", 
                   paste0(missingNormGenes$ensembl_id, collapse = ", "),
                   ". Their counts will have been set to zero prior to processing and predicting HR status. The  model's accuracy may have been slightly reduced; we recommend also running the plotMissingness() function."))
    
    missingNormWeightedGenes <- joinGenes |> 
      dplyr::filter(normWeightedGene == "normWeightedGene") |>
      dplyr::filter(is.na(inputStatus))  
    message(paste0("The following genes used in normalising library size, that are also weighted for IdentifiHR predictions, are missing from the input: ", 
                   paste0(missingNormWeightedGenes$ensembl_id, collapse = ", "),
                   ". Their counts will have been set to zero prior to processing and predicting HR status. The  model's accuracy will have been reduced; we recommend also running the plotMissingness() function."))
    
    return(joinGenes)
    
  }
  
}
