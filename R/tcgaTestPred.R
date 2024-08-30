#' IdentifiHR: TCGA HGSC testing cohort HR status predictions
#'
#' HR status predictions in the TCGA HGSC testing cohort.
#'
#' @docType data
#' @format A gene expression counts table with 52125 rows, containing ensembl gene identifiers as rownames, and 10 columns, representing distinct samples.
#'  \describe{
#'   \item{Sample}{A character of sample identifiers taken from the TCGA HGSC tetsing cohort}
#'   \item{hrPrediction}{A character vector indicating the discrete HR status of a sample, being deficient ("HRD") or proficient ("HRP")}
#'   \item{predictionProbability}{A numeric vector of the probability that a given sample is HRD}
#' }
#' 
#' @source Created during the testing and validation of the IdentifiHR model.
#' 
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
"tcgaTestPred"