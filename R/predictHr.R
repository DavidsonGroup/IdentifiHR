#' Predict HR status using processed gene expresssion counts
#'
#' The "predictHr" function utilises IdentifiHR, a pre-trained penalised logistic regression classifier, to predict homologous recombination status. While IdentifiHR was trained to classify high-grade serous ovarian carcinomas, it can be applied to other cancer types, though the accuracy of classification will be reduced.
#'
#' @param y A numeric matrix of gene expression abundances from the output of the "processCounts()" function, with samples presented in rows and genes presented in columns.
#' @return A data frame containing the sample identifier, the predicted HR status and the probability that a given sample if HRP.
#' \itemize{
#'   \item Sample - A character of sample identifiers taken from the input.
#'   \item hrPrediction - A character vector indicating the discrete homologous recombination (HR) status of a sample, being deficient ("HRD") or proficient ("HRP").
#'   \item predictionProbability - A numeric vector of the probability that a given sample is HRD.
#' }
#' @importFrom dplyr rename
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column
#' @importFrom stats predict
#' @import glmnet 
#' @export
#'
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
#' @examples
#' #See GitHub wiki

predictHr <- function(y) {
  
  bestLambda <- modelIdentifiHR$lambda[99]
  
  print("Predicting HR status using transformed and scaled gene expression counts.")
  # Make predictions using the trained model
  predictionHr <- predict(modelIdentifiHR,
                          newx = t(y),
                          s = bestLambda,
                          type = "class") |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "Sample") |>
    dplyr::rename(hrPrediction = s1)
  
  predictedProb <- predict(modelIdentifiHR,
                           newx = t(y),
                           s = bestLambda,
                           type = "response") |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "Sample") |>
    dplyr::rename(predictionProbability = s1)
  predictedProb$predictionProbability <- (1 -  predictedProb$predictionProbability)
  predictionHrDf <- left_join(predictionHr, predictedProb, by = "Sample")
  
  return(predictionHrDf)
  
}

