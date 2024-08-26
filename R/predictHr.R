#' IdentifiHR
#'
#' The "predictHr" function utilises IdentifiHR, a pre-trained penalised logistic regression classifier, to predict homologous recombination status. While IdentifiHR was trained to classify high-grade serous ovarian carcinomas, it can be applied to other cancer types, though the accuracy of classification will be reduced.
#'
#' @param y A numeric matrix of gene expression counts, with samples presented in rows and genes presented in columns.
#' @return A data frame containing the sample identifier, the predicted HR status and the probability that a given sample if HRP.
#' @import dplyr, stats
#' @export
#'
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
#'
#' @examples
#' ## NOT RUN
#' # Load packages
#' # library(dplyr)
#' # library(tidyverse)
#' # library(stats)
#' # library(edgeR)
#' # Load data
#' # data(tcgaOvScaled)
#' # Predict HR status with IdentifiHR
#' # predictHr(tcgaOvScaled, scaled = TRUE)

predictHr <- function(y) {
  
  print("Predicting HR status using scaled gene expression counts")
  # Make predictions using the trained model
  predictionHr <- stats::predict(modelIdentifiHR,
                                 newx = y,
                                 s = "lambda.min",
                                 type = "class") |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "Sample") |>
    dplyr::rename(., hrPrediction = lambda.min)
  predictedProb <- stats::predict(modelIdentifiHR,
                                  newx = y,
                                  s = "lambda.min",
                                  type = "response") |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "Sample") |>
    dplyr::rename(., predictionProb = lambda.min)
  predictionHrDf <- left_join(predictionHr, predictedProb, by = "Sample")
  return(predictionHrDf)
  
}

predictHr(countsCpmZ)
