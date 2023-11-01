#' IdentifiHR
#'
#' The "predictHr" function utilises IdentifiHR, a pre-trained penalised logistic regression classifier, to predict homologous recombination status. While IdentifiHR was trained to classify high-grade serous ovarian carcinomas, it can be applied to other cancer types, though the accuracy of classification will be reduced.
#'
#' @param y A numeric matrix of gene expression counts, with samples presented in rows and genes presented in columns.
#' @param logCpm Logical. Have the raw gene expression counts been transformed into log2-counts-per-million? Default is TRUE, being that the counts matrix has been transformed. If FALSE, a log2-counts-per-million transformation will be applied.
#' @param scaled Logical. Has z-score scaling been performed? Default is TRUE, being that the counts matrix has been scaled. If FALSE, z-score scaling of each gene will be performed.
#' @return A data frame containing the sample identifier, the predicted HR status and the probability that this prediction is accurate.
#' @export
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

predictHr <- function(y,
                      logCpm,
                      scaled) {

  if(logCpm == "TRUE") {

    if(scaled == "TRUE") {

      print("Counts have already been logCpm transformed")
      print("Predicting HR status using scaled gene expression counts")

      # Make predictions using the trained model
      predictionHr <- stats::predict(lmHrSig,
                                     newx = y,
                                     s = "lambda.min",
                                     type = "class") %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Sample") %>%
        dplyr::rename(., hrPrediction = lambda.min)

      predictedProb <- stats::predict(lmHrSig,
                                      newx = y,
                                      s = "lambda.min",
                                      type = "response") %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Sample") %>%
        dplyr::rename(., predictionProb = lambda.min)

      predictionHrDf <- left_join(predictionHr, predictedProb, by = "Sample")

      return(predictionHrDf)
    }

  }

  if(logCpm == "TRUE") {

    if(scaled == "FALSE") {

      print("Counts have already been logCpm transformed")
      print("Z-score scaling gene expression counts and then predicting HR status using scaled gene expression counts")

      # Z-score scale each row of the counts matrix
      scaledCounts <- scale(y,
                            center = TRUE,
                            scale = TRUE)

      # Make predictions using the trained model
      predictionHr <- stats::predict(lmHrSig,
                                     newx = y,
                                     s = "lambda.min",
                                     type = "class") %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Sample") %>%
        dplyr::rename(., hrPrediction = lambda.min)

      predictedProb <- stats::predict(lmHrSig,
                                      newx = y,
                                      s = "lambda.min",
                                      type = "response") %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Sample") %>%
        dplyr::rename(., predictionProb = lambda.min)

      predictionHrDf <- left_join(predictionHr, predictedProb, by = "Sample")

      return(predictionHrDf)
    }

  }
  if(logCpm == "FALSE") {

    if(scaled == "TRUE") {

      print("Counts are being logCpm transformed")
      
      y <- cpm(y, log = TRUE)

      print("Predicting HR status using scaled gene expression counts")

      # Make predictions using the trained model
      predictionHr <- stats::predict(lmHrSig,
                                     newx = y,
                                     s = "lambda.min",
                                     type = "class") %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Sample") %>%
        dplyr::rename(., hrPrediction = lambda.min)

      predictedProb <- stats::predict(lmHrSig,
                                      newx = y,
                                      s = "lambda.min",
                                      type = "response") %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Sample") %>%
        dplyr::rename(., predictionProb = lambda.min)

      predictionHrDf <- left_join(predictionHr, predictedProb, by = "Sample")

      return(predictionHrDf)
    }

  }

  if(logCpm == "FALSE") {

    if(scaled == "FALSE") {

      print("Counts are being logCpm transformed")
      
      y <- cpm(y, log = TRUE)
      
      print("Z-score scaling gene expression counts and then predicting HR status using scaled gene expression counts")

      # Z-score scale each row of the counts matrix
      scaledCounts <- scale(y,
                            center = TRUE,
                            scale = TRUE)

      # Make predictions using the trained model
      predictionHr <- stats::predict(lmHrSig,
                                     newx = y,
                                     s = "lambda.min",
                                     type = "class") %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Sample") %>%
        dplyr::rename(., hrPrediction = lambda.min)

      predictedProb <- stats::predict(lmHrSig,
                                      newx = y,
                                      s = "lambda.min",
                                      type = "response") %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Sample") %>%
        dplyr::rename(., predictionProb = lambda.min)

      predictionHrDf <- left_join(predictionHr, predictedProb, by = "Sample")

      return(predictionHrDf)
    }

  }

  if(logCpm == "FALSE") {

    if(scaled == "TRUE") {

      print("Cannot transform data if it has already been scaled. Please use un-scaled counts data as an input.")
    }

  }
}

test <- predictHr(testDataDe, scaled = TRUE)
library(edge)
