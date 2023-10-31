#' IdentifiHR
#'
#' The "predictHr" function utilises IdentifiHR, a pre-trained penalised logistic regression classifier, to predict homologous recombination status. While IdentifiHR was trained to classify high-grade serous ovarian carcinomas, it can be applied to other cancer types, though the accuracy of classification will be reduced.
#'
#' @param y A matrix of gene expression counts, with samples presented in rows and genes presented in columns.
#' @param scaled Logical. Has z-score scaling been performed? Default is TRUE, being that the counts matrix has been scaled. If FALSE, z-score scaling of each gene will be performed.
#' @param IdentifiHR The trained HR status predictor
#' @param
#' @return Some output
#' @export
#'
#' @examples
#' ## NOT RUN
#' # data()
#' # predictHr()
#'

predictHr <- function(y,
                      scaled) {

  if(scaled == "TRUE") {

    print("Predicting HR status using scaled gene expression counts")

    # Make predictions using the trained model
    predictionHr <- predict(lmHrSig,
                            newx = y,
                            s = "lambda.min",
                            type = c("class"))

    return(predictionHr)
  }

  if(scaled == "FALSE") {

    print("Z-score scaling gene expression counts and then predicting HR status using scaled gene expression counts")

    # Z-score scale each row of the counts matrix
    scaledCounts <- scale(y,
                          center = TRUE,
                          scale = TRUE)

    # Make predictions using the trained model
    predictionHr <- predict(lmHrSig,
                            newx = scaledCounts,
                            s = "lambda.min",
                            type = c("class"))

    return(predictionHr)
  }

}
predictHr(testDataDe, scaled = TRUE)

