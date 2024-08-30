#' Plot the probability that a sample is HRD, as predicted by IdentifiHR
#'
#' The "plotProb()" function allows you to visualise the probability of the input sample/s being "HRD", and to compare this against the probability distribution of the TCGA HGSC testing cohort analysed in our manuscript.
#'
#' @param y A data frame outputted by the "predictHr()" function, including columns "Sample", "hrPrediction" and "predictionProbability".
#' @return A plot of the input cohort and the TCGA test cohort, against the probability of each sample being HRD.
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @export
#'
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
#' @examples
#' data(rawCounts)
#' processedCounts <- processCounts(y = rawCounts, geneIds = "ENSEMBL")
#' predictions <- predictHr(processedCounts)
#' plotProb(predictions)

plotProb <- function(y) {
  
  y$cohort <- "input"
  joinCohort <- rbind(y, tcgaTestPred)
  plot(x = factor(joinCohort$cohort),
       y = joinCohort$predictionProbability, 
       xlab = "Cohort",
       ylab = "Probability of being HRD", 
       col = "white")
  points(x = factor(joinCohort$cohort),
         y = joinCohort$predictionProbability, 
         col = ifelse(joinCohort$hrPrediction == "HRD", "tomato", "darkturquoise"),
         pch = 19,
         cex = 0.8)
  abline(h = 0.5, col = "black", lwd = 1, lty = 2)
  legend("bottomleft", 
         legend = c("HRD", "HRP"),
         col = c( "tomato", "darkturquoise"),
         pch = 20)
  
}

utils::globalVariables(c("tcgaTestPred"))
