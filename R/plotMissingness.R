#' Plot which model genes are missing from input data
#'
#' The "plotMissingness()" function allows you to visualise which weighted genes, if any, required for the IdentifiHR model are missing in your dataset. There are 209 gene required to predict HR status by IdentifiHR, and if any are missing, the model's accuracy will be reduced. 
#'
#' @param y A data frame of all genes required by IdentifiHR, as outputted by the "interrogateMissingness()" function.
#' @return A plot of the beta coefficients of the genes that contribute to prediction in IdentifiHR, coloured by whether they were "missing" or "present" in the input data.
#' @export
#'
#' @author Ashley L Weir, \email{weir.a@@wehi.edu.au}
#' @examples
#' #See GitHub wiki

plotMissingness <- function(y) {
  
  # Order genes by their beta coefficient
  orderY <- y |>
    dplyr::filter(!is.na(betaCoef)) |>
    arrange(betaCoef)
  
  # Count how many genes used for normalisation and weighting were present in the input
  presentNWGene <- length(orderY$inputStatus[orderY$inputStatus == "present"])
  
  # Plot beta coefficients by the ordered 209 genes
  plot(-(orderY$betaCoef),
       xlab = "Genes used by IdentifiHR for HR status predicition",
       ylab = "Beta coeficient",
       main = paste0(presentNWGene, "/209 weighted genes present in input."),
       pch = 19,
       cex = 0.4,
       col = ifelse((orderY$inputStatus == "missing"), "red", "darkgrey"))  
  abline(h = 0, col = "black", lwd = 1, lty = 2)
  legend("bottomleft", c("missing", "present"),
         col = c("red", "darkgrey"), pch=20)
  
}
   

