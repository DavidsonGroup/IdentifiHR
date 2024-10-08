% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelSDGenesIdentifiHR.R
\docType{data}
\name{modelSDGenesIdentifiHR}
\alias{modelSDGenesIdentifiHR}
\title{IdentifiHR: standard deviation of 2604 DE genes}
\format{
A data frame with 2604 rows and 1 column, with ensembl identifiers given as rownames:
\describe{
  \item{sd}{A numeric vector of the standard deviations of the log2 counts-per-million values for genes required by the identifiHR model}
}
}
\source{
Created from training cohort for z score scaling.
}
\usage{
modelSDGenesIdentifiHR
}
\description{
Standard deviation of expression for model genes, to be used in z score scaling.
}
\author{
Ashley L Weir, \email{weir.a@wehi.edu.au}
}
\keyword{datasets}
