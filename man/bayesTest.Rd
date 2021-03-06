% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesTest.R
\name{bayesTest}
\alias{bayesTest}
\title{Bayes test of excess variant frequency with CADD-based prior}
\usage{
bayesTest(
  N.pair,
  N.trio,
  vData,
  cutoff = 0.2,
  ceiling = 2000,
  floor = 0.01,
  rate = 0.5,
  Q = 100
)
}
\arguments{
\item{N.pair}{Matrix of 0/1/2 counts for sister pairs}

\item{N.trio}{Matrix of 0/1/2/3 counts for sister trios}

\item{vData}{data.table of variant info (ExAC frequencies and CADD scores)}

\item{cutoff}{Exclude variants with ExAC frequencies above this cutoff}

\item{ceiling}{Ceiling for prior}

\item{floor}{Floor for prior}

\item{rate}{Decay rate for prior}

\item{Q}{Number of null simulations to run for q value calculation}
}
\description{
Bayes test of excess variant frequency with CADD-based prior
}
\examples{
N2 <- epcot_example$N2
N3 <- epcot_example$N3
vData <- epcot_example$vData
res <- bayesTest(N2, N3, vData)
head(res[order(res$FDR),])

}
