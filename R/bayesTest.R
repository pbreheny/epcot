#' Bayes test of excess variant frequency with CADD-based prior
#'
#' @param N.pair   Matrix of 0/1/2 counts for sister pairs
#' @param N.trio   Matrix of 0/1/2/3 counts for sister trios
#' @param vData    data.table of variant info (ExAC frequencies and CADD scores)
#' @param cutoff   Exclude variants with ExAC frequencies above this cutoff
#' @param ceiling  Ceiling for prior
#' @param floor    Floor for prior
#' @param rate     Decay rate for prior
#' @param Q        Number of null simulations to run for q value calculation
#'
#' @examples
#' N2 <- epcot_example$N2
#' N3 <- epcot_example$N3
#' vData <- epcot_example$vData
#' res <- bayesTest(N2, N3, vData)
#' head(res[order(res$FDR),])
#'
#' @export

bayesTest <- function(N.pair, N.trio, vData, cutoff=0.2, ceiling=2000, floor=0.01, rate=0.5, Q=100) {

  # Setup
  ind <- vData$ExAC < cutoff
  if (!sum(ind)) return(NA)
  N.pair <- N.pair[ind,,drop=FALSE]
  N.trio <- N.trio[ind,,drop=FALSE]
  vData <- vData[ind]
  npair <- apply(N.pair, 1, sum) #sum across the rows. gives # of total pairs
  ntrio <- apply(N.trio, 1, sum) #sum across the rows. gives # of total trios
  if (missing(vData)) vData <- get("vData", .GlobalEnv)
  CADD <- vData$CADD

  # Determine prior
  tau <- ceiling * exp(-rate * CADD) + floor

  # Calculate Exp
  npair <- apply(N.pair, 1, sum)
  ntrio <- apply(N.trio, 1, sum)

  # Calculate Prob
  p <- vData$ExAC
  q <- 1-p
  probPairs <- 1/4*(4*p*q^3) + 9/16*(4*p^2*q^2) + 2*p^2*q^2 + 4*p^3*q + p^4  #prob that both sisters in a pair have rare variant
  probTrios2 <- 3/8*(4*p*q^3) + 27/64*(4*p^2*q^2)  #prob that 2 of the 3 sisters in a trio have rare variant
  probTrios3 <- 1/8*(4*p*q^3) + 27/64*(4*p^2*q^2) + 2*p^2*q^2 + 4*p^3*q + p^4 #prob that all 3 of the sisters in a trio have rare variant

  # Exp, Obs
  Exp <- npair*probPairs + ntrio*(probTrios2+probTrios3)  #expected number of pairs w/ both sisters having variant and trios w/ at least 2 sisters having variant
  Obs <- N.pair[,3]+N.trio[,3]+N.trio[,4]
  Prob <- 1-pgamma(1, Obs + tau, Exp + tau)

  # qval
  if (Q > 0) {
    nv <- length(Exp)
    NN <- matrix(rpois(nv*Q, Exp), nv, Q)
    PP <- matrix(1-pgamma(1, NN + tau, Exp + tau), nv, Q)
    null <- ecdf(PP)
    FDR <- pmin(nv*(1-null(Prob))/rank(1-Prob), 1)
    return(data.frame(Obs, Exp, Prob, FDR))
  } else {
    return(data.frame(Obs, Exp, Prob))
  }
}
