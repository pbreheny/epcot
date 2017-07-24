bayesTest <- function(N.pair, N.trio, vData, cutoff=0.2) {
  ind <- vData$ExAC < 0.2
  if (!sum(ind)) return(NA)
  N.pair <- N.pair[ind,,drop=FALSE]
  N.trio <- N.trio[ind,,drop=FALSE]
  vData <- vData[ind]
  npair <- apply(N.pair, 1, sum) #sum across the rows. gives # of total pairs
  ntrio <- apply(N.trio, 1, sum) #sum across the rows. gives # of total trios
  if (missing(vData)) vData <- get("vData", .GlobalEnv)
  CADD <- vData$CADD
  prior <- 10^(-(CADD/10-3))

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

  # Return
  data.frame(Obs, Exp, Prob=1-pgamma(1, Obs + prior, Exp + prior))
}
