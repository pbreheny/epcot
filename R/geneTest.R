#' @export

geneTest <- function(Gene, N.pair, N.trio, vData, maf=.02, cadd=NULL, reqExAC=TRUE, plot=FALSE, floor=ifelse(reqExAC, 0, 1e-4), variants=FALSE) {
  if (missing(vData)) vData <- get("vData", .GlobalEnv)
  # Apply filters
  CADD <- vData$CADD
  if (is.null(cadd)) {
    passCADD <- rep(TRUE, nrow(vData))
  } else {
    passCADD <- CADD >= cadd
  }
  passMAF <- vData$ExAC < maf
  if (reqExAC) {
    passExAC <- vData$ExAC > 0
  } else {
    passExAC <- rep(TRUE, nrow(vData))
  }
  ind <- which(vData$Gene==Gene & passCADD & passMAF & passExAC)  #gives indices of the variants that pass the filters
  if (length(ind) == 0) return(NA)
  MAF <- pmax(floor, vData$ExAC[ind])
  CADD <- CADD[ind]
  NNpair <- N.pair[ind,,drop=FALSE] #variant info for pairs
  if (missing(N.trio)) {
    NNtrio <- matrix(0, nrow(NNpair), 4)
  } else {
    NNtrio <- N.trio[ind,,drop=FALSE] #variant info for trios
  }

  # Calculate Exp
  npair <- apply(NNpair, 1, sum) #sum across the rows. gives # of total pairs
  ntrio <- apply(NNtrio, 1, sum) #sum across the rows. gives # of total trios

  p <- MAF #minor allele freq
  q <- 1-MAF #major allele freq
  probPairs <- 1/4*(4*p*q^3) + 9/16*(4*p^2*q^2) + 2*p^2*q^2 + 4*p^3*q + p^4  #prob that both sisters in a pair have rare variant
  probTrios2 <- 3/8*(4*p*q^3) + 27/64*(4*p^2*q^2)  #prob that 2 of the 3 sisters in a trio have rare variant
  probTrios3 <- 1/8*(4*p*q^3) + 27/64*(4*p^2*q^2) + 2*p^2*q^2 + 4*p^3*q + p^4 #prob that all 3 of the sisters in a trio have rare variant

  Exp <- npair*probPairs + ntrio*(probTrios2+probTrios3)  #expected number of pairs w/ both sisters having variant and trios w/ at least 2 sisters having variant
  Obs <- NNpair[,3]+NNtrio[,3]+NNtrio[,4]

  # plot
  if (plot) {
    upr <- 2*max(sum(Obs), sum(Exp))
    n <- 0:upr
    bp <- barplot(dpois(n, sum(Exp)), border='white', names=n, las=1)
    abline(v=bp[n==sum(Obs),], col="#FF4E37FF", lwd=3)
  }

  # Return
  if (variants) {
    val <- data.frame(Obs=Obs, Exp, p=1-ppois(Obs-0.5, Exp))
  } else {
    val <- c(sum(Obs), sum(Exp), 1-ppois(sum(Obs)-0.5, sum(Exp)))
    names(val) <- c("Obs", "Exp", "p")
  }
  val
}
