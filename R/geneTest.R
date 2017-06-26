geneTest <- function(Gene, N, vData, maf=0.05, cadd=NULL, reqExAC=TRUE, plot=FALSE, floor=ifelse(reqExAC, 1e-4, 0), variants=FALSE) {
  # Apply filters
  CADD <- dplyr::coalesce(vData$CADD, 0)
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
  ind <- which(vData$Gene==Gene & passCADD & passMAF & passExAC)
  if (length(ind) == 0) return(NA)
  MAF <- pmax(floor, vData$ExAC[ind])
  CADD <- CADD[ind]
  NN <- N[ind,,drop=FALSE]

  # Calculate Exp
  n <- apply(NN, 1, sum)
  p <- MAF
  q <- 1-MAF
  prob <- 1/4*(4*p*q^3) + 9/16*(4*p^2*q^2) + 2*p^2*q^2 + 4*p^3*q + p^4
  Exp <- n*prob
  x <- sum(NN[,3])

  # plot
  if (plot) {
    n <- 0:ceiling(1.5*x)
    bp <- barplot(dpois(n, sum(Exp)), border='white', names=n, las=1)
    #x0 <- rnorm(10000, sum(Exp), sqrt(sum(Var)))
    #Hist(x0, xlim=range(c(x0, x)))
    abline(v=bp[n==x,], col=pal(2)[1], lwd=3)
  }
  #pnorm((x-sum(Exp))/sqrt(sum(Var)), lower.tail=FALSE)

  # Return
  if (variants) {
    val <- data.frame(Obs=NN[,3], Exp, p=1-ppois(NN[,3]-0.5, Exp))
  } else {
    val <- c(x, sum(Exp), 1-ppois(x-0.5, sum(Exp)))
    names(val) <- c("Obs", "Exp", "p")
  }
  val
}
