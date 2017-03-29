genePerm <- function(Gene, N, vData, cutoff=0.05, plot=FALSE, floor=1e-4) {
  ind <- which(vData$Gene==Gene & vData$ExAC < cutoff)
  if (length(ind) == 0) return(NA)
  maf <- pmax(floor, vData$ExAC[ind])
  NN <- N[ind,,drop=FALSE]
  x <- sum(NN[,3])
  x0 <- geneNull(maf)
  if (plot) {
    Hist(x0, xlim=range(c(x0, x)))
    abline(v=x, col=pal(2)[1], lwd=3)
  }
  mean(x0 >= x)
}
