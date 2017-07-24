varSummary <- function(v, N.pair, N.trio, vData) {
  if (missing(vData)) vData <- get("vData", .GlobalEnv)
  ind <- match(v, rownames(N.pair))
  NP <- N.pair[ind,,drop=FALSE]
  colnames(NP) <- gsub("N", "P", colnames(NP))
  V <- vData[ind]
  out <- cbind(V, RVIS=rvis(V$Gene), NP)
  if (!missing(N.trio)) {
    NT <- N.trio[ind,,drop=FALSE]
    colnames(NT) <- gsub("N", "T", colnames(NT))
    out <- cbind(out, NT)
  }
  res <- bayesTest(NP, NT, V)
  out <- as.data.frame(cbind(out, res))
  out[,-which(names(out)=="Type")]
}
