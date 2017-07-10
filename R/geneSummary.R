geneSummary <- function(Gene, N.pair, N.trio, vData) {
  if (missing(vData)) vData <- get("vData", .GlobalEnv)
  ind <- which(vData$Gene==Gene)
  if (length(ind) == 0) return(NA)
  NP <- N.pair[ind,,drop=FALSE]
  colnames(NP) <- gsub("N", "P", colnames(NP))
  out <- cbind(vData[ind,], RVIS=rvis(Gene), NP)
  if (!missing(N.trio)) {
    NT <- N.trio[ind,,drop=FALSE]
    colnames(NT) <- gsub("N", "T", colnames(NT))
    out <- cbind(out, NT)
  }
  res <- geneTest(Gene, N.pair, N.trio, vData=vData, variants=TRUE, reqExAC=FALSE, maf=1)
  out <- as.data.frame(cbind(out, res))
  out[,-which(names(out)=="Type")]
}
