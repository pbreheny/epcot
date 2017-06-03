showGene <- function(Gene, N, vData, cutoff=0.05) {
  ind <- which(vData$Gene==Gene & vData$ExAC < cutoff)
  cbind(vData[ind,], N[ind,,drop=FALSE])
}
