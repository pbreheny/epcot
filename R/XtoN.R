XtoN <- function(X, fam) {
  fam$V3 <- splitleft(fam$V3, "-")
  fam$V1 <- splitleft(fam$V1, "_")
  #nonPairs <- which(fam$V1 %in% names(which(table(fam$V1)!=2)))
  #npID <- fam[nonPairs, "V3"]
  Pairs <- which(fam$V1 %in% names(which(table(fam$V1)==2)))
  pID <- fam[Pairs, "V3"]

  ## Convert to 0/1/2 counts
  N <- matrix(NA, nrow(X), 3, dimnames=list(rownames(X), paste0("N", 0:2)))
  pb <- txtProgressBar(1, nrow(X), style=3)
  for (i in 1:nrow(X)) {
    famCounts <- tapply(X[i,pID]!=0, subset(fam, V3 %in% pID)$V1, sum)
    N[i,1] <- sum(famCounts==0)
    N[i,2] <- sum(famCounts==1)
    N[i,3] <- sum(famCounts==2)
    setTxtProgressBar(pb, i)
  }
  N
}
