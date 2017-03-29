XtoN <- function(X, fam) {
  Fam <- splitleft(fam[,1], "_")
  nonPairs <- which(Fam %in% names(which(table(Fam)!=2)))

  ## Convert to 0/1/2 counts
  N <- matrix(NA, nrow(X), 3, dimnames=list(rownames(X), 0:2))
  pb <- txtProgressBar(1, nrow(X), style=3)
  for (i in 1:nrow(X)) {
    famCounts <- tapply(X[i,-nonPairs]!=0, Fam[-nonPairs], sum)
    N[i,1] <- sum(famCounts==0)
    N[i,2] <- sum(famCounts==1)
    N[i,3] <- sum(famCounts==2)
    setTxtProgressBar(pb, i)
  }
  N
}
