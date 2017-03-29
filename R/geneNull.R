geneNull <- function(maf, N=1000, fams=93) {
  X <- matrix(NA, N, length(maf))

  for (i in 1:N) {
    for (j in 1:length(maf)) {
      p <- rbinom(fams, 2, maf[j])
      m <- rbinom(fams, 2, maf[j])
      d1 <- rbinom(93, 1, p/2) + rbinom(93, 1, m/2)
      d2 <- rbinom(93, 1, p/2) + rbinom(93, 1, m/2)
      X[i,j] <- sum(d1>0 & d2>0)
    }
  }
  apply(X, 1, sum)
}
