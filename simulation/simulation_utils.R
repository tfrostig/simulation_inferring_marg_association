
# Create genetic data with optional scaling
# Parameters:
# - n: number of samples
# - cov.mat: covariance matrix
# - maf.vec: minor allele frequency vector
# - scale: whether to scale the output data
createGene <- function(n, cov.mat, maf.vec, scale = TRUE) {
  p <- ncol(cov.mat)
  x <- mvrnorm(n, rep(0, p), cov.mat)
  if (scale) x <- scale(x)
  return(x)
}

# Generate correlation matrix with constant off-diagonal elements
# Parameters:
# - p: matrix dimension
# - rho: correlation value
corr2 <- function(p, rho) {
  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  return(Sigma)
}

# Calculate FDR (False Discovery Rate)
# Parameters:
# - sig.ind: indices of significant findings
# - false.null.ind: indices of true positives
calcFDR <- function(sig.ind, false.null.ind) {
  if (length(sig.ind) == 0) return(0)
  sum(!(sig.ind %in% false.null.ind)) / length(sig.ind)
}

# Calculate Power (True Positive Rate)
# Parameters:
# - sig.ind: indices of significant findings
# - false.null.ind: indices of true positives
calcPower <- function(sig.ind, false.null.ind) {
  if (length(false.null.ind) == 0) return(0)
  sum(sig.ind %in% false.null.ind) / length(false.null.ind)
}

# Summarize p-values to calculate FDR and power
# Parameters:
# - pval: vector of p-values
# - qu: significance threshold
# - false.null.ind: indices of true positives
summarisePval <- function(pval, qu, false.null.ind) {
  sig.ind <- which(pval < qu)
  if (length(sig.ind) == 0) {
    fdr <- 0
    power <- 0
  } else {
    fdr <- calcFDR(sig.ind, false.null.ind)
    power <- calcPower(sig.ind, false.null.ind)
  }
  return(c(FDR = fdr, power = power))
}



## Utility function to debug functions
TEST <- function() {
  eff.dim   <<- c(1,20)
  n.org     <<- 10000
  n.ref     <<- 500
  h         <<- 0.005
  qu        <<- 0.05
  #cov.mat  <<- corr2(9, 0.8)
  cov.mat   <<- corr2(20, 0)
  maf.vec   <<- runif(20)
  k.tag     <<- 8
  thres.val <<- qnorm(0.99995, 0, 1 / sqrt(n.org))
  rep.snp   <<- length(maf.vec) %/% 2
}

