
# normal_dist_sim.R
# Normal Distribution Simulation

library(MASS)
library(dplyr)
library(magrittr)
library(doParallel)
library(doRNG)
library(ECCCM)
library(CovTools)

# Load simulation utilities
source("simulation/simulation_utils.R")

# Normal Distribution Simulation Function
# Simulates data from multivariate normal distribution and applies different estimators
# Parameters:
# - eff.dim: Indices of causal SNPs
# - n.org: Number of samples in original study
# - n.ref: Number of samples in reference panel
# - h: Heritability
# - qu: Significance threshold
# - cov.mat: Covariance matrix
simNormal <- function(eff.dim, n.org, n.ref, h, qu, cov.mat) {
  cas.dim <- length(eff.dim)
  p <- ncol(cov.mat)

  # Create model coefficients
  beta.vec <- rep(0, p)
  effect.size <- rep(1, p)  # Constant effect size
  beta.vec[eff.dim] <- effect.size[eff.dim]

  # Generate data from multivariate normal
  x.org <- scale(mvrnorm(n.org, rep(0, ncol(cov.mat)), cov.mat))
  x.ref <- scale(mvrnorm(n.ref, rep(0, ncol(cov.mat)), cov.mat))

  # Generate outcome
  if (length(eff.dim) > 0) {
    true.y <- drop(rbind(x.org, scale(x.ref)) %*% beta.vec)
    unexplain.var <- (1 - h) / h
    sigma <- sqrt(unexplain.var)
    y <- true.y + rnorm(n.org + n.ref, mean = 0, sd = sqrt(unexplain.var))
  } else {
    sigma <- 1
    y <- rnorm(n.org + n.ref, mean = 0, sd = sigma)
  }

  # Scale outcome
  sd.y.org <- sd(y)
  y.org <- scale(y[1:n.org])
  sigma <- sigma / sd.y.org

  # Full estimator (oracle)
  inv.o <- solve((t(x.org) %*% x.org) / (n.org - 1))
  beta.org <- solve(t(x.org) %*% x.org) %*% t(x.org) %*% y.org
  var.beta.org <- sigma^2 * diag(inv.o) / n.org
  test.oracle <- testCoef(drop(beta.org), drop(var.beta.org), method = 'BH')

  # Marginal regression coefficients
  beta.report <- drop(solve(diag(apply(x.org, 2, var))) %*% t(x.org) %*% y.org / n.org)

  # Process reference panel
  cov.mat <- cov(x.ref)
  cov.mat.inv <- solve(cov.mat)
  cov.list.r <- list('cov' = cov.mat, 'omega' = cov.mat.inv)

  # Naive estimation
  marg.to.joint <- ECCCM::marginalToJoint(
    marg.beta.hat = beta.report,
    n.o = n.org,
    cor.r = cov2cor(cov.list.r$cov),
    sigma = 1
  )

  test.naive.df <- ECCCM::testCoef(
    est.beta = marg.to.joint$est.beta.hat,
    var.beta = marg.to.joint$naive.var.beta.hat,
    method = 'BH'
  )

  # Empirical estimator
  test.emp.all.est <- ECCCM:::analyzeRef(
    marg.beta.hat = beta.report,
    x.r = x.ref,
    n.o = n.org,
    method.filter = 'BH',
    sigma.method = 'estimate',
    qu = qu,
    to.diag = FALSE
  )$test.correct

  # Gaussian estimator
  test.gauss.emp.all.est <- ECCCM:::analyzeRefGauss(
    marg.beta.hat = beta.report,
    ld.mat = cov.list.r$cov,
    n.o = n.org,
    n.r = n.ref,
    method.filter = 'BH',
    sigma.method = 'estimate',
    qu = qu
  )$test.correct

  # Collate and summarize results
  stat.list <- list(
    'emperical.all' = test.emp.all.est,
    'emperical.gauss' = test.gauss.emp.all.est,
    'naive' = test.naive.df,
    'oracle' = test.oracle
  )

  pval.list <- lapply(stat.list, '[[', 4)
  result.vec <- unlist(lapply(pval.list, summarisePval, qu = qu, false.null.ind = eff.dim))

  return(c(result.vec, 'sigma' = sigma))
}

# Main simulation runner
# Parameters:
# - output_dir: Directory to save results
# - n_iter: Number of simulation iterations
# - parallel: Whether to use parallel processing
# - n_cores: Number of cores for parallel processing
runNormalSim <- function(output_dir = ".",
                         n_iter = 1000,
                         parallel = TRUE,
                         n_cores = 8,
                         seed = 9) {
  # Define parameters
  para.df <- expand.grid(
    'n.org' = c(10000, 20000, 50000, 100000, 200000),
    'n.ref' = c(500, 1000, 5000, 10000),
    'h' = c(0.0005, 0.0025, 0.005, 0.01, 0.05),
    'rho' = c(0.8),
    'p' = c(20),
    'eff.dim' = list(c(1, 5))
  )

  para.df <- para.df %>%
    mutate(cas.snp = sapply(para.df$eff.dim, length))

  # Setup parallel processing
  if (parallel) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    pack.names <- c('MASS', 'dplyr', 'magrittr', 'ECCCM', 'CovTools')
  }

  # Initialize results dataframe
  res.df <- NULL

  # Run simulations for each parameter combination
  for (j in 1:nrow(para.df)) {
    para.vec <- para.df[j, ]
    eff.dim <- unlist(para.vec$eff.dim)
    cov.mat <- corr2(para.vec$p, para.vec$rho)

    # Run simulation iterations
    if (parallel) {
      temp.result <- foreach(
        i = 1:n_iter,
        .packages = pack.names,
        .options.RNG = seed,
        .export = c("simNormal", "testCoef", "summarisePval", "corr2", "calcFDR", "calcPower"),
        .combine = rbind
      ) %dorng% {
        simNormal(
          eff.dim = eff.dim,
          n.org = para.vec$n.org,
          n.ref = para.vec$n.ref,
          h = para.vec$h,
          qu = 0.05,
          cov.mat = cov.mat
        )
      }
    } else {
      temp.result <- matrix(NA, n_iter, 9)
      for (i in 1:n_iter) {
        temp.result[i, ] <- simNormal(
          eff.dim = eff.dim,
          n.org = para.vec$n.org,
          n.ref = para.vec$n.ref,
          h = para.vec$h,
          qu = 0.05,
          cov.mat = cov.mat
        )
      }
    }

    # Calculate mean and SD for results
    mean.result <- apply(temp.result, 2, mean, na.rm = TRUE)
    sd.results <- apply(temp.result, 2, sd, na.rm = TRUE)

    # Append to results dataframe
    res.df <- rbind(res.df, mean.result, sd.results)

    cat('Done with', j, 'out of', nrow(para.df), '\n')
    write.csv(res.df, file.path(output_dir, paste0('COV_RESULTS_NORMAL_', Sys.Date(), '.csv')))
  }

  # Clean up parallel cluster
  if (parallel) {
    stopCluster(cl)
  }

  # Save final results
  para.df <- para.df %>%
    mutate(eff.dim = sapply(para.df$eff.dim, paste, collapse = ','))

  full.res <- data.frame(
    para.df %>% slice(rep(1:n(), each = 2)),
    res.df
  )

  write.csv(full.res, file.path(output_dir, paste0('FULL_RESULTS_NORMAL.csv')))
}
