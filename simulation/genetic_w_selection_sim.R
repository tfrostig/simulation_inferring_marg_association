
# genetic_w_selection_sim.R
# PSAT Genetic Simulation

library(MASS)
library(dplyr)
library(magrittr)
library(doParallel)
library(doRNG)
library(ECCCM)
library(CovTools)
library(PSATinference)

# Load simulation utilities
source("simulation/simulation_utils.R")
# Number of genes considered for threshold (to replicate results from manuscript, change to 20000)
GENE_NUMBER = 20000
# PSAT Genetic Simulation Function
# Simulates genetic data with selection and applies various post-selection inference methods
# Parameters:
# - eff.dim: Indices of causal SNPs
# - n.org: Number of samples in  original study
# - n.ref: Number of samples in reference panel
# - h: Heritability
# - qu: Significance threshold
# - cov.mat: Covariance matrix
# - maf.vec: Minor allele frequency vector
# - rep.snp: Index of SNP used for selection
# - thres.val: Threshold value for selection
simGeneticPSAT <- function(eff.dim, n.org, n.ref, h, qu, cov.mat, maf.vec,
                           rep.snp = ncol(cov.mat) %/% 2,
                           thres.val = qnorm(1 - qu/GENE_NUMBER, 0, 1 / sqrt(n.org)) ) {
  # Define result names
  result.name <- c(
    "Oracle.Corrected_polyhedral.symmetric.FDR", "Oracle.Corrected_polyhedral.symmetric.power",
    "Oracle.Corrected_naive.FDR", "Oracle.Corrected_naive.power",
    "Oracle.Corrected_global.FDR", "Oracle.Corrected_global.power",
    "Oracle.Corrected_hybrid.symmetric.FDR", "Oracle.Corrected_hybrid.symmetric.power",

    "MLE.Corrected_polyhedral.symmetric.FDR", "MLE.Corrected_polyhedral.symmetric.power",
    "MLE.Corrected_naive.FDR", "MLE.Corrected_naive.power",
    "MLE.Corrected_global.FDR", "MLE.Corrected_global.power",
    "MLE.Corrected_hybrid.symmetric.FDR", "MLE.Corrected_hybrid.symmetric.power",

    "Corrected_polyhedral.symmetric.FDR", "Corrected_polyhedral.symmetric.power",
    "Corrected_naive.FDR", "Corrected_naive.power",
    "Corrected_global.FDR", "Corrected_global.power",
    "Corrected_hybrid.symmetric.FDR", "Corrected_hybrid.symmetric.power",

    "Naive_polyhedral.symmetric.FDR", "Naive_polyhedral.symmetric.power",
    "Naive_naive.FDR", "Naive_naive.power",
    "Naive_global.FDR", "Naive_global.power",
    "Naive_hybrid.symmetric.FDR", "Naive_hybrid.symmetric.power",

    "Oracle_polyhedral.symmetric.FDR", "Oracle_polyhedral.symmetric.power",
    "Oracle_naive.FDR", "Oracle_naive.power",
    "Oracle_global.FDR", "Oracle_global.power",
    "Oracle_hybrid.symmetric.FDR", "Oracle_hybrid.symmetric.power"
  )

  print(thres.val)
  cas.dim <- length(eff.dim)
  p       <- ncol(cov.mat)

  # Create model coefficients
  beta.vec    <- rep(0, p)
  effect.size <- rep(1, p) ## Constant effect size
  beta.vec[eff.dim] <- effect.size[eff.dim]

  # Generate data
  x.org <- createGene(n.org, cov.mat, maf.vec, scale = TRUE)
  x.ref <- createGene(n.ref, cov.mat, maf.vec, scale = FALSE)

  # Define selection vector
  e.snp      <- rep(0, p)
  e.snp[rep.snp] <- 1

  # Generate outcome with selection until threshold is met
  number.to.pass <- 0
  select.val <- 0
  keep_val <- c()

  true.y <- drop(rbind(x.org, scale(x.ref)) %*% beta.vec)
  unexplain.var <- ((1 - h) / h) * (cas.dim > 0)
  sigma <- sqrt(unexplain.var)

  transpose_rep_x <- t(x.org[ ,rep.snp])

  while (select.val^2 < thres.val^2) {
    # Generate outcome
    y <- true.y + rnorm(n.org + n.ref, mean = 0, sd = sigma)

    # Scale outcome
    y.org    <- scale(y[1:n.org])

    # Compute marginal coefficients
    select.val  <- transpose_rep_x %*% y.org / n.org

    # Count selection attempts
    number.to.pass <- number.to.pass + 1
  }

  cat("number.to.pass:", number.to.pass, "\n")
  beta.report <- t(x.org) %*% y.org / n.org

  # Adjust sigma for scaling
  sd.y.org <- sd(y[1:n.org])
  sigma <- sigma / sd.y.org

  # Compute oracle estimator
  inv.o <- solve((t(x.org) %*% x.org) / (n.org - 1))
  beta.org <- solve(t(x.org) %*% x.org) %*% t(x.org) %*% y.org
  var.beta.org <- sigma^2 * diag(inv.o) / n.org
  test.oracle <- testCoef(drop(beta.org), drop(var.beta.org), method = 'BH')

  # Compute reference panel correlation
  cor.mat <- cor(x.ref)
  cor.mat.inv <- solve(cor.mat)
  cor.list.r <- list('cor' = cor.mat, 'cor_omega' = cor.mat.inv)

  # Joint coefficient transformation
  marg.to.joint <- ECCCM::marginalToJoint(
    marg.beta.hat = beta.report,
    n.o = n.org,
    cor.r = cor.list.r$cor,
    sigma = 1
  )

  if (select.val^2 > thres.val^2) {
    # Selection matrices
    k.o <- 1/n.org^2 * t(x.org) %*% x.org %*% e.snp %*% t(e.snp) %*% t(x.org) %*% x.org
    k.r <- cor.mat %*% e.snp %*% t(e.snp) %*% cor.mat

    # Run PSAT analyses with different variance estimators
    post.selection.naive.val <- PSATQuadratic(
      y = marg.to.joint$est.beta.hat,
      cov.mat = (1/n.org) * cor.mat.inv,
      K = k.r,
      test.type = 'symmetric',
      pval.type = c('polyhedral', 'naive', 'hybrid'),
      est.type = c('naive', 'mle'),
      ci.type = 'naive',
      alpha = qu,
      threshold = thres.val^2
    )

    # MLE variance correction
    correct.var.mle <- estimateAdditionalVariance(
      coeff = post.selection.naive.val$Point.estimation$mle,
      x.ref = x.ref,
      n.org = n.org,
      qu = qu,
      cor.list = cor.list.r,
      is.marg = FALSE
    )

    # Oracle variance correction
    correct.var.oracle <- estimateAdditionalVariance(
      coeff = beta.vec / sd.y.org,
      x.ref = x.ref,
      n.org = n.org,
      qu = qu,
      cor.list = cor.list.r,
      is.marg = FALSE
    )

    # Empirical variance correction
    correct.var.old <- estimateAdditionalVariance(
      coeff = beta.report,
      x.ref = x.ref,
      n.org = n.org,
      qu = qu,
      cor.list = cor.list.r,
      is.marg = TRUE
    )

    # PSAT with MLE correction
    post.selection.empirical.mle <- PSATQuadratic(
      y = marg.to.joint$est.beta.hat,
      cov.mat = correct.var.mle$var,
      K = k.r,
      test.type = 'symmetric',
      pval.type = c('polyhedral', 'naive', 'hybrid'),
      est.type = c('naive'),
      ci.type = 'naive',
      alpha = qu,
      threshold = thres.val^2
    )

    # PSAT with oracle correction
    post.selection.empirical.oracle <- PSATQuadratic(
      y = marg.to.joint$est.beta.hat,
      cov.mat = correct.var.oracle$var,
      K = k.r,
      test.type = 'symmetric',
      pval.type = c('polyhedral', 'naive', 'hybrid'),
      est.type = c('naive'),
      ci.type = 'naive',
      alpha = qu,
      threshold = thres.val^2
    )

    # PSAT with empirical correction
    post.selection.empirical.val <- PSATQuadratic(
      y = marg.to.joint$est.beta.hat,
      cov.mat = correct.var.old$var,
      K = k.r,
      test.type = 'symmetric',
      pval.type = c('polyhedral', 'naive', 'hybrid'),
      est.type = c('naive'),
      ci.type = 'naive',
      alpha = qu,
      threshold = thres.val^2
    )

    # Naive PSAT (different parameterization)
    post.selection.naive.val <- PSATQuadratic(
      y = marg.to.joint$est.beta.hat,
      cov.mat = drop(correct.var.old$sigma^2/n.org) * cor.mat.inv,
      K = k.r,
      test.type = 'symmetric',
      pval.type = c('polyhedral', 'naive', 'hybrid'),
      est.type = c('naive'),
      ci.type = 'naive',
      alpha = qu,
      threshold = thres.val^2
    )

    # Oracle PSAT
    post.selection.oracle.val <- PSATQuadratic(
      y = drop(beta.org),
      cov.mat = drop(sigma)^2 * inv.o / n.org,
      K = k.o,
      test.type = 'symmetric',
      pval.type = c('polyhedral', 'naive', 'hybrid'),
      est.type = 'naive',
      ci.type = 'naive',
      alpha = qu,
      threshold = thres.val^2
    )

    # Combine all p-values
    pval.df <- bind_cols(
      post.selection.empirical.oracle$Pvalues,
      post.selection.empirical.mle$Pvalues,
      post.selection.empirical.val$Pvalues,
      post.selection.naive.val$Pvalues,
      post.selection.oracle.val$Pvalues
    )

    # Adjust for multiple testing
    pval.df <- apply(pval.df, 2, p.adjust, method = 'BH') %>% as.data.frame()

    # Calculate FDR and power
    result.vec <- unlist(lapply(pval.df, summarisePval, qu = qu, false.null.ind = eff.dim))
  }

  if (!exists("result.vec")) {
    result.vec <- rep(NA, 40)
  }

  names(result.vec) <- result.name
  return(c(result.vec,
           'is_pass' = (select.val^2 > thres.val^2),
           'statistic' = select.val,
           'sigma' = sigma,
           'number.to.pass' = number.to.pass))
}

# Main simulation runner
# Parameters:
# - output_dir: Directory to save results
# - n_iter: Number of simulation iterations
# - parallel: Whether to use parallel processing
# - n_cores: Number of cores for parallel processing
runGeneticPSATSim <- function(output_dir = ".",
                              n_iter = 10,
                              parallel = TRUE,
                              n_cores = 8) {
  # Define parameters
  para.df <- expand.grid(
    'n.org' = c(10000),
    'n.ref' = c(500, 1000, 2000),
    'h' =c(0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2),
    'rho' = c(0.7, 0.85, 0.95),
    'p' = c(20),
    'eff.dim' = list(c(1, 20), c(9, 10, 11)),
    'g.sparse' = c('rare'),
    'n.iter' = n_iter
  )

  para.df <- para.df %>%
    mutate(cas.snp = sapply(para.df$eff.dim, length))

  # Generate MAF vectors
  maf.vec.1 <- rbeta(max(para.df$p), 1, 2) / 2
  maf.vec.1 <- pmax(maf.vec.1, 0.05)  # Ensure no nulls
  maf.vec.2 <- runif(max(para.df$p), 0.05, 0.5)

  # Setup parallel processing
  if (parallel) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    pack.names <- c('MASS', 'dplyr', 'magrittr', 'ECCCM', 'CovTools', 'PSATinference')
  }

  # Initialize results dataframe
  res.df <- NULL

  # Run simulations for each parameter combination
  for (j in 1:nrow(para.df)) {
    para.vec <- para.df[j, ]
    eff.dim <- unlist(para.vec$eff.dim)
    cov.mat <- corr2(para.vec$p, para.vec$rho)

    # Select MAF vector based on sparsity
    if (para.vec$g.sparse == 'rare') {
      maf.vec <- maf.vec.1
    } else {
      maf.vec <- maf.vec.2
    }
    maf.vec.temp <- maf.vec[1:para.vec$p]

    # Run simulation iterations
    if (parallel) {
      temp.result <- foreach(
        i = 1:n_iter,
        .packages = pack.names,
        .options.RNG = 9,
        .combine = rbind,
        .export = c("simGeneticPSAT", "createGene", "testCoef", "summarisePval",
                    "corr2", "calcFDR", "calcPower", "estimateAdditionalVariance",
                    "GENE_NUMBER")
      ) %dorng% {
        simGeneticPSAT(
          eff.dim = eff.dim,
          n.org = para.vec$n.org,
          n.ref = para.vec$n.ref,
          h = para.vec$h,
          qu = 0.05,
          cov.mat = cov.mat,
          maf.vec = maf.vec.temp,
          thres.val = qnorm(1 - 0.05 / GENE_NUMBER, 0, para.vec$n.org^(-0.5))
        )
      }
    } else {
      temp.result <- matrix(NA, n_iter, 44)  # Adjusted for additional outputs
      for (i in 1:n_iter) {
        temp.result[i, ] <- simGeneticPSAT(
          eff.dim = eff.dim,
          n.org = para.vec$n.org,
          n.ref = para.vec$n.ref,
          h = para.vec$h,
          qu = 0.05,
          cov.mat = cov.mat,
          maf.vec = maf.vec.temp,
          thres.val = qnorm(1 - 0.05 / GENE_NUMBER, 0, para.vec$n.org^(-0.5)))
      }
    }
    # Calculate mean and SD for results
    print(temp.result)
    mean.result <- apply(temp.result, 2, mean, na.rm = TRUE)
    sd.results <- apply(temp.result, 2, sd, na.rm = TRUE)

    # Append to results dataframe
    res.df <- rbind(res.df, mean.result, sd.results)

    cat('Done with', j, 'out of', nrow(para.df), '\n')
    write.csv(res.df, file.path(output_dir, paste0('COV_RESULTS_PSAT_', Sys.Date(), '.csv')))
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

  write.csv(full.res, file.path(output_dir, paste0('FULL_RESULTS_PSAT.csv')))

  return(full.res)
}


# Estimate additional variance for coefficients
# coeff: Coefficient estimates
# x.ref: Reference panel
# n.org: Original sample size
# qu: Significance threshold
# cor.list: List with correlation matrix and its inverse
# is.marg: Boolean, are coefficients marginal
# Returns: List with variance and sigma
estimateAdditionalVariance <- function(coeff, x.ref, n.org, qu, cor.list, is.marg) {
  add.var <- 0
  test.emp.all.est <- ECCCM:::analyzeRef(marg.beta.hat = coeff,
                                         x.r = x.ref,
                                         n.o = n.org,
                                         method.filter = 'BH',
                                         sigma.method = 'estimate',
                                         qu = qu,
                                         to.diag = FALSE,
                                         is.marg = is.marg)
  if (!is.null(test.emp.all.est$add.var)) {
    add.var <- test.emp.all.est$add.var
  }
  correct.var <- test.emp.all.est$sigma^2 * (cor.list$cor_omega / n.org) + add.var
  return(list('var' = correct.var,
              'sigma' = test.emp.all.est$sigma))
}
