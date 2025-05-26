
# Calculate conditional coefficient distribution given marginal coefficients
# Parameters:
# - marg_coefs: Marginal regression coefficients
# - x_ref: Reference panel data
# - n_org: Original sample size
# - inds_in_question: Indices to test
# - inds_to_condition: Indices to condition on
# - var_type: Variance type ('conservative' or 'estimate')
conditonal_coef_dist_given_marg_coefs <- function(marg_coefs,
                                                  x_ref,
                                                  n_org,
                                                  inds_in_question,
                                                  inds_to_condition = NULL,
                                                  var_type = 'conservative') {
  n_r <- nrow(x_ref)
  x_2 <- x_ref[, inds_in_question, drop = FALSE]
  inv_x_2 <- solve(t(x_2) %*% x_2 / (nrow(x_2) - 1))
  beta_2 <- marg_coefs[inds_in_question]

  # Case with no conditioning
  if (is.null(inds_to_condition)) {
    cond_est <- inv_x_2 %*% beta_2
    cond_var <- inv_x_2 / n_org
    p_value <- 2 * (1 - pnorm(abs(cond_est), sd = sqrt(diag(cond_var))))
    return(data.frame(cond_est = cond_est, cond_var = diag(cond_var), p_value = p_value))
  }

  # Extract matrices for conditioning
  x_1 <- x_ref[, inds_to_condition, drop = FALSE]
  inv_x_1 <- solve(t(x_1) %*% x_1 / (nrow(x_1) - 1))
  beta_1 <- marg_coefs[inds_to_condition]

  # Check for overlap
  if (any(inds_in_question %in% inds_to_condition)) {
    stop('Indices in question overlap with conditional indices')
  }

  # Cross-products
  t_x2_x1 <- t(x_2) %*% x_1 / n_r

  # Estimate conditional effects
  cond_est <- n_r / n_org * inv_x_2 %*% beta_2 -
    n_r / n_org * inv_x_2 %*% t_x2_x1 %*% inv_x_1 %*% beta_1

  # Calculate variance
  var_resid_y <- 1  # Conservative estimate

  # Compute the conditional variance
  cond_var <- var_resid_y * inv_x_2 / n_org -
    (var_resid_y / n_org) * inv_x_2 %*% t_x2_x1 %*% inv_x_1 %*% t(t_x2_x1) %*% inv_x_2

  # Calculate p-values
  p_value <- 2 * (1 - pnorm(abs(cond_est), sd = sqrt(diag(cond_var))))

  # Return results
  return(data.frame(cond_est = cond_est, cond_var = diag(cond_var), p_value = p_value))
}


# Forward selection algorithm
# Parameters:
# - marg_coefs: Marginal regression coefficients
# - x_ref: Reference panel data
# - n_org: Original sample size
# - var_type: Variance estimation method
# - max_steps: Maximum number of steps
# - p_value_threshold: P-value threshold for selection
# - selected_indices: Initially selected indices
forward_selection <- function(marg_coefs,
                              x_ref,
                              n_org,
                              var_type = 'conservative',
                              max_steps = NULL,
                              p_value_threshold = 0.05,
                              selected_indices = NULL) {
  n_pred <- length(marg_coefs)
  scaled_x_ref <- scale(x_ref)
  remaining_indices <- 1:n_pred  # all indices are initially available

  # Initialize selected indices
  if (is.null(selected_indices)) {
    selected_indices <- c()
  } else {
    remaining_indices <- setdiff(remaining_indices, selected_indices)
  }

  # Storage for results
  best_cond_estimates <- c()
  best_cond_vars <- c()
  best_cond_p_values <- c()

  # Set default max_steps
  if (is.null(max_steps)) {
    max_steps <- n_pred
  }

  # Forward selection loop
  for (step in 1:max_steps) {
    best_p_value <- Inf
    best_index <- NULL
    best_cond_var <- NULL
    best_cond_est <- NULL

    # Test adding each remaining index
    for (index in remaining_indices) {
      result <- conditonal_coef_dist_given_marg_coefs(
        marg_coefs = marg_coefs,
        x_ref = scaled_x_ref,
        n_org = n_org,
        inds_in_question = index,
        inds_to_condition = selected_indices,
        var_type = var_type
      )

      current_p_value <- result$p_value[length(result$p_value)]

      # Update if better p-value found
      if (current_p_value < best_p_value) {
        best_p_value <- current_p_value
        best_index <- index
        best_cond_var <- result$cond_var[length(result$p_value)]
        best_cond_est <- result$cond_est[length(result$p_value)]
      }
    }

    # Add best index to selected indices
    selected_indices <- c(selected_indices, best_index)

    # Jointly fit the selected variables
    joint_fit <- ECCCM::marginalToJoint(
      marg.beta.hat = marg_coefs[selected_indices],
      n.o = n_org,
      cor.r = cov(scaled_x_ref)[selected_indices, selected_indices, drop = FALSE],
      sigma = 1
    )

    # Add indices to joint fit
    joint_fit['index'] <- selected_indices
    joint_fit['p_value'] <- 2 * (1 - pnorm(abs(joint_fit$est.beta.hat / sqrt(joint_fit$naive.var.beta.hat))))

    # Filter based on p-value threshold
    selected_indices <- joint_fit$index[joint_fit$p_value < p_value_threshold]

    # Update remaining indices
    remaining_indices <- setdiff(remaining_indices, selected_indices)

    # Store results
    best_cond_estimates <- c(best_cond_estimates, best_cond_est)
    best_cond_vars <- c(best_cond_vars, best_cond_var)
    best_cond_p_values <- c(best_cond_p_values, best_p_value)

    # Stop if no improvement found
    if (best_p_value > p_value_threshold) {
      break
    }
  }

  # Return results
  results <- list(
    selected_indices = selected_indices,
    best_cond_estimates = joint_fit$est.beta.hat[joint_fit$p_value < p_value_threshold],
    cond_var = joint_fit$naive.var.beta.hat[joint_fit$p_value < p_value_threshold],
    p_value = joint_fit$p_value[joint_fit$p_value < p_value_threshold]
  )

  return(results)
}
