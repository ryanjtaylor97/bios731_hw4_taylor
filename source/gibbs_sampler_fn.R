################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Gaussian Mixture Gibbs Sampler functions for BIOS 731 HW 4
################################################################################

##### Gaussian K-Mixture Gibbs Sampler Function
### Takes data vector (y), number of groups (K), and mean variance (sigma2)
### Returns iteration group means (means) and group assignments (group_labels)
gibbs_k_gauss <- function(
    y, K,
    mu0, sigma2,
    iters = 1e4, update_pct = 0.1
){

  tic.clear()
  tic.clearlog()

  ### Start time tracking
  tic("Total")

  tic("Setup")

  ### Define constants for algorithm

  # Iterations to update
  iters_update <- round(seq(update_pct * iters, iters,
                            length.out = 1/update_pct))

  # Determine dimensions
  N <- length(y)

  ### Initialize parameters

  # Parameters to be saved
  mu <- rep_len(mu0, K)
  c_mx <- matrix(0, N, K)

  # Probabilities of group assignments
  c_likes <- matrix(1/K, N, K)

  # Posterior calculations for moments group means
  mu_expect <- mu
  mu_var <- mu
  n_c_vec <- rep(0, K)
  y_c_sum <- rep(0, K)

  ### Create objects to save output

  # Group means
  mean_out <- array(NA, dim = c(K, iters))

  # Group assignments (updated first)
  labs_out <- array(NA, dim = c(N, K, iters))

  toc(log = T, quiet = T)

  # Start loop
  for(iter in 1:iters){

    if(iter %in% iters_update){
      message(paste("On iteration", iter, "at", Sys.time()))
      }

    ### (1) Update group assignments
    tic("Update C")

    # Calculate likelihoods of each group assignment
    for(k in 1:K){
      c_likes[, k] <- dnorm(y, mean = mu[k])
    }

    # Draw new group assignments based on these likelihoods (non-normalized)
    c_mx <- t(
      apply(c_likes, MARGIN = 1,
                  FUN = function(z){ rmultinom(1, 1, z) },
                  simplify = T)
    )

    toc(log = T, quiet = T)

    ### (2) Update group means
    tic("Update Mu")

    # Calculate sum of group memberships
    n_c_vec <- colSums(c_mx)

    # Calculate sum of data in each group
    y_c_sum <- crossprod(c_mx, y)

    # Posterior variances
    mu_var <- 1 / (n_c_vec + 1 / sigma2)

    # Posterior expectations
    mu_expect <- y_c_sum * mu_var

    # Draw new group means
    mu <- rnorm(K, mean = mu_expect, sd = sqrt(mu_var))

    # Sort to avoid label switching
    mu <- sort(mu)

    toc(log = T, quiet = T)

    ### Save results
    mean_out[, iter] <- mu
    labs_out[, , iter] <- c_mx
  }

  # Save timing
  toc(log = T, quiet = T)

  time_log <- tic.log(format = F)

  timing <- bind_rows(
    lapply(time_log,
           function(x){
             tibble(step = x$msg,
                    time = x$toc - x$tic)
           }
    ))

  # Save output
  return(list(
    "iters" = iters,
    "means" = mean_out,
    "group_labels" = labs_out,
    "time" = timing
  ))
}

##### Gibbs Sampler Summary Function
### Takes output from Gaussian K-Mixture Model above and desired burn-in %
### Returns posterior estimates of group means and group assignments

summ_gibbs_mix <- function(model, burn_pct = 0.2){

  iters_total <- model$iters

  iters_burn <- ceiling(model$iters * burn_pct)

  # Take posterior means of group means
  mu_mean <- apply(model$means[, iters_burn:iters_total],
                   1,
                   FUN = function(x) matrix(
                     c(mean(x),
                       var(x),
                       quantile(x, 0.025),
                       quantile(x, 0.975)),
                     ncol = 1
                   ),
                   simplify = T)

  # Convert to formatted table
  mu_mean_tbl <- data.frame(mu_mean) %>%
    mutate(statistic = c("mean", "variance", "2.5%", "97.5%"),
           .before = everything()) %>%
    rename_with(.cols = starts_with("X"),
                ~str_replace(., "X", "group_"))

  # Take posterior means of group labels
  labels_pct <- apply(model$group_labels[,, iters_burn:iters_total],
                      c(1,2), mean,
                      simplify = T)

  labels_tbl <- data.frame(labels_pct) %>%
    mutate(row_num = 1:n(), .before = everything()) %>%
    rename_with(.cols = starts_with("X"),
                ~str_replace(., "X", "group_"))

  # Identify most common group label for each observation
  group_assignments <- apply(labels_pct, 1, which.max,
                             simplify = T)

  return(list(
    "mean_summary" = mu_mean_tbl,
    "group_labels" = group_assignments,
    "group_pcts" = labels_tbl
  ))
}
