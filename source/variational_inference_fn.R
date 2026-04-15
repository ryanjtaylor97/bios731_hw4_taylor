################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Variational Inference for Bayesian K-Mixture Model in BIOS 731
################################################################################

##### Variational Inference Function
### Takes data (y), number of groups (K), and mean variance (sigma2)
### Returns estimated group means and group membership likelihoods
cavi_k_gauss <- function(
    y, K,
    mu0, s2_0, sigma2,
    max_iters = 1e3, update_pct = 0.1, tol = 1e-6
){

  tic.clear()
  tic.clearlog()

  ### Start time tracking
  tic("Total")

  tic("Setup")

  ### Define constants for algorithm

  # Iterations to update
  iters_update <- round(seq(0, max_iters, length.out = 1/update_pct))

  # Determine dimensions
  N <- length(y)

  ### Initialize parameters

  # Parameters to be saved
  m_vec <- rep_len(mu0, K)
  s2_vec <- rep(1, K)
  psi_mx <- matrix(dnorm(rnorm(N*K)), N, K)

  ### Create objects to save output

  # ELBO storage
  elbo_vec <- rep(NA, max_iters)
  elbo_likes <- rep(NA, K)
  converge_count <- 0

  toc(log = T, quiet = T)

  # Start loop
  for(iter in 1:max_iters){

    if(iter %in% iters_update){
      message(paste("On iteration", iter, "at", Sys.time()))
    }

    ### (1) Update group assignment probabilities
    tic("Update Psi")

    for(k in 1:K){
      psi_mx[, k] <- exp( m_vec[k] * y - 0.5 * (s2_vec[k] + m_vec[k]^2) )
    }

    # Normalize these probabilities to sum to 1
    psi_sums <- apply(psi_mx, 1, sum, simplify = T)

    psi_mx <- sweep(psi_mx, 1, psi_sums, FUN = "/")

    toc(log = T, quiet = T)

    ### (2) Update group means
    tic("Update mean q")

    # Calculate denominator common to both moments
    mu_denom <- 1 / sigma2 + colSums(psi_mx)

    # Calculate group mean vector
    m_vec <- crossprod(psi_mx, y) / mu_denom

    # Calculate group mean variance vector
    s2_vec <- 1 / mu_denom

    toc(log = T, quiet = T)

    ### Calculate ELBO
    tic("ELBO")

    # Start with likelihood, since it is complicated
    for(k in 1:K){
      elbo_likes[k] <- sum(
        -0.5 * psi_mx[,k] * (
          s2_vec[k] + (y - m_vec[k])^2
        )
      )
    }

    # Sum along with other parts
    elbo_vec[iter] <-
      sum(elbo_likes) -
      (1 / (2 * sigma2)) * sum(m_vec^2 + s2_vec) -
      sum(psi_mx * log(psi_mx)) +
      0.5 * sum( log(2 * pi * s2_vec) )

    toc(log = T, quiet = T)

    ### Check convergence
    if (iter > 1 && abs(elbo_vec[iter] - elbo_vec[iter-1]) < tol){
      converge_count <- converge_count + 1
    } else{
      converge_count <- 0
    }

    if(converge_count >= 10){ break }

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

  return(list(
    "means" = m_vec,
    "variances" = s2_vec,
    "group_probs" = psi_mx,
    "elbo" = elbo_vec[seq_len(iter)],
    "time" = timing
  ))
}


