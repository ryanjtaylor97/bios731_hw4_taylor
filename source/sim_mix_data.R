################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Simulate K-Mixture data
################################################################################

sim_k_gauss <- function(n, mus, variance){

  K <- length(mus)

  # Simulate group assignments
  c_t <- rmultinom(n, 1, prob = rep(1/K, K))

  # Multiply group assignments by means to get individual means
  c_mu <- crossprod(c_t, mus)

  y <- rnorm(n, mean = c_mu, sd = sqrt(variance))

  return(y)
}
