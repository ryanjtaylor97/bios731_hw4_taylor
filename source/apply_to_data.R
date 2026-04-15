################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run Gibbs Sampler and CAVI Function on Old Faithful Data
################################################################################


# Load data ---------------------------------------------------------------

data("faithful")

# Scale so that the data have unit variance
scale_factor <- sd(faithful$waiting)

# Extract objects for algorithms
y_vec <- faithful$waiting / scale_factor

K_expected <- 2


# Run Gibbs Sampler -------------------------------------------------------

set.seed(978)

# Run 1 Gibbs chain on Old Faithful data
of_gibbs_100 <- gibbs_k_gauss(
  y = y_vec,
  K = K_expected,
  mu0 = 70 / scale_factor,
  sigma2 = 100
)

of_gibbs_100_summ <- summ_gibbs_mix(of_gibbs_100)

set.seed(617)

# Run second Gibbs chain with lower proposed variance
of_gibbs_10 <- gibbs_k_gauss(
  y = y_vec,
  K = K_expected,
  mu0 = 70 / scale_factor,
  sigma2 = 10
)

of_gibbs_10_summ <- summ_gibbs_mix(of_gibbs_10)

set.seed(508)

# Run third Gibbs chain with even lower proposed variance
of_gibbs_1 <- gibbs_k_gauss(
  y = y_vec,
  K = K_expected,
  mu0 = 70 / scale_factor,
  sigma2 = 1
)

of_gibbs_1_summ <- summ_gibbs_mix(of_gibbs_1)


# Run CAVI ----------------------------------------------------------------

# Run CAVI function on Old Faithful data
of_cavi <- cavi_k_gauss(
  y = y_vec,
  K = K_expected,
  mu0 = c(50, 80) / scale_factor,
  s2_0 = c(10, 10) / scale_factor,
  sigma2 = 100,
  tol = 1e-6
)
