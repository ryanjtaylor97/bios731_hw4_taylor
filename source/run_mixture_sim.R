################################################################################
### AUTHOR: Ryan Taylor
### PURPOSE: Run Simulations for Bayesian K-mixture model
################################################################################


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tictoc))


wd = getwd()

if(substring(wd, 2, 6) == "Users"){
  doLocal = TRUE
}else{
  doLocal = FALSE
}

# Source necessary functions
source(here::here("source", "sim_mix_data.R"))
source(here::here("source", "gibbs_sampler_fn.R"))
source(here::here("source", "variational_inference_fn.R"))

# Set simulation parameters
nsim = 500

n_data = c(100, 1000, 10000)

mu_set <- c(0, 5, 10, 20)



## define number of simulations and parameter scenario
if(doLocal) {
  scenario = 2
  n_scen = n_data[scenario]
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
  n_scen = n_data[scenario]
}

# generate a random seed for each simulated dataset
seed = floor(runif(nsim, 1, 500*4))
results = as.list(rep(NA, nsim))

for(i in 1:nsim){
  set.seed(seed[i])

  ####################
  # simulate data
  simdata = sim_k_gauss(n = n_scen,
                        mus = mu_set,
                        variance = 1)

  ### Apply methods

  # Gibbs
  sim_gibbs <- gibbs_k_gauss(
    y = simdata,
    K = length(mu_set),
    mu0 = 10,
    sigma2 = 100
  )

  # CAVI
  sim_cavi <- cavi_k_gauss(
    y = simdata,
    K = length(mu_set),
    mu0 = 10,
    s2_0 = 10,
    sigma2 = 100,
    tol = 1e-6
  )

  #### Get results

  gibbs_res <- summ_gibbs_mix(sim_gibbs)$mean_summary %>%
    mutate(method = "Gibbs")

  cavi_res <- matrix(
    c(sim_cavi$means,
      sim_cavi$variances,
      sim_cavi$means - 1.96 * sqrt(sim_cavi$variances),
      sim_cavi$means + 1.96 * sqrt(sim_cavi$variances)),
    nrow = 4, byrow = T
  ) %>%
    data.frame() %>%
    mutate(statistic = c("mean", "variance", "2.5%", "97.5%"),
           .before = everything()) %>%
    rename_with(.cols = starts_with("X"),
                ~str_replace(., "X", "group_")) %>%
    mutate(method = "CAVI")

  estimates <- bind_rows(gibbs_res, cavi_res) %>%
    mutate(n = n_scen, seed = seed[i])

  results[[i]] <- estimates
}

## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("results"), Date), showWarnings = FALSE)

filename = paste0(here::here("results", Date), "/", scenario, ".RDA")
save(results,
     file = filename)
