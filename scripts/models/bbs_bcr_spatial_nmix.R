---
title: "Integrated Spatial Dynamic N-mixture model"
subtitle: Modeling abundance of birds from eBird, BBS, and IMBCR
author: Alex Baecher
format:
  html:
    theme: cosmo
    margin: 10px
    toc: true
    toc-location: left
    toc-depth: 3
    smooth-scroll: true
    highlight-style: github 
---
#` Code for an integrated spatial dynamic N-mixture model using eBird, BBS, and IMBCR data. 
#` This script runs a spatial model across a stratified grid of ~16k cells in the western USA. 
#` Currently, this model runs extremely slow. Consider applying a [paralellization](https://r-nimble.org/examples/parallelizing_NIMBLE.html) approach.
#` Basic process and observation submodels were adapted from [Zhao et al. MEE 2024](https://doi.org/10.1111/2041-210X.14368). 

# Load packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())
pacman::p_load(
  tidyverse
  , here
  , sf
  , spdep 
  , fields
  , nimble
  , mgcv
  , mvgam
  # , 
  # , 
)


# Gridded study region ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load 20km x 20km gridded sf object of study area 
strat <- read_rds("data/region/strat.rds")

# Loading BBS data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load BBS bird data stratified by the `strat` object 
bbs_strat <- read_rds("data/bbs/species/strat/ybma.rds") %>%
  mutate(is_observed = !is.na(total_count) & total_count > 0) %>%   # flag observed and unobserved sites 
  arrange(cell_id, time_step) 

## Filter stratified BBS data to only OBSERVED sites
bbs_strat_obs <- bbs_strat %>%
  filter(is_observed == T) %>% 
  group_by(cell_id) %>%
  mutate(obs_site_id = cur_group_id()) %>%                          # create new site ID 
  ungroup() %>%
  group_by(route_id) %>%
  mutate(route_site_id = cur_group_id()) %>%                        # create sampling route ID
  select(
    state, route_id, route_site_id, obs_site_id, year, time_step, 
    count1, count2, count3, count4, count5, 
    count6, count7, count8, count9, count10, 
    total_count, obs_id, new, nyears_site, geometry
  )



ggplot() + 
  geom_sf(data = strat) + 
  geom_sf(data = bbs_strat_obs, aes(fill = total_count)) + 
  theme_bw()


# Create spatial knots for prediction with reduced rank Gaussian process ~~~~~~~~~~~~~

## Get coordinates for ALL sites in strat (this is what we'll predict to)
strat_coords <- strat %>% 
  st_centroid() %>%
  st_coordinates()

## Get observed site coordinates (for knot selection)
strat_obs_coords <- bbs_strat_obs %>%     # add EBD and BCR data here! 
  st_centroid() %>%
  st_coordinates()

## Select knot locations based on observed sites
n_knots <- 200  # Adjust based on your spatial domain
set.seed(123)  # for reproducibility
knot_indices <- kmeans(strat_obs_coords, centers = n_knots)$cluster
knot_coords <- strat_obs_coords[!duplicated(knot_indices), ]

## Calculate distance between all sites and knots (for distance-decay calculation)
dist_all_knots <- rdist(strat_coords, knot_coords)  # ALL sites to knots


## Function to compute covariance with fixed parameters (this was for a previous spatial structure)

# compute_cov <- function(dist_mat, sigma2, phi) {
#   sigma2 * exp(-dist_mat / phi)
# }

# ## Use initial values for `sigma2_s` and `phi_s`
# sigma2_init <- 1
# phi_init <- 50  # Adjust based on your spatial scale

# ## Compute covariance matrices
# sigma_knots <- compute_cov(dist_knots, sigma2_init, phi_init)
# chol_sigma_knots <- chol(sigma_knots)

# cov_all_knots <- compute_cov(dist_all_knots, sigma2_init, phi_init)
# c_all_knots <- cov_all_knots %*% solve(sigma_knots)  # Projection matrix for ALL sites


# Create model inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Observed data

data <- list(
  bbs_cnt = bbs_strat_obs %>%
    select(count1:count10) %>% 
    st_drop_geometry() %>%
    as.matrix()
  # , EBD data
  # , BCR data
)


## Constants

constants <- list(
  ### Universal constants
  nsite = nrow(strat)   # total sites in gridded region

  ### Gaussian processes constants
  # , mu_knots = rep(0, n_knots)
  # , chol_sigma_knots = chol_sigma_knots
  # , c_all_knots = c_all_knots  # Projection matrix for all sites

  ### Alternative spatial structure
  , n_knots = n_knots
  # , dist_knots = dist_knots             # alternative spatial structure
  , dist_all_knots = dist_all_knots       # alternative spatial structure 
  ### BBS constants
  , bbs_max_years = max(bbs_strat_obs$time_step)
  , bbs_nobs = nrow(bbs_strat_obs)
  , bbs_sites = bbs_strat_obs$obs_site_id
  , bbs_time_step = bbs_strat_obs$time_step
  , bbs_nroute = length(unique(bbs_strat_obs$route_site_id))
  , bbs_route = bbs_strat_obs$route_site_id
  , bbs_nobser = length(unique(bbs_strat_obs$obs_id))
  , bbs_obser = bbs_strat_obs$obs_id
  , bbs_nstop = ncol(data$bbs_cnt)
  , bbs_stops = seq(-2, 2, length.out = ncol(data$bbs_cnt))
  , bbs_new = bbs_strat_obs$new
  ### EBD constants
  ### BCR constants
)


## Initialize N matrix for ALL sites

Ni <- matrix(
  10, 
  constants$nsite, 
  constants$bbs_max_years)


### Only set higher initial values for observed sites

for(i in unique(bbs_strat_obs$obs_site_id)){
  site_data <- bbs_strat_obs[bbs_strat_obs$obs_site_id == i, ]
  if(nrow(site_data) > 0) {
    for(j in 1:nrow(site_data)){
      time_step <- site_data$time_step[j]
      obs_counts <- data$bbs_cnt[which(
        bbs_strat_obs$obs_site_id == i & bbs_strat_obs$time_step == time_step), ]
      if(length(obs_counts) > 0){
        Ni[i, time_step] <- max(max(obs_counts, na.rm = TRUE), 10) * 3
      }
    }
  }
}


## Initialize z for all sites

zi <- ifelse(rowSums(Ni) == 10 * constants$bbs_max_years, 0, 1)


## Initial values 

inits <- list(
  # Universal inits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process model inits 
    zeta = 0.5
  , z = zi
  , N = Ni
  , beta_lambda0 = 0.5
  , beta_rho = -0.05  
  
  # Spatial Gaussian process inits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  , w = rep(0, n_knots)
  , sigma2_s = 1
  , phi_s = 50
  # , tau2_s = 0.1                 # spatial nugget (site-level variation)
  # , s = rep(0, constants$nsite)  # Spatial effects for ALL sites

  # Observation submodels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # BBS inits 
  , chi_mu = -5
  , chi_route_sd = 0.5
  , chi_obser_sd = 0.5
  , chi_new = 0
  , chi_stop_beta = c(-1, 0.1)
  , chi_route_epsilon = rep(0, constants$bbs_nroute)
  , chi_obser_epsilon = rep(0, constants$bbs_nobser)
  , log_chi = matrix(-5, constants$bbs_nobs, constants$bbs_nstop)
    # EBD constants
  , omega_beta = c(-6,0.1,0.5,0.2,0), 
  , omega_locat_sd = 0.02, 
  , omega_obser_sd = 1, 
  , omega_zeta = 0.3, 
  , omega_z = ifelse(ebird_cnt>0,1,0), 
  , omega_delta = 0.01, 
    # BCR constants
)

# Write model code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

code <- nimbleCode({ 
  # Priors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process submodel priors: 
  zeta ~ dunif(0, 1)                                           # true abundance prior
  beta_lambda0 ~ dnorm(0, sd = 10)                             # abundance fixed effect prior
  beta_rho ~ dnorm(0, sd = 10)                                 # pop growth fixed effect prior

    # Observation process 
  chi_mu ~ dnorm(0, sd = 10)                                   # BBS observation: mean of effect
  chi_new ~ dnorm(0, sd = 10)                                  # BBS observation: new observer effect
  chi_route_sd ~ dgamma(.01, .01)                              # BBS observation: route-level sd
  chi_obser_sd ~ dgamma(.01, .01)                              # BBS observation: observer-level sd
  chi_stop_beta[1] ~ dnorm(0, sd = 10)                         # BBS observation: stop-level fixed effect
  chi_stop_beta[2] ~ T(dnorm(0, sd = 10), 0, )
  
  # Spatial process priors
  sigma2_s ~ dinvgamma(2, 1)                                   # spatial random effect sd
  phi_s ~ dunif(0, 500)                                        # spatial range of knot interpolation

  # Spatial structure: sparse Gaussian process: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Fix: Loop over knots to assign univariate priors
  for(k in 1:n_knots) {
    w[k] ~ dnorm(0, sd = sqrt(sigma2_s))                       # latent spatial random effect at each knot
  }

  # Interpolate to all sites using inverse distance weighting
  for(i in 1:nsite) {
    for(k in 1:n_knots) {
      weight[i,k] <- exp(-dist_all_knots[i,k] / phi_s)
    }
    weight_sum[i] <- sum(weight[i, 1:n_knots])
    s[i] <- inprod(weight[i, 1:n_knots], w[1:n_knots]) / weight_sum[i]
  }

  # Abundance model: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(i in 1:nsite) {                                          
    z[i] ~ dbinom(zeta, 1)                                     # true abundance 
    lambda0[i] <- exp(beta_lambda0 + s[i])                     # expected abundance (also fixed the indexing here!)
    N[i,1] ~ dpois(lambda0[i] * z[i] + 0.1 * (1 - z[i]))       # t=1: Changed N[i] to N[i,1]
  }
  
  # Compute rho once (it's constant across sites)
  rho_val <- exp(beta_rho)                                     # compute rho

  for(i in 1:nsite) { 
    for (t in 2:bbs_max_years) {                               # for t > 1
      N[i,t] ~ dpois(
        N[i, t-1] * rho_val * z[i] + 0.1 * (1 - z[i]))
    }
  }

  # Observation submodels: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # BBS
  for (j in 1:bbs_nroute) {
    chi_route_epsilon[j] ~ dnorm(0, sd = chi_route_sd) # route-level variation
  }

  for (j in 1:bbs_nobser) {
    chi_obser_epsilon[j] ~ dnorm(0, sd = chi_obser_sd) # observer-level variation
  }
  
  for (j in 1:bbs_nstop) {
    chi_stop_sd[j] <- exp(chi_stop_beta[1] + chi_stop_beta[2] * bbs_stops[j])
  }
  
  for (k in 1:bbs_nobs) {
    log_chi_mu[k] <- chi_mu + 
      chi_route_epsilon[bbs_route[k]] + 
      chi_obser_epsilon[bbs_obser[k]] + 
      chi_new * bbs_new[k]
    
    for (j in 1:bbs_nstop) {
      log_chi[k,j] ~ dnorm(log_chi_mu[k], sd = chi_stop_sd[j])
      chi[k,j] <- exp(log_chi[k,j])
      bbs_cnt[k,j] ~ dpois(N[bbs_sites[k], bbs_time_step[k]] * chi[k,j])
    }
  }

    # EBD observation model
    # BCR observation model 
})


# Create model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model <- nimbleModel(
  code, 
  constants = constants, 
  data = data, 
  inits = inits, 
  calculate = FALSE)
write_rds(model, "data/bbs/model/nmix_gpp_model_object.rds")
model <- read_rds("data/bbs/model/nmix_gpp_model_object.rds")


# Check initializeation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model$initializeInfo()


## Calculate initial log probability 

cat("Initial log probability:", model$calculate(), "\n")


# Configure MCMC and compile movel ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Configure MCMC
mcmc_conf <- configureMCMC(model)

## Change samplers
mcmc_conf$removeSamplers(c('beta_lambda0', 'beta_rho', 'chi_stop_beta', 'w'))
mcmc_conf$addSampler(target = c('beta_lambda0'), type = "RW_block")
mcmc_conf$addSampler(target = c('beta_rho'), type = "RW_block")
mcmc_conf$addSampler(target = c('chi_stop_beta'), type = "RW_block")
mcmc_conf$addSampler(target = c('w'), type = "RW_block")

## Add monitors
mcmc_conf$addMonitors(c("N")) 

mcmc <- buildMCMC(mcmc_conf)

compiled <- compileNimble(
  model, 
  mcmc, 
  showCompilerOutput = F, 
  resetFunctions = T
); write_rds(compiled, "data/bbs/model/nmix_gpp_compiled.rds")
compiled <- read_rds("data/bbs/model/nmix_gpp_compiled.rds")


# Fit model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fit <- runMCMC(
    compiled$mcmc
  , nchains = 1          # for later model runs, set to 3 to assess convergence 
  , niter = 2000         # enough iters to see if parameters are moving toward sensible values; set to 10k later
  , nburnin = 1000       # have of the inter for now; set to 5k later
  , summary = TRUE
  , samples = TRUE
  , progressBar = T
); write_rds(fit, "data/bbs/model/nmix_gpp_test.rds")



# Extract MCMC posterior samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N_summary <- fit$summary$all.chains  %>%
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  filter(str_detect(parameter, "^N\\[")) %>%
  mutate(
    # Extract the numbers between brackets
    site = as.integer(str_extract(parameter, "(?<=\\[)\\d+")),
    year = as.integer(str_extract(parameter, "(?<=, )\\d+(?=\\])"))
  ); head(N_summary, 10)

N_summary_strat <- N_summary %>%
    left_join(
        strat %>%
            mutate(site = row_number()) %>%
            select(site, cell_id, x),
        by = "site") %>%
  st_as_sf()

write_rds(N_summary_strat, "data/bbs/model/nmix_gpp_posterior_summary.rds")

N_summary_strat %>%
  filter(year == 10) %>%
  ggplot() +
  geom_sf(aes(fill = log(Mean + 1))) +
  scale_fill_viridis_c() +
  labs(title = "Predicted Abundance - Year 1") +
  theme_minimal()


