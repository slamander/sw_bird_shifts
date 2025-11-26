pacman::p_load(
  tidyverse
  , here
  , sf
  , nimble
  # , 
  # , 
)
# ============================================================================
# Hierarchical Spatiotemporal N-mixture Model with Data Integration
# ============================================================================

# NIMBLE model code
model_code <- nimbleCode({
  
  # =========================================================================
  # PRIORS
  # =========================================================================
  
  # --- Abundance model parameters ---
  # Intercept and covariate effects on abundance
  beta0 ~ dnorm(0, sd = 10)           # Abundance intercept
  beta_habitat ~ dnorm(0, sd = 10)    # Habitat covariate effect
  beta_temp ~ dnorm(0, sd = 10)       # Temperature effect (temporal)
  
  # Spatial random effects
  sigma_spatial ~ dunif(0, 10)        # SD of spatial random effects
  tau_spatial <- 1 / (sigma_spatial^2)
  
  for(s in 1:n_sites) {
    eta_spatial[s] ~ dnorm(0, tau = tau_spatial)
  }
  
  # Temporal random effects
  sigma_temporal ~ dunif(0, 10)       # SD of temporal random effects
  tau_temporal <- 1 / (sigma_temporal^2)
  
  for(t in 1:n_years) {
    eta_temporal[t] ~ dnorm(0, tau = tau_temporal)
  }
  
  # Spatiotemporal interaction (optional but recommended)
  sigma_st ~ dunif(0, 10)
  tau_st <- 1 / (sigma_st^2)
  
  for(s in 1:n_sites) {
    for(t in 1:n_years) {
      eta_st[s, t] ~ dnorm(0, tau = tau_st)
    }
  }
  
  # --- Detection model parameters ---
  
  # Data source 1: Standardized surveys (distance + removal)
  alpha0_std ~ dnorm(0, sd = 10)      # Detection intercept
  alpha_distance ~ dnorm(0, sd = 10)  # Distance decay effect
  alpha_removal ~ dnorm(0, sd = 10)   # Removal sampling effect
  alpha_obs1 ~ dnorm(0, sd = 10)      # Observer effect
  
  # Data source 2: Opportunistic data (observer error)
  alpha0_opp1 ~ dnorm(0, sd = 10)     # Detection intercept
  alpha_effort ~ dnorm(0, sd = 10)    # Effort/duration effect
  
  # Observer random effects for opportunistic data
  sigma_observer ~ dunif(0, 5)
  tau_observer <- 1 / (sigma_observer^2)
  
  for(o in 1:n_observers) {
    eta_observer[o] ~ dnorm(0, tau = tau_observer)
  }
  
  # Data source 3: Opportunistic data (nested spatial design)
  alpha0_opp2 ~ dnorm(0, sd = 10)     # Detection intercept
  alpha_spatial_scale ~ dnorm(0, sd = 10)  # Spatial scale effect
  
  # Nested spatial random effects (e.g., region within site)
  sigma_nested ~ dunif(0, 5)
  tau_nested <- 1 / (sigma_nested^2)
  
  for(r in 1:n_regions) {
    eta_nested[r] ~ dnorm(0, tau = tau_nested)
  }
  
  # =========================================================================
  # STATE PROCESS: Abundance model
  # =========================================================================
  
for(s in 1:n_sites) {
   for(t in 1:n_years) {
  # Log-linear model for expected abundance
    log(lambda[s, t]) <- exp(
      beta0 + 
      beta_habitat * habitat[s] +
      beta_temp * temperature[t] +
      eta_spatial[s] +
      eta_temporal[t] +
      eta_st[s, t]
    
    # True abundance (latent state)
    N[s, t] ~ dpois(lambda[s, t])
  }# t
}# s 
  
  # =========================================================================
  # OBSERVATION PROCESS 1: Standardized surveys
  # =========================================================================
  
  for(i in 1:n_obs_std) {
    
    # Detection probability with distance decay and removal
    logit(p_std[i]) <- alpha0_std + 
                       alpha_distance * distance[i] +
                       alpha_removal * removal_period[i] +
                       alpha_obs1 * observer_cov[i]
    
    # Observed counts (binomial observation model)
    y_std[i] ~ dbin(p_std[i], N[site_std[i], year_std[i]])
  }
  
  # =========================================================================
  # OBSERVATION PROCESS 2: Opportunistic data (observer error)
  # =========================================================================
  
  for(i in 1:n_obs_opp1) {
    
    # Detection probability with effort and observer effects
    logit(p_opp1[i]) <- alpha0_opp1 + 
                        alpha_effort * effort[i] +
                        eta_observer[observer_id[i]]
    
    # Observed counts
    y_opp1[i] ~ dbin(p_opp1[i], N[site_opp1[i], year_opp1[i]])
  }
  
  # =========================================================================
  # OBSERVATION PROCESS 3: Opportunistic data (nested spatial design)
  # =========================================================================
  
  for(i in 1:n_obs_opp2) {
    
    # Detection probability with nested spatial effects
    logit(p_opp2[i]) <- alpha0_opp2 + 
                        alpha_spatial_scale * spatial_scale[i] +
                        eta_nested[region_id[i]]
    
    # Observed counts
    y_opp2[i] ~ dbin(p_opp2[i], N[site_opp2[i], year_opp2[i]])
  }
  
  # =========================================================================
  # DERIVED QUANTITIES (optional)
  # =========================================================================
  
  # Total abundance across all sites by year
  for(t in 1:n_years) {
    total_abundance[t] <- sum(N[1:n_sites, t])
  }
  
  # Mean detection probability by data source
  mean_p_std <- mean(p_std[1:n_obs_std])
  mean_p_opp1 <- mean(p_opp1[1:n_obs_opp1])
  mean_p_opp2 <- mean(p_opp2[1:n_obs_opp2])
})

# ============================================================================
# DATA PREPARATION EXAMPLE
# ============================================================================

# This is a template - replace with your actual data
constants <- list(
  n_sites = 100,          # Number of spatial sites
  n_years = 10,           # Number of years
  n_observers = 25,       # Number of unique observers
  n_regions = 20,         # Number of nested regions
  
  # Sample sizes for each data source
  n_obs_std = 500,        # Standardized survey observations
  n_obs_opp1 = 1000,      # Opportunistic observations (observer error)
  n_obs_opp2 = 800        # Opportunistic observations (nested spatial)
)

# Data list (you'll need to populate these with your actual data)
data <- list(
  # Observed counts from each data source
  y_std = rep(NA, constants$n_obs_std),
  y_opp1 = rep(NA, constants$n_obs_opp1),
  y_opp2 = rep(NA, constants$n_obs_opp2),
  
  # Covariates for abundance model
  habitat = rnorm(constants$n_sites),
  temperature = rnorm(constants$n_years),
  
  # Covariates and indices for standardized surveys
  distance = runif(constants$n_obs_std, 0, 100),
  removal_period = sample(1:3, constants$n_obs_std, replace = TRUE),
  observer_cov = rnorm(constants$n_obs_std),
  site_std = sample(1:constants$n_sites, constants$n_obs_std, replace = TRUE),
  year_std = sample(1:constants$n_years, constants$n_obs_std, replace = TRUE),
  
  # Covariates and indices for opportunistic data 1
  effort = runif(constants$n_obs_opp1, 1, 10),
  observer_id = sample(1:constants$n_observers, constants$n_obs_opp1, replace = TRUE),
  site_opp1 = sample(1:constants$n_sites, constants$n_obs_opp1, replace = TRUE),
  year_opp1 = sample(1:constants$n_years, constants$n_obs_opp1, replace = TRUE),
  
  # Covariates and indices for opportunistic data 2
  spatial_scale = rnorm(constants$n_obs_opp2),
  region_id = sample(1:constants$n_regions, constants$n_obs_opp2, replace = TRUE),
  site_opp2 = sample(1:constants$n_sites, constants$n_obs_opp2, replace = TRUE),
  year_opp2 = sample(1:constants$n_years, constants$n_obs_opp2, replace = TRUE)
)

# ============================================================================
# INITIAL VALUES
# ============================================================================

# Function to generate initial values
inits <- function() {
  list(
    # Abundance parameters
    beta0 = rnorm(1, 0, 1),
    beta_habitat = rnorm(1, 0, 0_5),
    beta_temp = rnorm(1, 0, 0_5),
    
    # Detection parameters
    alpha0_std = rnorm(1, 0, 1),
    alpha_distance = rnorm(1, -0_5, 0_2),
    alpha_removal = rnorm(1, 0, 0_5),
    alpha_obs1 = rnorm(1, 0, 0_5),
    alpha0_opp1 = rnorm(1, 0, 1),
    alpha_effort = rnorm(1, 0_5, 0_2),
    alpha0_opp2 = rnorm(1, 0, 1),
    alpha_spatial_scale = rnorm(1, 0, 0_5),
    
    # Random effects variances
    sigma_spatial = runif(1, 0_1, 2),
    sigma_temporal = runif(1, 0_1, 2),
    sigma_st = runif(1, 0_1, 1),
    sigma_observer = runif(1, 0_1, 1),
    sigma_nested = runif(1, 0_1, 1),
    
    # Latent abundance (initialize with reasonable values)
    N = matrix(rpois(constants$n_sites * constants$n_years, 10), 
               nrow = constants$n_sites, ncol = constants$n_years)
  )
}

# ============================================================================
# MODEL COMPILATION AND MCMC SETUP
# ============================================================================

# Create the NIMBLE model
nimble_model <- nimbleModel(code = model_code,
                            constants = constants,
                            data = data,
                            inits = inits())

# Compile the model
C_model <- compileNimble(nimble_model)

# Configure MCMC
mcmc_conf <- configureMCMC(nimble_model, 
                          monitors = c('beta0', 'beta_habitat', 'beta_temp',
                                      'alpha0_std', 'alpha_distance', 'alpha_removal',
                                      'alpha0_opp1', 'alpha_effort',
                                      'alpha0_opp2', 'alpha_spatial_scale',
                                      'sigma_spatial', 'sigma_temporal', 'sigma_st',
                                      'sigma_observer', 'sigma_nested',
                                      'total_abundance',
                                      'mean_p_std', 'mean_p_opp1', 'mean_p_opp2'),
                          thin = 10,
                          print = TRUE)

# Build MCMC
mcmc_build <- buildMCMC(mcmc_conf)

# Compile MCMC
C_mcmc <- compileNimble(mcmc_build, project = nimble_model)

# ============================================================================
# RUN MCMC
# ============================================================================

# Run the MCMC (adjust niter, nburnin, nchains as needed)
samples <- runMCMC(
  C_mcmc,
  niter = 50000,
  nburnin = 10000,
  nchains = 3,
  thin = 10,
  samplesAsCodaMCMC = TRUE,
  summary = TRUE)

# ============================================================================
# EXAMINE RESULTS
# ============================================================================

# View summary statistics
print(samples$summary)

# Check convergence
library(coda)
gelman_diag(samples$samples)

# Trace plots
plot(samples$samples[, c('beta0', 'alpha0_std', 'sigma_spatial')])

# Posterior distributions
library(ggplot2)
library(bayesplot)
mcmc_dens(samples$samples, pars = c('beta0', 'beta_habitat', 'beta_temp'))
mcmc_dens(samples$samples, pars = c('alpha0_std', 'alpha0_opp1', 'alpha0_opp2'))
