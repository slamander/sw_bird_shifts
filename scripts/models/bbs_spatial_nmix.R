rm(list = ls())
pacman::p_load(
  tidyverse
  , here
  , sf
  , spdep 
  , nimble
  # , 
  # , 
)

region <- read_rds("data/region/region.rds")

# load gridded count data
bbs <- read_rds("data/bbs/species/strat/pija.rds") %>%
  mutate(year = as.numeric(year)) %>%
  left_join(
    data.frame(
      route_id = unique(.$route_id), 
      site = seq(1:length(unique(.$route_id))))) %>%
  arrange(route_id, year) %>%
  group_by(route_id) %>%
  mutate(
    time_step = row_number(),
    nyears_site = n()) %>%
  ungroup() %>%
  mutate(
    obs_id = as.numeric(factor(obs_n, levels = unique(obs_n)))
  )

ggplot() +
  # geom_sf(data = region, fill = NA) +
  geom_sf(data = bbs, aes(fill = total_count)) +
  scale_fill_continuous(na.value = "transparent") +
  theme_void()
ggsave("figures/bbs_strat_pija.jpeg")
  
# reduce to only unique sites
bbs_sites <- bbs %>%
  group_by(cell_id) %>%
  summarize(total_count = sum(total_count))

cat(
  "Number of grid sites with observations", 
  sum(is.na(bbs_sites$total_count) == F), 
  "\n")

cat(
  "Number of grid sites without observations", 
  sum(is.na(bbs_sites$total_count)), 
  "\n")

# create adjacency matrix
neighbors <- bbs_sites %>%
  st_centroid() %>%
  st_coordinates() %>%
  knearneigh(k = 1) %>%
  knn2nb() %>%
  make.sym.nb()

weights <- nb2listw(neighbors, style = "B", zero.policy = TRUE)
car_adjacency <- as.carAdjacency(neighbors, weights$weights)

# Observed data
data <- list(
  bbs_cnt = bbs %>%
    select(count1:count10) %>% 
    st_drop_geometry() %>%
    as.matrix()
)

# Constants
constants <- list(
    nsite = nrow(bbs_sites)
  , max_years = max(bbs$nyears_site)
  , nobs_bbs = nrow(bbs)
  , sites_bbs = bbs$site
  , years_site = bbs$time_step
  , nroute_bbs = length(unique(bbs$site))
  , route_bbs = bbs$site
  , nobser_bbs = length(unique(bbs$obs_id))
  , obser_bbs = bbs$obs_id
  , nstop_bbs = ncol(data$bbs_cnt)
  , stops_bbs = seq(-2, 2, length.out = ncol(data$bbs_cnt))
  , new_bbs = bbs$new

# CAR constants
  , adj = car_adjacency$adj         # Vector of neighbor indices
  , weights = car_adjacency$weights # Vector of weights corresponding to adj
  , num = car_adjacency$num         # Number of neighbors for each site
  , L = length(car_adjacency$adj)   # Total number of neighbor pairs
)

# Model code
code <- nimbleCode({ 

  # Priors
  zeta ~ dunif(0, 1)                            # true counts
  beta_lambda0 ~ dnorm(0, sd = 10)              # state process intercept
  beta_rho ~ dnorm(0, sd = 10)                  # growth rate intercept
  chi_mu ~ dnorm(0, sd = 10)                    # observation process mean
  chi_route_sd ~ dgamma(.01, .01)               # observation process route-level sd
  chi_obser_sd ~ dgamma(.01, .01)               # observation process observer-level sd
  chi_new ~ dnorm(0, sd = 10)                   # observation process new observer
  chi_stop_beta[1] ~ dnorm(0, sd = 10)          # observation process route replicate 1
  chi_stop_beta[2] ~ T(dnorm(0, sd = 10), 0, )  # observation process route replicate 2
  tau ~ dgamma(0.001, 0.001)                    # CAR precision prior

  # CAR process
  s[1:nsite] ~ dcar_normal(
    adj[1:L], 
    weights[1:L], 
    num[1:nsite], 
    tau, 
    zero_mean = 1)

  # Process model
    # for t = 1
  for (i in 1:nsite) {
    z[i] ~ dbinom(zeta, 1)
    lambda0[i] <- exp(beta_lambda0)
    N[i,1] ~ dpois(lambda0[i] * z[i] + 0.1 * (1 - z[i]) + s[i]) 
    
    # for t > 1
    for (t in 2:max_years) {
      rho[i, t-1] <- exp(beta_rho)
      N[i,t] ~ dpois(N[i, t-1] * rho[i,t-1] * z[i] + 0.1 * (1 - z[i]))
    }
  }

  # Observation model
  for (k in 1:nroute_bbs) {
    chi_route_epsilon[k] ~ dnorm(0, sd = chi_route_sd)
  }
    
  for (k in 1:nobser_bbs) {
    chi_obser_epsilon[k] ~ dnorm(0, sd = chi_obser_sd)
  }
    
  for (j in 1:nstop_bbs) {
    chi_stop_sd[j] <- exp(
      chi_stop_beta[1] + chi_stop_beta[2] * stops_bbs[j])
  }
    
  for (k in 1:nobs_bbs) {
    log_chi_mu[k] <- 
      chi_mu + 
      chi_route_epsilon[route_bbs[k]] + 
      chi_obser_epsilon[obser_bbs[k]] + 
      chi_new * new_bbs[k]
    
    for (j in 1:nstop_bbs) {
      log_chi[k,j] ~ dnorm(log_chi_mu[k], sd = chi_stop_sd[j])
      chi[k,j] <- exp(log_chi[k,j])
      bbs_cnt[k,j] ~ dpois(N[sites_bbs[k], years_site[k]] * chi[k,j])
    }
  }
})

# Initialize N matrix - same as before
Ni <- matrix(10, constants$nsite, constants$max_years)
for(i in 1:constants$nsite){
  site_data <- bbs[bbs$site == i, ]
  for(j in 1:nrow(site_data)){
    time_step <- site_data$time_step[j]
    obs_counts <- data$bbs_cnt[which(bbs$site == i & bbs$time_step == time_step), ]
    if(length(obs_counts) > 0){
      Ni[i, time_step] <- max(max(obs_counts, na.rm = TRUE), 10) * 3
    }
  }
}

zi <- ifelse(rowSums(Ni) == 10 * constants$max_years, 0, 1)

inits <- list(
    zeta = 0.5
  , beta_lambda0 = 0.5
  , beta_rho = -0.05
  , chi_mu = -5
  , chi_route_sd = 0.5
  , chi_obser_sd = 0.5
  , chi_new = 0
  , chi_stop_beta = c(-1, 0.1)
  , N = Ni
  , z = zi
  , chi_route_epsilon = rep(0, constants$nroute_bbs)
  , chi_obser_epsilon = rep(0, constants$nobser_bbs)
  , log_chi = matrix(-5, constants$nobs_bbs, constants$nstop_bbs)
  , tau = 1
  , s = rep(0, nrow(bbs_sites))
)

# Create model
model <- nimbleModel(
  code, 
  constants = constants, 
  data = data, 
  inits = inits)

# Check what's not initialized:
model$initializeInfo()

# Calculate initial log probability
cat("Initial log probability:", model$calculate(), "\n")

# Configure MCMC
mcmc_conf <- configureMCMC(model)

mcmc_conf$removeSamplers(c('beta_lambda0', 'beta_rho', 'chi_stop_beta'))

mcmc_conf$addSampler(target = c('beta_lambda0'), type = "RW_block")
mcmc_conf$addSampler(target = c('beta_rho'), type = "RW_block")
mcmc_conf$addSampler(target = c('chi_stop_beta'), type = "RW_block")

mcmc_conf$addMonitors(c('N', 'z', 's')) 

mcmc <- buildMCMC(mcmc_conf)

compiled <- compileNimble(
  model, 
  mcmc, 
  showCompilerOutput = TRUE, 
  resetFunctions = TRUE
)

fit <- runMCMC(
    compiled$mcmc
  , nchains = 3
  , niter = 10000
  , nburnin = 5000
  , thin = 10
  , summary = TRUE
  , samples = TRUE
  , progressBar = T
)

write_rds(fit, "data/bbs/model/nmix.rds")

fit$summary

data_nodes <- model$getNodeNames(dataOnly = TRUE)
parent_nodes <- model$getParents(data_nodes, stochOnly = TRUE)  # `getParents` is new in nimble 0.11.0
sim_nodes <- model$getDependencies(parent_nodes, self = FALSE)

model$simulate(sim_nodes, includeData = TRUE)

model$N
