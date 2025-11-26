pacman::p_load(
  tidyverse 
  , here
  , nimble
  , rstan
  , auk
  , sf
  # ,
  # , 
)

obs_file <- here("data/ebd/ebd/observations_clean.txt")
chk_file <- here("data/ebd/ebd/checklists_clean.txt")

crs <- read_rds("data/region/crs.rds")

# grsg_raw <- obs_file %>%
#   auk_ebd() %>%
#   auk_species(species = "Greater Sage-Grouse") %>%
#   auk_bbox(bbox = c(-128.71584, 24.26827, -98.81304, 53.97169)) %>%
#   auk_filter(
#     file = "data/ebd/filters/grsg.txt", 
#     overwrite = TRUE) %>%
#   read_ebd() %>%
#   st_as_sf(
#     coords = c("longitude", "latitude"), 
#     crs = st_crs(4326)) %>%
#   st_transform(crs) %>%
#   select(
#     checklist_id, taxonomic_order, category, common_name, scientific_name, 
#     observation_count, breeding_code, breeding_category, behavior_code, country, 
#     state, county, bcr_code, locality_id, locality_type, observation_date, 
#     time_observations_started, observer_id, observation_type, effort_area_ha, 
#     group_identifier, geometry)
# write_rds(grsg_raw, "data/ebd/species/point/examples/grsg.rds")
grsg <- read_rds("data/ebd/species/point/examples/grsg.rds")

grsg_sum <- grsg %>%
  group_by(state) %>%
  summarize(
    N_ids = n_distinct(locality_id), 
    N_rep_ids = sum(duplicated(locality_id)), 

    mean_rep_ids = mean(table(locality_id)),
    med_rep_ids = median(table(locality_id)),
    max_rep_ids = max(table(locality_id))
  ); grsg_sum

plot(density(table(grsg$locality_id)))

ebird_sum <- ebird %>%
  summarize(
    N_loc = n_distinct(Location), 
    N_fish_id = n_distinct(fish_id), 
    N_FishID = n_distinct(FishID), 

    N_rep_loc = sum(duplicated(Location)), 
    N_rep_fish_id = sum(duplicated(fish_id)), 
    N_rep_FishID = sum(duplicated(FishID)), 

    mean_rep_loc = sum(duplicated(Location)/N_loc),
    mean_fish_id = sum(duplicated(fish_id)/N_fish_id),
    mean_FishID = sum(duplicated(FishID)/N_FishID)
  ); ebird_sum

fcode <- nimbleCode({

  # Priors
  ## binomial presences parameter
  zeta ~ dunif(0, 1)

  ## abundance covariates
  for (k in 1:4) {
    beta_lambda0[k] ~ dnorm(0, sd = 10)
    beta_rho[k] ~ dnorm(0, sd = 10)
  } # k
  
  ## ebird observation parameters
  omega_locat_sd ~ dgamma(.01, .01)
  omega_obser_sd ~ dgamma(.01, .01)
  omega_zeta ~ dunif(0, 1)
  omega_delta ~ dgamma(.01, .01)

  # Process model
  for (i in 1:nsite) {

    z[i] ~ dbinom(zeta, 1)

    lambda0[i] <- exp(
      beta_lambda0[1] + 
      beta_lambda0[2] * gra[i,1] + 
      beta_lambda0[3] * tem[i,1] + 
      beta_lambda0[4] * tem[i,1] * tem[i,1])

    N[i,1] ~ dpois(lambda0[i] * z[i] + 5e-3 * (1 - z[i])) # 5e-3 -- scaling parameter? 

    for (t in 2:nyear) {
      rho[i,t-1] <- exp(
        beta_rho[1] + 
        beta_rho[2] * gra[i,t] + 
        beta_rho[3] * tem[i,t] + 
        beta_rho[4] * tem[i,t] * tem[i,t])

      N[i,t] ~ dpois(N[i,t-1] * rho[i,t-1] * z[i] + 5e-3 * (1 - z[i]))

    } # t
  } # i

  # Observation model, eBird
  for (k in 1:nlocat_ebird) {
    omega_locat_epsilon[k] ~ dnorm(0, sd=omega_locat_sd)
  } # k

  for (k in 1:nobser_ebird) {
    omega_obser_epsilon[k] ~ dnorm(0, sd=omega_obser_sd)
  } # k

  for (k in 1:nobs_ebird) {

    omega_z[k] ~ dbinom(omega_zeta, 1)

    omega[k] <- exp(
      omega_beta[1] + 
      omega_beta[2] * type_ebird[k] + 
      omega_beta[3] * duration_ebird[k] + 
      omega_beta[4] * distance_ebird[k] + 
      omega_beta[5] * n_obsver_ebird[k] + 
      omega_locat_epsilon[locat_ebird[k]] + 
      omega_obser_epsilon[obser_ebird[k]])

    ebird_cnt[k] ~ dpois(N[sites_ebird[k], 
      years_ebird[k]] * omega[k] * omega_z[k] + 
        omega_delta * (1 - omega_z[k])
    )

    } # j
  } # k
) # nimbleCode
