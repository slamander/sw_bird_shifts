library(boot) 

#==============
# Basic values
#==============
nsite <- 60   # number of sites
nyear <- 20   # number of years
nreps <- 3    # number of replicates
ncovs <- 1    # number of environmental covariates
npcvs <- 2    # number of observational covariates
nspec <- 2    # number of species

# intercept and slopes for initial population size
beta0 <- matrix(c(5, 50, .5, .5), nspec, ncovs+1)

# standard deviation of stochasticity in initial population size
sigma0 <- rep(.1, nspec)

# intercept and slopes for local population growth rate
beta_phi <- matrix(c(0, 0, -.5, -.5, .5, -.5, 0, .5), nspec, nspec+ncovs+1)

# decay parameter for movement
kappa <- c(2, 2)

# standard deviation of stochasticity in population growth
sigma <- rep(.1, nspec)

# intercept and slopes for detection probability
beta_pobs <- matrix(c(
  logit(.67), 
  logit(.67), 
  .5, 
  .5), 
  nspec, 
  npcvs+1)

# standard deviation of stochasticity in detection probability
sigma_pobs <- rep(.1, nspec)

#===============
# Simulate data
#===============
### Site locations
lat <- runif(nsite, 0, 5)
lon <- runif(nsite, 0, 5)

### Distance between sites
dist <- matrix(, nsite, nsite)
for (i in 1:nsite) {
  for (j in 1:nsite) {
    dist[i,j] <- sqrt((lat[i] - lat[j]) ^ 2 + (lon[i] - lon[j]) ^ 2)
  } # j
} # i

### Environmental covariates
x_mean <- matrix(rnorm(nsite * ncovs, 0, 1), nsite, ncovs)
x <- array(, dim=c(nsite, nyear, ncovs))
for (i in 1:nsite) {
  for (s in 1:ncovs) {
    x[i,,s] <- rnorm(nyear, x_mean[i,s], .2)
  } # s
} # i

### Metapopulation dynamics
lambda0 <- matrix(, nsite, nspec)
for (s in 1:nspec) {
  lambda0[,s] <- exp(
    cbind(1, x[,1,]) %*% beta0[s,] + rnorm(nsite, 0, sigma0[s])
  )
} # s

N <- array(, dim=c(nsite, nyear, nspec))
for (s in 1:nspec) {
  N[,1,s] <- rpois(nsite, lambda0[,s])
} # s

eta <- theta <- array(, dim=c(nsite, nsite, nspec))
for (s in 1:nspec) {
  eta[,,s] <- exp(-1 * kappa[s] * dist)
  theta[,,s] <- eta[,,s] / rowSums(eta[,,s])
} # s

phi <- lambda <- array(, dim=c(nsite, nyear-1, nspec))
for (t in 2:nyear) {
  for (s in 1:nspec) {
    phi[,t-1,s] <- exp(
      cbind(1, (N[,t-1,]- lambda0) / lambda0, x[,t,] - x[,1,]) %*% 
      beta_phi[s,])
    lambda[,t-1,s] <- ((N[,t-1,s] * phi[,t-1,s]) %*% 
      theta[,,s]) * exp(rnorm(nsite, 0, sigma[s]))
    N[,t,s] <- rpois(nsite, lambda[,t-1,s])
  } # s
} # t

### Count data
w <- array(
  rnorm(nsite * nyear * npcvs * nreps, 0, 1), 
  dim=c(nsite, nyear, npcvs, nreps))
pobs <- y <- array(0, dim=c(nsite, nyear, nspec, nreps))
for (i in 1:nsite) {
  for (t in 1:nyear) {
    for (s in 1:nspec) {
      pobs[i,t,s,] <- inv.logit(cbind(1, t(w[i,t,,])) %*% 
        beta_pobs[s,] + rnorm(nreps, 0, sigma_pobs[s]))
      y[i,t,s,] <- rbinom(nreps, N[i,t,s], pobs[i,t,s,])
    } # s
  } # t
} # i