# JAGS model run for snow fences example

library(rjags)
library(mcmcplots)
library(broom.mixed)
library(tidyverse)
load.module('dic')

mscript <- "model{
  for(i in 1:N){
    # Likelihood
    y[i] ~ dpois(theta[i]*x[i])
    
    # Prior for theta
    theta[i] ~ dgamma(alpha, beta)
  }
  
  # Root node priors for alpha and beta
  alpha ~ dgamma(2, 1)
  beta ~ dexp(1)
  
  # Calculated popn mean and std
  pmean <- alpha/beta
  pstd <- sqrt(alpha/pow(beta, 2))
}"

# Data
sf <- data.frame(y = c(138, 91, 132, 123, 173, 124, 109, 154, 138, 134),
                 x = c(72, 50, 55, 60, 78, 63, 54, 70, 80, 68))

datlist <- list(y = sf$y,
                x = sf$x,
                N = nrow(sf))

# Initials
mean(sf$y/sf$x) # sample mean

# target value of 0.5, 5, 20
initslist <- list(list(alpha = 0.5, beta = 1),
                  list(alpha = 5, beta = 1),
                  list(alpha = 40, beta = 2))

# Model
jm <- jags.model(file = textConnection(mscript),
                 data = datlist,
                 inits = initslist,
                 n.chains = 3)
update(jm, 1000)

jm_coda <- coda.samples(model = jm,
                        variable.names = c("theta", 
                                           "alpha", "beta",
                                           "pmean", "pstd"),
                        n.iter = 3000)

# model diagnostics
traplot(jm_coda,
        parms = c("theta", "alpha", "beta"))

mcmcplot(jm_coda, 
         parms = c("alpha", "beta"))

# Re-run with thinning by 20
jm <- jags.model(file = textConnection(mscript),
                 data = datlist,
                 inits = initslist,
                 n.chains = 3)
jm_coda <- coda.samples(model = jm,
                        variable.names = c("theta", 
                                           "alpha", "beta",
                                           "pmean", "pstd"),
                        n.iter = 3000*20,
                        thin = 20)
mcmcplot(jm_coda,
         parms = c("alpha", "beta"))
nrow(jm_coda[[1]])

# Rhat, or Gelman diagnostic
gelman.diag(jm_coda, multivariate = FALSE)

# catepillar plots
caterplot(jm_coda, parms = "theta",
          reorder = FALSE)

caterplot(jm_coda, parms = "pmean",
          reorder = FALSE)

caterplot(jm_coda, parms = "pstd",
          reorder = FALSE)

# Summarizing posterior
summary(jm_coda)[[1]]

head(jm_coda[[1]])

tidyMCMC(jm_coda,
         conf.int = TRUE) %>%
  filter(grepl("theta", term)) %>%
  ggplot(aes(x = term, y = estimate)) +
  geom_pointrange(aes(ymin = conf.low,
                      ymax = conf.high))

# effective number of parameters
dic.samples(jm, n.iter = 3000)
