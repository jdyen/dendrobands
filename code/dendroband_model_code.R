# R code to fit Bayesian hierarchical model to tree growth data

# optional: set working directory
# setwd("PATH/TO/DIR")

# load packages
if (!require(rstan)) {
  install.packages("rstan")
  library(rstan)
}

# load data
data_set <- get(load("./data/dendroband_data.R"))

# compile stan model
stan_mod <- stan_model(file = "./code/dendro_mod.stan")

# sample from stan model
mod <- sampling(object = stan_mod,
                data = data_set,
                iter = 10000,
                chains = 4,
                thin = 2,
                init = 'random',
                control = list(adapt_delta = 0.9,
                               max_treedepth = 30),
                cores = 4)

# save fitted model summary and data sets
out <- summary(mod, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975))$summary