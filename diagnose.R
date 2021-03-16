library(tidyverse)
library(cmdstanr)
library(bayesplot)

color_scheme_set("darkgray")

dist_fit <- read_rds("data/stan_analysis_data/dist_fit.rds")

# dist_fit <- dir("data/stan_analysis_data/", "dist_fit37.*csv", full.names = TRUE) %>% 
#   read_cmdstan_csv(variables = c("raw_u_sd", "structural_cluster_obs_v", "wtp_sigma", "wtp_value_utility"))

# post <- dist_fit$post_warmup_draws
post <- dist_fit$draws(c("group_dist_mean", "group_dist_sd"))
np <- nuts_params(dist_fit)

mcmc_parcoord(post, np = np)

# c("raw_u_sd[1]", "raw_u_sd[2]", "raw_u_sd[3]",  "wtp_sigma", "wtp_value_utility") %>% 
c("raw_u_sd[1]", "raw_u_sd[2]", "raw_u_sd[3]", "structural_cluster_obs_v[100]",  "wtp_sigma", "wtp_value_utility") %>% 
  mcmc_pairs(post[,, .], np = np, pars = ., off_diag_args = list(size = 0.75))

mcmc_pairs(post, np = np, off_diag_args = list(size = 0.75))
