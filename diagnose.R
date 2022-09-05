library(tidyverse)
library(cmdstanr)
library(bayesplot)

color_scheme_set("blue")

dist_prior <- dir("data/stan_analysis_data/", pattern = r"{dist_prior57_STRUCTURAL_LINEAR_U_SHOCKS-\d\.csv}", full.names = TRUE) %>% 
  as_cmdstan_fit()

dist_fit <- dir("data/stan_analysis_data/", pattern = r"{dist_fit57_STRUCTURAL_LINEAR_U_SHOCKS-\d\.csv}", full.names = TRUE) %>% 
  as_cmdstan_fit()

post <- dist_fit %>% 
  posterior::as_draws_array()

np <- dist_fit %>% 
  nuts_params()

# mcmc_parcoord(post) #, np = np)

# c("raw_u_sd[1]", "raw_u_sd[2]", "raw_u_sd[3]",  "wtp_sigma", "wtp_value_utility") %>% 
# c("raw_u_sd[1]", "raw_u_sd[2]", "raw_u_sd[3]", "structural_cluster_obs_v[100]",  "wtp_sigma", "wtp_value_utility") %>% 
# c("raw_u_sd[2]", "raw_u_sd[3]", "wtp_sigma", "beta_calendar_effect", "base_mu_rep", "wtp_value_utility") %>% 
# c("raw_u_sd[1]", "raw_u_sd[2]", "raw_u_sd[3]", "base_mu_rep", "wtp_value_utility") %>% 
c("raw_u_sd[1]", "base_mu_rep", "wtp_value_utility") %>% 
  mcmc_pairs(
    post[,, .], 
    np = np, 
    pars = ., 
    off_diag_args = list(size = 0.75), 
    np_style = pairs_style_np(div_size = 2),
    #transformations = list("raw_u_sd[1]" = "log")
    transformations = "log"
  )

# new_post[1:1000,,] %>% mcmc_parcoord()

list(prior = dist_prior, post = dist_fit) %>% 
  map(posterior::as_draws_df) %>% 
  map_dfr(tidybayes::spread_draws, wtp_value_utility, .id = "fit_type") %>% 
  ggplot() +
  tidybayes::geom_halfeyeh(aes(wtp_value_utility, fill = fit_type), alpha = 0.25)

list(prior = dist_prior, post = dist_fit) %>% 
  map(posterior::as_draws_df) %>% 
  map_dfr(tidybayes::spread_draws, raw_u_sd[treatment], .id = "fit_type") %>% 
  ggplot() +
  tidybayes::geom_halfeyeh(aes(raw_u_sd, fill = fit_type), alpha = 0.25)
