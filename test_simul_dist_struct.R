#!/usr/bin/Rscript

if (!interactive()) withr::with_dir("..", source(file.path("packrat", "init.R")))

library(magrittr)
library(tidyverse)
library(parallel)
library(pbmcapply)
library(HRW)
library(splines2)
library(nleqslv)
library(rstan)

options(mc.cores = 12)
rstan_options(auto_write = TRUE)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

num_clusters <- 140
cluster_size <- 50
dist_range <- c(0, 2500) # meters

# beta <- c(-0.05, 0, 0.15, 0.05)
beta <- rep(-0.05, 4)
beta_dist <- c(0.5, 0.0)
# mu <- c(0.05, 0.1, 0.12, 0.15)
mu <- rep(0.5, 4)

sigma_benefit <- 0.0
sigma_rep <- 0.0
sigma_v <- 1
sigma_cluster <- 0.0

rep_normal <- function(v, ...) dnorm(v, ...) / (pnorm(v, ...) * pnorm(v, ..., lower.tail = FALSE))
generate_v_cutoff_fixedpoint <- function(b, mu) {
  function(v_cutoff) {
    v_cutoff + b + mu * rep_normal(v_cutoff)
  }
}

sim_data <- tibble(
  cluster_id = seq(num_clusters),
  assigned_treatment = rep(as_factor(c("control", "ink", "calendar", "bracelet")), each = num_clusters %/% 4) %>% 
    sample(), 
  assigned_dist = runif(num_clusters, dist_range[1], dist_range[2])
) %>% 
  mutate(
    cluster_standard_dist = standardize(assigned_dist),
    cluster_dist_cost = beta_dist[1] * cluster_standard_dist + beta_dist[2] * cluster_standard_dist^2, 
    cluster_noise = rnorm(num_clusters, 0, sigma_cluster),
    cluster_net_benefit = beta[as.integer(assigned_treatment)] - cluster_dist_cost + cluster_noise,
    cluster_v_cutoff_fixedpoint = map2(cluster_net_benefit, mu[as.integer(assigned_treatment)], generate_v_cutoff_fixedpoint),
    cluster_v_cutoff = map2(cluster_v_cutoff_fixedpoint, cluster_net_benefit, 
                            ~ nleqslv(x = - ..2, fn = ..1)) %>% 
      map_dbl(pluck, "x"),
    cluster_delta_v = rep_normal(cluster_v_cutoff),
    cluster_rep = mu[as.integer(assigned_treatment)] * cluster_delta_v,
    hh_obs = rerun(num_clusters, tibble(hh_id = seq(cluster_size),
                                        benefit_noise = rnorm(cluster_size, 0, sigma_benefit),
                                        rep_noise = rnorm(cluster_size, 0, sigma_rep),
                                        v = rnorm(cluster_size, 0, sigma_v))) %>% 
      lst(hh_obs = ., cluster_net_benefit, cluster_rep, cluster_v_cutoff) %>% 
      pmap(function(hh_obs, cluster_net_benefit, cluster_rep, cluster_v_cutoff) {
        hh_obs %>% 
          mutate(y = cluster_net_benefit + cluster_rep + benefit_noise + rep_noise + v > 0)
      })
  )

sim_data %>% 
  ggplot() +
  # geom_line(aes(cluster_v_cutoff, cluster_delta_v)) +
  # geom_line(aes(assigned_dist, cluster_delta_v)) +
  geom_line(aes(assigned_dist, cluster_v_cutoff)) +
  # geom_line(aes(cluster_v_cutoff, cluster_rep, color = assigned_treatment)) +
  NULL

num_interior_knots <- 100

Z_osullivan <- calculate_splines(sim_data$cluster_standard_dist, num_interior_knots = num_interior_knots, spline_type = "osullivan")

grid_dist <- get_spline_range(sim_data$cluster_standard_dist) %>% unname() %>% list_modify(length = 1001) %>% do.call(seq, .)
       
stan_data <- lst(
  num_obs = map_int(sim_data$hh_obs, nrow) %>% sum(),
  num_grid_obs = length(grid_dist),
  num_treatments = 4,
  num_clusters,
  obs_cluster_id = unnest(sim_data, hh_obs) %>% pull(cluster_id),
  cluster_assigned_treatment = sim_data$assigned_treatment,
  takeup = unnest(sim_data, hh_obs) %>% pull(y),
  cluster_standard_dist = sim_data$cluster_standard_dist,
  
  Z_splines_v = Z_osullivan,
  num_knots_v = ncol(Z_splines_v),
  
  u_splines_v_sigma_sd = 1,
  mu_rep_sd = 0.05,
 
  num_excluded_clusters = 0,
  excluded_clusters = array(dim = 0),
 
  use_binomial = FALSE,
  use_cluster_effects = FALSE,
  use_mu_cluster_effects = FALSE,
  suppress_shocks = TRUE,
  use_cost_k_restrictions = FALSE,
  use_cost_model = cost_model_types["param_linear"], 
  use_private_incentive_restrictions = FALSE,
  use_salience_effect = FALSE,
  suppress_reputation = FALSE,
  simulate_new_data = FALSE,
  num_v_mix = 1,
) %>%  
  list_modify(!!!enum2stan_data(cost_model_types)) %>% 
  map_if(is.factor, as.integer)

dist_model <- stan_model(file.path("stan_models", "dist_struct_fixedpoint.stan"))
      
test_sim_fit <- sampling(
  dist_model,
  iter = 4000,
  thin = 1,
  chains = 4,
  control = lst(max_treedepth = 12, adapt_delta = 0.9),
  save_warmup = FALSE,
  refresh = 100,
  data = stan_data,
  init = generate_initializer(4, 
                              num_clusters,
                              structural_type = 1, 
                              num_mix = 1, 
                              use_cluster_effects = stan_data$use_cluster_effects,
                              use_mu_cluster_effects = stan_data$use_mu_cluster_effects,
                              restricted_private_incentive = FALSE,
                              semiparam = stan_data$use_cost_model == cost_model_types["semiparam"],
                              suppress_reputation = stan_data$suppress_reputation))

save(test_sim_fit, sim_data, stan_data, file = file.path("stan_analysis_data", "test_sim_dist.RData"))