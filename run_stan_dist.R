#!/usr/bin/Rscript

if (!interactive()) withr::with_dir("..", source(file.path("packrat", "init.R")))

script_options <- docopt::docopt(
"Usage:
  run_stan_dist fit [--no-save --sequential --chains=<chains> --iter=<iter> --thin=<thin> --force-iter --models=<models> --outputname=<output file name> --update-output --no-cluster-effects --no-dist-cluster-effects --simulate-new]
  run_stan_dist cv [--folds=<number of folds> --no-save --sequential --chains=<chains> --iter=<iter> --thin=<thin> --force-iter --models=<models> --outputname=<output file name> --update-output --no-cluster-effects --no-dist-cluster-effects]
  
Options:
  --folds=<number of folds>  Cross validation folds [default: 10]
  --chains=<chains>  Number of Stan chains [default: 4]
  --iter=<iter>  Number of (warmup + sampling) iterations [default: 8000]
  --thin=<thin>  Thin samples [default: 1]",

  args = if (interactive()) "fit --sequential --models=9" else commandArgs(trailingOnly = TRUE) 
) 

library(magrittr)
library(tidyverse)
library(parallel)
library(pbmcapply)
library(HRW)
library(splines2)
library(rstan)
library(loo)

options(mc.cores = 12)
rstan_options(auto_write = TRUE)

folds <- as.integer(script_options$folds %||% 10)
chains <- as.integer(script_options$chains)
iter <- as.integer(script_options$iter)
output_name <- if (!is_null(script_options$outputname)) { script_options$outputname } else if (script_options$fit) { "dist_fit" } else { "dist_kfold" }
output_file_name <- file.path("stan_analysis_data", str_c(output_name, ".RData"))
use_cluster_effects <- !(script_options$`no-cluster-effects` %||% FALSE)
use_dist_cluster_effects <- !(script_options$`no-dist-cluster-effects` %||% FALSE)
models_to_run <- if (!is_null(script_options$models)) c(str_split(script_options$models, ",", simplify = TRUE)) %>% as.integer()
thin_by <- as.integer(script_options$thin)
model_file <- file.path("stan_models", script_options$`model-file`)
simulate_new_data <- as.logical(script_options$`simulate-new`)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

# Splines -----------------------------------------------------------------

num_interior_knots <- 100

cluster_standard_dist <- distinct(monitored_nosms_data, cluster_id, standard_cluster.dist.to.pot) %>% 
  arrange(cluster_id) %>% 
  pull(standard_cluster.dist.to.pot)

Z_osullivan <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, spline_type = "osullivan")
Z_i_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, spline_type = "i-spline")
Z_b_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, spline_type = "b-spline")

grid_dist <- get_spline_range(cluster_standard_dist) %>% unname() %>% list_modify(length = 1001) %>% do.call(seq, .)

Z_grid_osullivan <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, splines_for = grid_dist, spline_type = "osullivan")
Z_grid_i_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, splines_for = grid_dist, spline_type = "i-spline")
Z_grid_b_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, splines_for = grid_dist, spline_type = "b-spline")

treatment_dist <- distinct(monitored_nosms_data, assigned.treatment, standard_cluster.dist.to.pot) %>% 
  group_by(assigned.treatment) %>% 
  group_map(~ pull(.x, "standard_cluster.dist.to.pot"))

grid_dist2 <- treatment_dist %>% 
  map(get_spline_range) %>% 
  map(unname) %>% 
  map(list_modify, length = 1001) %>% 
  map(~ do.call(seq, .))

Z_grid_osullivan2 <- map2(treatment_dist, grid_dist2, ~ calculate_splines(.x, num_interior_knots = num_interior_knots, splines_for = .y, spline_type = "osullivan")) 

# Models ------------------------------------------------------------------

num_treatments <- n_distinct(monitored_nosms_data$assigned.treatment)
num_clusters <- n_distinct(monitored_nosms_data$cluster_id)

generate_initializer <- function(base_init = function() lst(beta_cluster_sd = abs(rnorm(num_treatments))), 
                                 structural_type = 0,
                                 num_mix = 1,
                                 restricted_private_incentive = FALSE,
                                 semiparam = FALSE,
                                 suppress_reputation = FALSE) {
  base_list <- base_init()
  
  if (structural_type > 0 || !is_empty(base_list)) {
    function() {
      if (structural_type > 0) {
        num_beta_param <- if (restricted_private_incentive) num_treatments - 1 else num_treatments
        # num_beta_param <- num_treatments
        
        base_list %>% 
          list_modify(
            mu_rep_raw = if (suppress_reputation) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
            mu_rep = if (suppress_reputation) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
            # structural_beta = rnorm(num_treatments),
            dist_cost_k_raw = if (semiparam) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
            dist_cost_k = if (semiparam) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
            # structural_beta_cluster = if (use_cluster_effects) matrix(rnorm(num_clusters * num_beta_param), nrow = num_clusters, ncol = num_beta_param) else array(dim = 0),
            # structural_beta_cluster_raw = if (use_cluster_effects) matrix(rnorm(num_clusters * num_beta_param), nrow = num_clusters, ncol = num_beta_param) else array(dim = c(0, num_beta_param)),
            # structural_beta_cluster_sd = if (use_cluster_effects) abs(rnorm(num_beta_param, 0, 0.25)) else array(dim = 0),
            structural_beta_cluster = if (use_cluster_effects) matrix(rnorm(num_clusters * num_treatments), nrow = num_clusters, ncol = num_treatments) else array(dim = 0),
            structural_beta_cluster_raw = if (use_cluster_effects) matrix(rnorm(num_clusters * num_treatments), nrow = num_clusters, ncol = num_treatments) else array(dim = c(0, num_treatments)),
            structural_beta_cluster_sd = if (use_cluster_effects) abs(rnorm(num_treatments, 0, 0.25)) else array(dim = 0),
            structural_cluster_takeup_prob = matrix(rbeta(num_clusters * num_mix, 10, 10), nrow = num_mix),
            lambda_v_mix = rep(1 /num_mix, num_mix),
            v_mix_mean = as.array(0.1),
            
            # mu_rep_raw = rep(0.05, num_treatments),
            # mu_rep = rep(0.1, num_treatments),
            structural_beta = rep(0.2, num_beta_param),
            # dist_cost_k_raw = rep(0.05, num_treatments),
            # dist_cost_k = rep(0.05, num_treatments),
            # structural_beta_cluster = if (use_cluster_effects) matrix(rep(0.2, num_clusters * num_treatments), nrow = num_clusters, ncol = num_treatments) else array(dim = 0),
            # structural_beta_cluster_raw = if (use_cluster_effects) matrix(rep(0.2, num_clusters * num_treatments), nrow = num_clusters, ncol = num_treatments) else array(dim = c(0, 4)),
            # structural_beta_cluster_sd = if (use_cluster_effects) rep(0.1, num_treatments) else array(dim = 0),
            # structural_cluster_takeup_prob = rep(0.5, num_clusters) 
          )
                      
      } else {
        return(base_list) 
      }
    }
  } else return(NULL)
}

models <- lst(
  NO_DIST = lst(model_type = 1,
                model_file = "dist_splines3.stan"),
  
  DISCRETE_DIST = lst(model_type = 2),
  
  LINEAR_DIST = lst(model_type = 3,
                    use_dist_cluster_effects = FALSE,
                    model_file = "dist_splines3.stan"),
  
  QUADRATIC_NONLINEAR_DIST = lst(model_type = 4,
                                 use_dist_cluster_effects = FALSE,
                                 model_file = "dist_splines3.stan"),
  
  CUBIC_NONLINEAR_DIST = lst(model_type = 5,
                             use_dist_cluster_effects = FALSE),
  
  SEMIPARAM_NONLINEAR_DIST_OSULLIVAN = lst(model_type = 6,
                                           model_file = "dist_splines3.stan",
                                           use_dist_cluster_effects = FALSE,
                                           init = generate_initializer(structural_type = 0)),
                                           # init = function() lst(u_splines_cluster_v_sigma = abs(rnorm(num_clusters, 0, 0.25)))),
  
  SEMIPARAM_NONLINEAR_DIST_BSPLINE = lst(model_type = 7, 
                                         model_file = "dist_splines3.stan",
                                         Z_splines_v = Z_b_spline,
                                         num_knots_v = ncol(Z_splines_v),
                                         Z_grid_v = Z_grid_b_spline,
                                         iter = 16000,
                                         init = function() lst(u_splines_cluster_v_sigma = abs(rnorm(num_clusters, 0, 0.25)))),
  
  SEMIPARAM_NONLINEAR_DIST_ISPLINE = lst(model_type = 8,
                                         model_file = "dist_splines3.stan",
                                         Z_splines_v = Z_i_spline,
                                         num_knots_v = ncol(Z_splines_v),
                                         Z_grid_v = Z_grid_i_spline,
                                         init = function() lst(u_splines_cluster_v_sigma = abs(rnorm(num_clusters, 0, 0.25)))),
  
  GP = lst(model_type = 9,
           model_file = "dist_splines3.stan",
           use_dist_cluster_effects = FALSE,
           init = generate_initializer(structural_type = 0)),
  
  STRUCTURAL = lst(model_type = 10,
                   model_file = "dist_struct_fixedpoint.stan",
                   num_v_mix = 1,
                   use_semiparam_cost = 1,
                   use_private_incentive_restrictions = 0,
                   suppress_reputation = 0,
                   simulate_new_data,
                   iter = 4000,
                   init = generate_initializer(structural_type = 1, 
                                               num_mix = num_v_mix, 
                                               restricted_private_incentive = use_private_incentive_restrictions,
                                               semiparam = use_semiparam_cost,
                                               suppress_reputation = suppress_reputation)),
  
  STRUCTURAL_RESTRICTED_PRIVATE = 
    lst(model_type = 10,
        model_file = "dist_struct_fixedpoint.stan",
        control = lst(max_treedepth = 12),
        num_v_mix = 1,
        use_semiparam_cost = 1,
        use_private_incentive_restrictions = 1,
        suppress_reputation = 0,
        simulate_new_data,
        iter = 4000,
        init = generate_initializer(structural_type = 1, 
                                    num_mix = num_v_mix, 
                                    restricted_private_incentive = use_private_incentive_restrictions,
                                    semiparam = use_semiparam_cost,
                                    suppress_reputation = suppress_reputation)),
  
  STRUCTURAL_NO_REPUTATION = 
    lst(model_type = 10,
        model_file = "dist_struct_fixedpoint.stan",
        num_v_mix = 1,
        use_semiparam_cost = 1,
        use_private_incentive_restrictions = 1,
        suppress_reputation = 1,
        simulate_new_data,
        iter = 4000,
        init = generate_initializer(structural_type = 1, 
                                    num_mix = num_v_mix, 
                                    restricted_private_incentive = use_private_incentive_restrictions,
                                    semiparam = use_semiparam_cost,
                                    suppress_reputation = suppress_reputation)),
  
  STRUCTURAL_PARAM = lst(model_type = 10,
                   model_file = "dist_struct_fixedpoint.stan",
                   num_v_mix = 1,
                   use_semiparam_cost = 0,
                   use_private_incentive_restrictions = 0,
                   suppress_reputation = 0,
                   simulate_new_data,
                   iter = 4000,
                   init = generate_initializer(structural_type = 1, 
                                               num_mix = num_v_mix, 
                                               restricted_private_incentive = use_private_incentive_restrictions,
                                               semiparam = use_semiparam_cost,
                                               suppress_reputation = suppress_reputation)),
  
  STRUCTURAL_PARAM_RESTRICTED_PRIVATE = 
    lst(model_type = 10,
        model_file = "dist_struct_fixedpoint.stan",
        num_v_mix = 1,
        use_semiparam_cost = 0,
        use_private_incentive_restrictions = 1,
        suppress_reputation = 0,
        simulate_new_data,
        iter = 4000,
        init = generate_initializer(structural_type = 1, 
                                    num_mix = num_v_mix, 
                                    restricted_private_incentive = use_private_incentive_restrictions,
                                    semiparam = use_semiparam_cost,
                                    suppress_reputation = suppress_reputation)),
  
  STRUCTURAL_PARAM_NO_REPUTATION = 
    lst(model_type = 10,
        model_file = "dist_struct_fixedpoint.stan",
        num_v_mix = 1,
        use_semiparam_cost = 0,
        use_private_incentive_restrictions = 1,
        suppress_reputation = 1,
        simulate_new_data,
        iter = 4000,
        init = generate_initializer(structural_type = 1, 
                                    num_mix = num_v_mix, 
                                    restricted_private_incentive = use_private_incentive_restrictions,
                                    semiparam = use_semiparam_cost,
                                    suppress_reputation = suppress_reputation)),
  
  LINEAR_DIST_CE = lst(model_type = 3,
                       model_file = "dist_splines3.stan",
                       use_dist_cluster_effects = TRUE,
                       init = function() lst(dist_beta_cluster_linear_sd = abs(rnorm(num_treatments, 0, 0.25)))),
  
  QUADRATIC_NONLINEAR_DIST_CE = lst(model_type = 4,
                                    model_file = "dist_splines3.stan",
                                    use_dist_cluster_effects = TRUE,
                                    init = function() lst(dist_beta_cluster_linear_sd = abs(rnorm(num_treatments, 0, 0.25)),
                                                          dist_beta_cluster_quadratic_sd = abs(rnorm(num_treatments, 0, 0.25)))),
  
  CUBIC_NONLINEAR_DIST_CE = lst(model_type = 5,
                                model_file = "dist_splines3.stan",
                                use_dist_cluster_effects = TRUE,
                                init = function() lst(dist_beta_cluster_linear_sd = abs(rnorm(num_treatments, 0, 0.25)),
                                                      dist_beta_cluster_quadratic_sd = abs(rnorm(num_treatments, 0, 0.25)),
                                                      dist_beta_cluster_cubic_sd = abs(rnorm(num_treatments, 0, 0.25)))),
  
)

# Stan Run ----------------------------------------------------------------

stan_data <- lst(
  num_obs = nrow(monitored_nosms_data),
  num_grid_obs = length(grid_dist),
  num_treatments = 4,
  num_clusters = n_distinct(monitored_nosms_data$cluster_id),
  obs_cluster_id = monitored_nosms_data$cluster_id,
  cluster_assigned_treatment = distinct(monitored_nosms_data, cluster_id, assigned.treatment) %>% 
    arrange(cluster_id) %>% 
    pull(assigned.treatment),
  takeup = monitored_nosms_data$dewormed,
  cluster_standard_dist = distinct(monitored_nosms_data, cluster_id, standard_cluster.dist.to.pot) %>% 
    arrange(cluster_id) %>% 
    pull(standard_cluster.dist.to.pot),
  cluster_assigned_dist_group_treatment = distinct(monitored_nosms_data, cluster_id, assigned.treatment, dist.pot.group) %>% 
    arrange(cluster_id) %>% 
    group_indices(assigned.treatment, dist.pot.group),
  num_dist_group_treatments = n_distinct(cluster_assigned_dist_group_treatment),
  grid_dist,
  grid_dist2,
  
  Z_splines_v = Z_osullivan,
  num_knots_v = ncol(Z_splines_v),
  
  Z_grid_v = Z_grid_osullivan,
  Z_grid_v2 = Z_grid_osullivan2,
  
  u_splines_v_sigma_sd = 1,
 
  num_excluded_clusters = 0,
  excluded_clusters = array(dim = 0),
 
  is_structural = 0, 
  use_binomial = 0,
  use_cluster_effects,
  use_dist_cluster_effects,
  use_cost_k_restrictions = 1,
) %>% 
  list_modify(!!!map(models, pluck, "model_type") %>% set_names(~ str_c("MODEL_TYPE_", .)))

models <- if (!is_null(models_to_run)) {
  models[models_to_run]
} else {
  models[-c(7,8)]
} 

if (script_options$fit) {
  # dist_fit <- models[-c(7,8)] %>% stan_list(stan_data)
  dist_fit <- models %>% stan_list(stan_data)
  
  if (script_options$`update-output`) {
    new_dist_fit <- dist_fit
  
    dist_fit <- tryCatch({
      load(output_file_name)
      
      list_modify(dist_fit, !!!new_dist_fit)
    })  
  }
} else if (script_options$cv) {
  dist_kfold <- models %>% stan_list(stan_data)
  
  if (script_options$`update-output`) {
    new_dist_kfold <- dist_kfold
  
    dist_kfold <- tryCatch({
      load(output_file_name)
      
      list_modify(dist_kfold, !!!new_dist_kfold)
    })  
  }
}

if (!script_options$`no-save`) {
  if (script_options$fit) {
    save(dist_fit, models, monitored_nosms_data, cluster_analysis_data, grid_dist, extract_sim_level, extract_obs_fit_level, extract_sim_diff, standardize, unstandardize, stan_data,
         file = output_file_name)
  } else if (script_options$cv) {
    save(dist_kfold, models, monitored_nosms_data, 
         file = output_file_name)
  }
}

