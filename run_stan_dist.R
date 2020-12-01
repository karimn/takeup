#!/usr/bin/Rscript

script_options <- docopt::docopt(
"Usage:
  run_stan_dist fit [--no-save --sequential --chains=<chains> --iter=<iter> --thin=<thin> --force-iter --models=<models> --outputname=<output file name> --update-output --predict-prior]
  run_stan_dist cv [--folds=<number of folds> --no-save --sequential --chains=<chains> --iter=<iter> --thin=<thin> --force-iter --models=<models> --outputname=<output file name> --update-output]
  
Options:
  --folds=<number of folds>  Cross validation folds [default: 10]
  --chains=<chains>  Number of Stan chains [default: 4]
  --iter=<iter>  Number of (warmup + sampling) iterations [default: 8000]
  --thin=<thin>  Thin samples [default: 1]",

  # args = if (interactive()) "fit --sequential --outputname=dist_fit28 --update-output" else commandArgs(trailingOnly = TRUE) 
  args = if (interactive()) "fit --sequential --outputname=dist_fit29 --models=5" else commandArgs(trailingOnly = TRUE) 
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

folds <- as.integer(script_options$folds %||% 10) # CV k-folds
chains <- as.integer(script_options$chains) # Stan chains
iter <- as.integer(script_options$iter) # Stan iterations
output_name <- if (!is_null(script_options$outputname)) { script_options$outputname } else if (script_options$fit) { "dist_fit" } else { "dist_kfold" }
output_file_name <- file.path("data", "stan_analysis_data", str_c(output_name, ".RData"))
models_to_run <- if (!is_null(script_options$models)) c(str_split(script_options$models, ",", simplify = TRUE)) %>% as.integer()
thin_by <- as.integer(script_options$thin)
predict_prior <- as.logical(script_options$predict_prior) # Prior prediction run

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

analysis_data <- monitored_nosms_data

# Splines -----------------------------------------------------------------

# These are used for semiparameteric models. Currently, we are not using any semiparameteric models.

num_interior_knots <- 100

cluster_standard_dist <- distinct(analysis_data, cluster_id, standard_cluster.dist.to.pot) %>% 
  arrange(cluster_id) %>% 
  pull(standard_cluster.dist.to.pot)

Z_osullivan <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, spline_type = "osullivan")
Z_i_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, spline_type = "i-spline")
Z_b_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, spline_type = "b-spline")

grid_dist <- get_spline_range(cluster_standard_dist) %>% unname() %>% list_modify(length = 1001) %>% do.call(seq, .)

Z_grid_osullivan <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, splines_for = grid_dist, spline_type = "osullivan")
Z_grid_i_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, splines_for = grid_dist, spline_type = "i-spline")
Z_grid_b_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, splines_for = grid_dist, spline_type = "b-spline")

treatment_dist <- distinct(analysis_data, assigned.treatment, standard_cluster.dist.to.pot) %>% 
  group_by(assigned.treatment) %>% 
  group_map(~ pull(.x, "standard_cluster.dist.to.pot"))

grid_dist2 <- treatment_dist %>% 
  map(get_spline_range) %>% 
  map(unname) %>% 
  map(list_modify, length = 1001) %>% 
  map(~ do.call(seq, .))

Z_grid_osullivan2 <- map2(treatment_dist, grid_dist2, ~ calculate_splines(.x, num_interior_knots = num_interior_knots, splines_for = .y, spline_type = "osullivan")) 

# Models ------------------------------------------------------------------

num_treatments <- n_distinct(analysis_data$assigned.treatment)
num_clusters <- n_distinct(analysis_data$cluster_id)
num_counties <- n_distinct(analysis_data$county)

# Models not using "takeup_struct.stan" can be ignored

struct_model_stan_pars <- c(
  "total_error_sd", "cluster_dist_cost", "structural_cluster_benefit_cost", "structural_cluster_obs_v", "structural_cluster_takeup_prob",
  "beta", "dist_beta_v", "dist_quadratic_beta_v", "mu_rep", "cluster_cf_benefit_cost", "mu_cluster_effects_raw", "mu_cluster_effects_sd", "cluster_mu_rep", # "linear_dist_cost", 
  "cluster_rep_benefit_cost", "sim_benefit_cost",
  "group_dist_mean", "group_dist_sd", "group_dist_mix",
  "dist_beta_county_raw", "dist_beta_county_sd")

reduced_model_stan_pars <- c(
  "structural_cluster_benefit_cost", "structural_cluster_obs_v", "structural_cluster_takeup_prob",
  "beta", "cluster_cf_benefit_cost", 
  "cluster_rep_benefit_cost")

models <- lst(
  STRUCTURAL_QUADRATIC = lst(
    model_file = "takeup_struct.stan",
    pars = struct_model_stan_pars,
    control = lst(max_treedepth = 12, adapt_delta = 0.99),
    use_binomial = FALSE,
    num_v_mix = 1,
    use_cost_model = cost_model_types["param_quadratic"],
    use_single_cost_model = TRUE,
    use_private_incentive_restrictions = TRUE,
    use_salience_effect = FALSE,
    use_cluster_effects = TRUE,
    use_county_effects = TRUE,
    use_param_dist_cluster_effects = FALSE,
    use_param_dist_county_effects = FALSE,
    use_mu_cluster_effects = FALSE,
    use_mu_county_effects = FALSE,
    use_shifting_v_dist = FALSE,
    suppress_reputation = FALSE,
    suppress_shocks = FALSE,
    generate_sim = FALSE,
    iter = 2000,
    thin = 1,
   
    # Priors 
    mu_rep_sd = 1,
    structural_beta_county_sd_sd = 0.25,
    structural_beta_cluster_sd_sd = 0.25,
    
    init = generate_initializer(
      num_treatments = num_treatments, 
      num_clusters = num_clusters,
      num_counties = num_counties,
      structural_type = 1, 
      num_mix = num_v_mix, 
      use_cluster_effects = use_cluster_effects,
      use_county_effects = use_county_effects,
      use_param_dist_cluster_effects = use_param_dist_cluster_effects,
      use_mu_cluster_effects = use_mu_cluster_effects,
      use_mu_county_effects = use_mu_county_effects,
      restricted_private_incentive = use_private_incentive_restrictions,
      cost_model_type = use_cost_model,
      use_single_cost_model = use_single_cost_model,
      num_knots = ncol(Z_osullivan),
      name_matched = FALSE,
      suppress_reputation = suppress_reputation)) %>% 
    list_modify(!!!enum2stan_data(cost_model_types)),
  
  # STRUCTURAL_SEMIPARAM = lst(
  #   model_type = 10,
  #   model_file = "takeup_struct.stan",
  #   control = lst(max_treedepth = 12, adapt_delta = 0.99),
  #   use_binomial = FALSE,
  #   num_v_mix = 1,
  #   use_cost_model = cost_model_types["semiparam"],
  #   use_single_cost_model = TRUE,
  #   use_private_incentive_restrictions = TRUE,
  #   use_salience_effect = FALSE,
  #   use_cluster_effects = TRUE,
  #   use_county_effects = TRUE,
  #   use_param_dist_cluster_effects = FALSE,
  #   use_param_dist_county_effects = FALSE,
  #   use_mu_cluster_effects = FALSE,
  #   use_mu_county_effects = FALSE,
  #   use_shifting_v_dist = FALSE,
  #   suppress_reputation = FALSE,
  #   suppress_shocks = TRUE,
  #   # simulate_new_data,
  #   iter = 1000,
  #   thin = 1,
  #   init = generate_initializer(
  #     num_treatments = num_treatments, 
  #     num_clusters = num_clusters,
  #     num_counties = num_counties,
  #     structural_type = 1, 
  #     num_mix = num_v_mix, 
  #     use_cluster_effects = use_cluster_effects,
  #     use_param_dist_cluster_effects = use_param_dist_cluster_effects,
  #     use_mu_cluster_effects = use_mu_cluster_effects,
  #     use_mu_county_effects = use_mu_county_effects,
  #     restricted_private_incentive = use_private_incentive_restrictions,
  #     cost_model_type = use_cost_model,
  #     use_single_cost_model = use_single_cost_model,
  #     num_knots = ncol(Z_osullivan),
  #     name_matched = FALSE,
  #     suppress_reputation = suppress_reputation)) %>% 
  #   list_modify(!!!enum2stan_data(cost_model_types)),
  
  STRUCTURAL_LINEAR = lst(
    model_file = "takeup_struct.stan",
    pars = struct_model_stan_pars,
    control = lst(max_treedepth = 12, adapt_delta = 0.99),
    use_binomial = FALSE,
    num_v_mix = 1,
    use_cost_model = cost_model_types["param_linear"],
    use_single_cost_model = TRUE,
    use_private_incentive_restrictions = TRUE,
    use_salience_effect = FALSE,
    use_cluster_effects = TRUE,
    use_county_effects = TRUE,
    use_param_dist_cluster_effects = FALSE,
    use_param_dist_county_effects = FALSE,
    use_mu_cluster_effects = FALSE,
    use_mu_county_effects = FALSE,
    use_shifting_v_dist = FALSE,
    suppress_reputation = FALSE,
    suppress_shocks = FALSE,
    generate_sim = FALSE,
    iter = 2000,
    thin = 1,
    
    mu_rep_sd = 1,
    structural_beta_county_sd_sd = 0.25,
    structural_beta_cluster_sd_sd = 0.25,
    
    init = generate_initializer(
      num_treatments = num_treatments, 
      num_clusters = num_clusters,
      num_counties = num_counties,
      structural_type = 1, 
      num_mix = num_v_mix, 
      use_cluster_effects = use_cluster_effects,
      use_county_effects = use_county_effects,
      use_param_dist_cluster_effects = use_param_dist_cluster_effects,
      use_mu_cluster_effects = use_mu_cluster_effects,
      use_mu_county_effects = use_mu_county_effects,
      restricted_private_incentive = use_private_incentive_restrictions,
      cost_model_type = use_cost_model,
      use_single_cost_model = use_single_cost_model,
      num_knots = ncol(Z_osullivan),
      name_matched = FALSE,
      suppress_reputation = suppress_reputation)) %>% 
    list_modify(!!!enum2stan_data(cost_model_types)),
  
  # REDUCED_FORM_LINEAR = lst(
  #   model_type = 10,
  #   model_file = "takeup_struct.stan",
  #   control = lst(max_treedepth = 12, adapt_delta = 0.99),
  #   use_binomial = FALSE,
  #   num_v_mix = 1,
  #   use_cost_model = cost_model_types["param_linear"],
  #   use_single_cost_model = FALSE,
  #   use_private_incentive_restrictions = FALSE,
  #   use_salience_effect = FALSE,
  #   use_cluster_effects = TRUE,
  #   use_county_effects = TRUE,
  #   use_param_dist_cluster_effects = FALSE,
  #   use_param_dist_county_effects = FALSE,
  #   use_mu_cluster_effects = FALSE,
  #   use_mu_county_effects = FALSE,
  #   use_shifting_v_dist = FALSE,
  #   suppress_reputation = TRUE,
  #   suppress_shocks = TRUE,
  #   # simulate_new_data,
  #   iter = 1000,
  #   thin = 1,
  #   init = generate_initializer(
  #     num_treatments = num_treatments, 
  #     num_clusters = num_clusters,
  #     num_counties = num_counties,
  #     structural_type = 1, 
  #     num_mix = num_v_mix, 
  #     use_cluster_effects = use_cluster_effects,
  #     use_param_dist_cluster_effects = use_param_dist_cluster_effects,
  #     use_mu_cluster_effects = use_mu_cluster_effects,
  #     use_mu_county_effects = use_mu_county_effects,
  #     restricted_private_incentive = use_private_incentive_restrictions,
  #     cost_model_type = use_cost_model,
  #     use_single_cost_model = use_single_cost_model,
  #     num_knots = ncol(Z_osullivan),
  #     name_matched = FALSE,
  #     suppress_reputation = suppress_reputation)) %>% 
  #   list_modify(!!!enum2stan_data(cost_model_types)),
  
  STRUCTURAL_QUADRATIC_SALIENCE = lst(
    model_type = 10,
    model_file = "takeup_struct.stan",
    control = lst(max_treedepth = 12, adapt_delta = 0.99),
    num_v_mix = 1,
    use_single_cost_model = TRUE,
    use_cost_model = cost_model_types["param_quadratic_salience"],
    use_private_incentive_restrictions = TRUE,
    use_salience_effect = TRUE,
    use_cluster_effects = TRUE,
    use_county_effects = TRUE,
    use_param_dist_cluster_effects = FALSE,
    use_param_dist_county_effects = FALSE,
    use_mu_cluster_effects = FALSE,
    use_mu_county_effects = FALSE,
    use_shifting_v_dist = FALSE,
    suppress_reputation = FALSE,
    suppress_shocks = FALSE,
    
    mu_rep_sd = 1,
    structural_beta_county_sd_sd = 0.25,
    structural_beta_cluster_sd_sd = 0.25,
    
    iter = 2000,
    thin = 1,
    init = generate_initializer(
      num_treatments = num_treatments, 
      num_clusters = num_clusters,
      num_counties = num_counties,
      structural_type = 1, 
      num_mix = num_v_mix, 
      use_cluster_effects = use_cluster_effects,
      use_param_dist_cluster_effects = use_param_dist_cluster_effects,
      use_mu_cluster_effects = use_mu_cluster_effects,
      restricted_private_incentive = use_private_incentive_restrictions,
      cost_model_type = use_cost_model,
      use_single_cost_model = TRUE,
      num_knots = ncol(Z_osullivan),
      name_matched = FALSE,
      suppress_reputation = suppress_reputation)) %>% 
    list_modify(!!!enum2stan_data(cost_model_types)),
  
  # STRUCTURAL_LINEAR_SALIENCE = lst(
  #   model_type = 10,
  #   model_file = "takeup_struct.stan",
  #   control = lst(max_treedepth = 12, adapt_delta = 0.99),
  #   num_v_mix = 1,
  #   use_single_cost_model = TRUE,
  #   use_cost_model = cost_model_types["param_linear_salience"],
  #   use_private_incentive_restrictions = TRUE,
  #   use_salience_effect = TRUE,
  #   use_cluster_effects = TRUE,
  #   use_county_effects = TRUE,
  #   use_param_dist_cluster_effects = FALSE,
  #   use_param_dist_county_effects = FALSE,
  #   use_mu_cluster_effects = FALSE,
  #   use_mu_county_effects = FALSE,
  #   use_shifting_v_dist = FALSE,
  #   suppress_reputation = FALSE,
  #   suppress_shocks = FALSE,
  #   # simulate_new_data,
  #   iter = 2000,
  #   thin = 1,
  #   init = generate_initializer(
  #     num_treatments = num_treatments, 
  #     num_clusters = num_clusters,
  #     num_counties = num_counties,
  #     structural_type = 1, 
  #     num_mix = num_v_mix, 
  #     use_cluster_effects = use_cluster_effects,
  #     use_param_dist_cluster_effects = use_param_dist_cluster_effects,
  #     use_mu_cluster_effects = use_mu_cluster_effects,
  #     restricted_private_incentive = use_private_incentive_restrictions,
  #     cost_model_type = use_cost_model,
  #     use_single_cost_model = TRUE,
  #     num_knots = ncol(Z_osullivan),
  #     name_matched = FALSE,
  #     suppress_reputation = suppress_reputation)) %>% 
  #   list_modify(!!!enum2stan_data(cost_model_types)),
  
  REDUCED_FORM_NO_RESTRICT = lst(
    model_type = 10,
    model_file = "takeup_reduced.stan",
    pars = reduced_model_stan_pars,
    control = lst(max_treedepth = 10, adapt_delta = 0.99),
    num_v_mix = 1,
    use_single_cost_model = FALSE,
    use_cost_model = cost_model_types["discrete"],
    use_private_incentive_restrictions = FALSE,
    use_salience_effect = FALSE,
    use_cluster_effects = TRUE,
    use_county_effects = TRUE,
    use_param_dist_cluster_effects = FALSE,
    use_param_dist_county_effects = FALSE,
    use_mu_cluster_effects = FALSE,
    use_mu_county_effects = FALSE,
    use_shifting_v_dist = FALSE,
    suppress_reputation = TRUE,
    suppress_shocks = TRUE,
    
    structural_beta_county_sd_sd = 1,
    structural_beta_cluster_sd_sd = 1,
    
    iter = 2000,
    thin = 1,
    init = generate_initializer(
      num_treatments = num_treatments, 
      num_clusters = num_clusters,
      num_counties = num_counties,
      structural_type = 1, 
      num_mix = num_v_mix, 
      use_cluster_effects = use_cluster_effects,
      use_county_effects = use_county_effects,
      use_param_dist_cluster_effects = use_param_dist_cluster_effects,
      use_mu_cluster_effects = use_mu_cluster_effects,
      restricted_private_incentive = use_private_incentive_restrictions,
      cost_model_type = use_cost_model,
      use_single_cost_model = TRUE,
      num_knots = ncol(Z_osullivan),
      name_matched = FALSE,
      suppress_reputation = suppress_reputation)) %>% 
    list_modify(!!!enum2stan_data(cost_model_types)),
  
  # REDUCED_FORM_NO_LEVELS = lst(
  #   model_type = 10,
  #   model_file = "takeup_struct.stan",
  #   control = lst(max_treedepth = 12, adapt_delta = 0.99),
  #   num_v_mix = 1,
  #   use_single_cost_model = TRUE,
  #   use_cost_model = cost_model_types["discrete"],
  #   use_private_incentive_restrictions = FALSE,
  #   use_salience_effect = FALSE,
  #   use_cluster_effects = FALSE,
  #   use_county_effects = FALSE,
  #   use_param_dist_cluster_effects = FALSE,
  #   use_param_dist_county_effects = FALSE,
  #   use_mu_cluster_effects = FALSE,
  #   use_mu_county_effects = FALSE,
  #   use_shifting_v_dist = FALSE,
  #   suppress_reputation = TRUE,
  #   suppress_shocks = TRUE,
  #   # simulate_new_data,
  #   iter = 1500,
  #   thin = 1,
  #   init = generate_initializer(
  #     num_treatments = num_treatments, 
  #     num_clusters = num_clusters,
  #     num_counties = num_counties,
  #     structural_type = 1, 
  #     num_mix = num_v_mix, 
  #     use_cluster_effects = use_cluster_effects,
  #     use_param_dist_cluster_effects = use_param_dist_cluster_effects,
  #     use_mu_cluster_effects = use_mu_cluster_effects,
  #     restricted_private_incentive = use_private_incentive_restrictions,
  #     cost_model_type = use_cost_model,
  #     use_single_cost_model = TRUE,
  #     num_knots = ncol(Z_osullivan),
  #     name_matched = FALSE,
  #     suppress_reputation = suppress_reputation)) %>% 
  #   list_modify(!!!enum2stan_data(cost_model_types)),
  
  STRUCTURAL_QUADRATIC_SALIENCE_NO_REP = lst(
    model_type = 10,
    model_file = "takeup_struct.stan",
    control = lst(max_treedepth = 12, adapt_delta = 0.99),
    num_v_mix = 1,
    use_single_cost_model = TRUE,
    use_cost_model = cost_model_types["param_quadratic_salience"],
    use_private_incentive_restrictions = TRUE,
    use_salience_effect = TRUE,
    use_cluster_effects = TRUE,
    use_county_effects = TRUE,
    use_param_dist_cluster_effects = FALSE,
    use_param_dist_county_effects = FALSE,
    use_mu_cluster_effects = FALSE,
    use_mu_county_effects = FALSE,
    use_shifting_v_dist = FALSE,
    suppress_reputation = TRUE,
    suppress_shocks = FALSE,
    
    mu_rep_sd = 1,
    structural_beta_county_sd_sd = 0.25,
    structural_beta_cluster_sd_sd = 0.25,
    
    iter = 4000,
    thin = 1,
    init = generate_initializer(
      num_treatments = num_treatments, 
      num_clusters = num_clusters,
      num_counties = num_counties,
      structural_type = 1, 
      num_mix = num_v_mix, 
      use_cluster_effects = use_cluster_effects,
      use_param_dist_cluster_effects = use_param_dist_cluster_effects,
      use_mu_cluster_effects = use_mu_cluster_effects,
      restricted_private_incentive = use_private_incentive_restrictions,
      cost_model_type = use_cost_model,
      use_single_cost_model = TRUE,
      num_knots = ncol(Z_osullivan),
      name_matched = FALSE,
      suppress_reputation = suppress_reputation)) %>% 
    list_modify(!!!enum2stan_data(cost_model_types)),
)

# Stan Run ----------------------------------------------------------------

stan_data <- lst(
  num_obs = nrow(analysis_data),
  num_grid_obs = length(grid_dist),
  num_treatments,
  is_name_matched = !analysis_data$monitored,
  num_clusters,
  num_counties,
  obs_cluster_id = analysis_data$cluster_id,
  obs_county_id = as.integer(analysis_data$county),
  cluster_county_id = analysis_data %>% 
    distinct(cluster_id, county) %>% 
    arrange(cluster_id) %>% 
    pull(county) %>% 
    as.integer(),
  cluster_assigned_treatment = distinct(analysis_data, cluster_id, assigned.treatment) %>% 
    arrange(cluster_id) %>% 
    pull(assigned.treatment),
  takeup = analysis_data$dewormed,
  
  cluster_standard_dist = distinct(analysis_data, cluster_id, standard_cluster.dist.to.pot) %>% 
    arrange(cluster_id) %>% 
    pull(standard_cluster.dist.to.pot),
  
  cluster_treatment_map = distinct(analysis_data, assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>% 
    arrange(assigned_dist_group, assigned_treatment),
  
  cluster_assigned_dist_group = distinct(analysis_data, cluster_id, dist.pot.group) %>% 
    arrange(cluster_id) %>% 
    pull(dist.pot.group),
  
  cluster_assigned_dist_group_treatment = distinct(analysis_data, cluster_id, assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>% 
    left_join(cluster_treatment_map %>% mutate(treatment_id = seq(n())), by = c("assigned_treatment", "assigned_dist_group")) %>% 
    arrange(cluster_id) %>% 
    pull(treatment_id),
  
  num_dist_group_treatments = n_distinct(cluster_assigned_dist_group_treatment),
  grid_dist,
  grid_dist2,
  
  num_discrete_dist = 2,
  
  small_grid_dist = grid_dist[(seq_along(grid_dist) %% 100) == 0],
  num_small_grid_obs = length(small_grid_dist),
  
  Z_splines_v = Z_osullivan,
  num_knots_v = ncol(Z_splines_v),
  
  Z_grid_v = Z_grid_osullivan,
  Z_grid_v2 = Z_grid_osullivan2,
  
 
  num_excluded_clusters = 0,
  excluded_clusters = array(dim = 0),

  use_name_matched_obs = FALSE, 
  is_structural = FALSE, 
  use_binomial = FALSE,
  use_single_cost_model = FALSE,
  use_cost_k_restrictions = TRUE,
  use_shifting_v_dist = FALSE,
  suppress_shocks = FALSE,
  cluster_log_lik = TRUE,
  generate_rep = TRUE,
  generate_sim = FALSE,
  predict_prior,
  
  thin = thin_by,
  
  # Priors
  
  u_splines_v_sigma_sd = 1,
  mu_rep_sd = 1,
  
  beta_control_sd = 5,
  beta_far_effect_sd = 1,
  beta_ink_effect_sd = 1,
  beta_calendar_effect_sd = 1,
  beta_bracelet_effect_sd = 1,
  
  structural_beta_county_sd_sd = 0.25,
  structural_beta_cluster_sd_sd = 0.25,
  
  analysis_data
) %>% 
  list_modify(!!!map(models, pluck, "model_type") %>% set_names(~ str_c("MODEL_TYPE_", .)))

models <- if (!is_null(models_to_run)) {
  models[models_to_run]
} else {
  models
}

if (script_options$fit) {
  dist_fit <- models %>% stan_list(stan_data)
  
  if (script_options$update_output) {
    new_dist_fit <- dist_fit
  
    dist_fit <- tryCatch({
      old_data_env <- new.env()
      load(output_file_name, envir = old_data_env)
      
      list_modify(old_data_env$dist_fit, !!!new_dist_fit)
    }, error = function(e) dist_fit)  
  }
} else if (script_options$cv) {
  dist_kfold <- models %>% stan_list(stan_data)
  
  if (script_options$update_output) {
    new_dist_kfold <- dist_kfold
  
    dist_kfold <- tryCatch({
      old_data_env <- new.env()
      load(output_file_name, envir = old_data_env)
      
      list_modify(old_data_env$dist_kfold, !!!new_dist_kfold)
    })  
  }
}

if (!script_options$no_save) {
  if (script_options$fit) {
    save(dist_fit, models, grid_dist, stan_data, file = output_file_name)
  } else if (script_options$cv) {
    save(dist_kfold, models, file = output_file_name)
  }
}

