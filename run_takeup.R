#!/usr/bin/Rscript

script_options <- docopt::docopt(
  stringr::str_glue("Usage:
  run_takeup.R takeup prior [--no-save --sequential --chains=<chains> --threads=<threads> --iter=<iter> --thin=<thin> --force-iter --models=<models> --outputname=<output file name> --update-output --cmdstanr --include-paths=<paths> --output-path=<path>]
  run_takeup.R takeup fit [--no-save --sequential --chains=<chains> --threads=<threads> --iter=<iter> --thin=<thin> --force-iter --models=<models> --outputname=<output file name> --update-output --cmdstanr --include-paths=<paths> --output-path=<path>]
  run_takeup.R takeup cv [--folds=<number of folds> --no-save --sequential --chains=<chains> --threads=<threads> --iter=<iter> --thin=<thin> --force-iter --models=<models> --outputname=<output file name> --update-output --cmdstanr --include-paths=<paths> --output-path=<path>]
  
  run_takeup.R beliefs fit [--no-save --chains=<chains> --iter=<iter> --outputname=<output file name> --include-paths=<paths> --output-path=<path>]
  
Options:
  --folds=<number of folds>  Cross validation folds [default: 10]
  --chains=<chains>  Number of Stan chains [default: 4]
  --threads=<threads>  Number of threads per chain [default: 1]
  --iter=<iter>  Number of (warmup + sampling) iterations [default: 8000]
  --thin=<thin>  Thin samples [default: 1]
  --include-paths=<paths>  Includes path for cmdstanr [default: .]
  --output-path=<path>  Where to save output files [default: {file.path('data', 'stan_analysis_data')}]
"),

  # args = if (interactive()) "fit --sequential --outputname=dist_fit28 --update-output" else commandArgs(trailingOnly = TRUE) 
  args = if (interactive()) "takeup prior --sequential --outputname=test --output-path=~/Code/takeup/data/stan_analysis_data --models=STRUCTURAL_LINEAR_U_SHOCKS --cmdstanr --include-paths=~/Code/takeup/stan_models --force-iter --iter=100" else commandArgs(trailingOnly = TRUE)
  # args = if (interactive()) "beliefs fit --chains=8 --outputname=test --output-path=~/Code/takeup/data/stan_analysis_data --include-paths=~/Code/takeup/stan_models --iter=1000" else commandArgs(trailingOnly = TRUE) 
) 

library(magrittr)
library(tidyverse)
library(parallel)
library(pbmcapply)
library(HRW)
library(splines2)
library(loo)

script_options %<>% 
  modify_at(c("chains", "iter", "threads"), as.integer) %>% 
  modify_at(c("models"), ~ c(str_split(script_options$models, ",", simplify = TRUE)))

if (script_options$cmdstanr || script_options$beliefs) {
  library(cmdstanr)
} else {
  library(rstan)
  
  rstan_options(auto_write = TRUE)
}

folds <- as.integer(script_options$folds %||% 10) # CV k-folds
chains <- as.integer(script_options$chains) # Stan chains
iter <- as.integer(script_options$iter) # Stan iterations
output_name <- if (!is_null(script_options$outputname)) { 
  script_options$outputname 
} else if (script_options$takeup) { 
  if (script_options$fit) { 
    "dist_fit" 
  } else { 
    "dist_kfold" 
  }
} else if (script_options$beliefs) {
  "beliefs"
}

output_file_name <- file.path(script_options$output_path, str_c(output_name, ".RData"))
thin_by <- as.integer(script_options$thin)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

# Data --------------------------------------------------------------------

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_nosms_data %>% 
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group)

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
  "structural_cluster_benefit", "cluster_linear_dist_cost", "cluster_quadratic_dist_cost",
  "beta", "dist_beta_v", "dist_quadratic_beta_v", "mu_rep", "cluster_cf_benefit_cost", "mu_cluster_effects_raw", "mu_cluster_effects_sd", "cluster_mu_rep", # "linear_dist_cost", 
  "cluster_rep_benefit_cost", "sim_benefit_cost",
  "group_dist_mean", "group_dist_sd", "group_dist_mix",
  "dist_beta_county_raw", "dist_beta_county_sd")

reduced_model_stan_pars <- c(
  "structural_cluster_benefit_cost", "structural_cluster_obs_v", "structural_cluster_takeup_prob",
  "beta", "cluster_cf_benefit_cost", 
  "cluster_rep_benefit_cost")

models <- lst(
  # STRUCTURAL_QUADRATIC = lst(
  #   model_file = "takeup_struct.stan",
  #   pars = struct_model_stan_pars,
  #   control = lst(max_treedepth = 12, adapt_delta = 0.99),
  #   use_binomial = FALSE,
  #   num_v_mix = 1,
  #   use_cost_model = cost_model_types["param_quadratic"],
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
  #   generate_sim = FALSE,
  #   iter = 4000,
  #   thin = 1,
  #  
  #   # Priors 
  #   mu_rep_sd = 1,
  #   structural_beta_county_sd_sd = 0.25,
  #   structural_beta_cluster_sd_sd = 0.25,
  #   
  #   init = generate_initializer(
  #     num_treatments = num_treatments, 
  #     num_clusters = num_clusters,
  #     num_counties = num_counties,
  #     structural_type = 1, 
  #     num_mix = num_v_mix, 
  #     use_cluster_effects = use_cluster_effects,
  #     use_county_effects = use_county_effects,
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
  
  STRUCTURAL_LINEAR_U_SHOCKS = lst(
    model_file = "takeup_struct.stan",
    pars = struct_model_stan_pars,
    control = lst(max_treedepth = 12, adapt_delta = 0.99),
    use_binomial = FALSE,
    num_v_mix = 1,
    use_cost_model = cost_model_types["param_linear"],
    use_single_cost_model = TRUE,
    use_private_incentive_restrictions = TRUE,
    use_salience_effect = FALSE,
    use_cluster_effects = FALSE,
    use_county_effects = FALSE,
    use_param_dist_cluster_effects = FALSE,
    use_param_dist_county_effects = FALSE,
    use_mu_cluster_effects = FALSE,
    use_mu_county_effects = FALSE, 
    use_restricted_mu = TRUE,
    use_shifting_v_dist = FALSE,
    use_u_in_delta = TRUE,
    use_wtp_model = TRUE,
    use_homoskedastic_shocks = TRUE,
    use_strata_levels = use_county_effects, # WTP
    suppress_reputation = FALSE,
    generate_sim = FALSE,
    iter = 800,
    thin = 1,
    alg_sol_f_tol = 0.001,
    alg_sol_max_steps = 1e9L,
    alg_sol_rel_tol = 0.0000001,

    # Priors
    mu_rep_sd = c(1.5, rep(1, num_treatments - 1)), # The first is control
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
    generate_sim = FALSE,
    iter = 4000,
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
  
  # STRUCTURAL_LINEAR = lst(
  #   model_file = "takeup_struct.stan",
  #   pars = struct_model_stan_pars,
  #   control = lst(max_treedepth = 12, adapt_delta = 0.99),
  #   use_binomial = FALSE,
  #   num_v_mix = 1,
  #   use_cost_model = cost_model_types["param_linear"],
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
  #   generate_sim = FALSE,
  #   iter = 2000,
  #   thin = 1,
  #   
  #   mu_rep_sd = 1,
  #   structural_beta_county_sd_sd = 0.25,
  #   structural_beta_cluster_sd_sd = 0.25,
  #   
  #   init = generate_initializer(
  #     num_treatments = num_treatments, 
  #     num_clusters = num_clusters,
  #     num_counties = num_counties,
  #     structural_type = 1, 
  #     num_mix = num_v_mix, 
  #     use_cluster_effects = use_cluster_effects,
  #     use_county_effects = use_county_effects,
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
  
  # STRUCTURAL_QUADRATIC_SALIENCE = lst(
  #   model_type = 10,
  #   model_file = "takeup_struct.stan",
  #   control = lst(max_treedepth = 12, adapt_delta = 0.99),
  #   num_v_mix = 1,
  #   use_single_cost_model = TRUE,
  #   use_cost_model = cost_model_types["param_quadratic_salience"],
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
  #   
  #   mu_rep_sd = 1,
  #   structural_beta_county_sd_sd = 0.25,
  #   structural_beta_cluster_sd_sd = 0.25,
  #   
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
  
  STRUCTURAL_LINEAR_SALIENCE = lst(
    model_type = 10,
    model_file = "takeup_struct.stan",
    control = lst(max_treedepth = 12, adapt_delta = 0.99),
    num_v_mix = 1,
    use_single_cost_model = TRUE,
    use_cost_model = cost_model_types["param_linear_salience"],
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
    # simulate_new_data,
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
  
  # STRUCTURAL_QUADRATIC_SALIENCE_NO_REP = lst(
  #   model_type = 10,
  #   model_file = "takeup_struct.stan",
  #   control = lst(max_treedepth = 12, adapt_delta = 0.99),
  #   num_v_mix = 1,
  #   use_single_cost_model = TRUE,
  #   use_cost_model = cost_model_types["param_quadratic_salience"],
  #   use_private_incentive_restrictions = TRUE,
  #   use_salience_effect = TRUE,
  #   use_cluster_effects = TRUE,
  #   use_county_effects = TRUE,
  #   use_param_dist_cluster_effects = FALSE,
  #   use_param_dist_county_effects = FALSE,
  #   use_mu_cluster_effects = FALSE,
  #   use_mu_county_effects = FALSE,
  #   use_shifting_v_dist = FALSE,
  #   suppress_reputation = TRUE,
  #   
  #   mu_rep_sd = 1,
  #   structural_beta_county_sd_sd = 0.25,
  #   structural_beta_cluster_sd_sd = 0.25,
  #   
  #   iter = 4000,
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
)


# WTP Stan Data -----------------------------------------------------------

wtp_stan_data <- analysis.data %>% 
  mutate(stratum = county) %>% 
  prepare_bayes_wtp_data(
    wtp.data,
    
    preference_value_diff = seq(-100, 100, 10), 
    num_preference_value_diff = length(preference_value_diff), 
    
    wtp_utility_df = 3,
    tau_mu_wtp_diff = 100,
    mu_wtp_df_student_t = 7,
    tau_sigma_wtp_diff = 50,
    sigma_wtp_df_student_t = 2.5
  )

# Treatment Details -------------------------------------------------------

treatment_formula <- ~ assigned_treatment * assigned_dist_group 

cluster_treatment_map = distinct(analysis_data, assigned_treatment, assigned_dist_group) %>% 
  arrange(assigned_dist_group, assigned_treatment) 

treatment_map_design_matrix <- cluster_treatment_map %>%
  modelr::model_matrix(treatment_formula)

# Beliefs Data ------------------------------------------------------------

analysis_data %<>% 
  nest_join(
    endline.know.table.data %>% 
      filter(fct_match(know.table.type, "table.A")),
    by = "KEY.individ", 
    name = "knowledge_data"
  ) %>% 
  mutate(
    map_dfr(knowledge_data, ~ {
      tibble(
        obs_know_person = sum(.x$num.recognized),
        obs_know_person_prop = mean(.x$num.recognized),
        knows_other_dewormed = sum(fct_match(.x$dewormed, c("yes", "no")), na.rm = TRUE),
        knows_other_dewormed_yes = sum(fct_match(.x$dewormed, "yes"), na.rm = TRUE),
        thinks_other_knows = sum(fct_match(.x$second.order, c("yes", "no")), na.rm = TRUE),
        thinks_other_knows_yes = sum(fct_match(.x$second.order, "yes"), na.rm = TRUE),
      )
    }
  ))

beliefs_ate_pairs <- cluster_treatment_map %>% 
  mutate(treatment_id = seq(n())) %>% {
  bind_rows(
    left_join(., filter(., fct_match(assigned_treatment, "control")), by = c("assigned_dist_group"), suffix = c("", "_control")) %>% 
      filter(assigned_treatment != assigned_treatment_control) %>% 
      select(treatment_id, treatment_id_control),
    
    left_join(., filter(., fct_match(assigned_dist_group, "close")), by = c("assigned_treatment"), suffix = c("", "_control")) %>% 
      filter(assigned_dist_group != assigned_dist_group_control) %>% 
      select(treatment_id, treatment_id_control),
  )
} %>%
  arrange(treatment_id, treatment_id_control) 

# Stan Data ---------------------------------------------------------------

stan_data <- lst(
  # Beliefs Model 
  beliefs_use_stratum_level = TRUE,
  beliefs_use_cluster_level = TRUE,
  beliefs_use_obs_level = TRUE,
  
  know_table_A_sample_size = 10,
  num_beliefs_obs = filter(analysis_data, obs_know_person > 0) %>% nrow(),
  beliefs_obs_index = mutate(analysis_data, obs_index = seq(n())) %>% 
    filter(obs_know_person > 0) %>% 
    pull(obs_index),
  
  num_recognized = filter(analysis_data, obs_know_person > 0) %>% pull(obs_know_person),
  num_knows_1ord = filter(analysis_data, obs_know_person > 0) %>% pull(knows_other_dewormed),
  num_knows_2ord = filter(analysis_data, obs_know_person > 0) %>% pull(thinks_other_knows),
  
  beliefs_treatment_map_design_matrix = treatment_map_design_matrix,
  
  beliefs_ate_pairs,
  num_beliefs_ate_pairs = nrow(beliefs_ate_pairs),
  
  # Take-up Model 
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
  
  cluster_treatment_map,
  
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
  use_u_in_delta = FALSE,
  multithreaded = script_options$threads > 1,
  cluster_log_lik = TRUE,
  generate_rep = FALSE,
  generate_sim = FALSE,
  fit_model_to_data = !script_options$prior,
  fit_wtp_model_to_data = !script_options$prior,
  cross_validate = script_options$cv,
  use_wtp_model = FALSE,
  use_homoskedastic_shocks = FALSE,
  
  thin = thin_by,
  
  alg_sol_f_tol = 1e-5,
  alg_sol_rel_tol = 1e-5, 
  alg_sol_max_steps = 1e6L,
  
  CALENDAR_TREATMENT_INDEX = which(fct_match(fct_unique(cluster_assigned_treatment), "calendar")),
  BRACELET_TREATMENT_INDEX = which(fct_match(fct_unique(cluster_assigned_treatment), "bracelet")),
  
  # Priors
  
  u_splines_v_sigma_sd = 1,
  mu_rep_sd = 1,
  
  dist_beta_v_sd = 0.25,
  
  beta_control_sd = 1,
  beta_far_effect_sd = 0.5,
  beta_ink_effect_sd = 0.5,
  beta_calendar_effect_sd = 0.5,
  beta_bracelet_effect_sd = 0.5,
  
  structural_beta_county_sd_sd = 0.25,
  structural_beta_cluster_sd_sd = 0.25,
  
  analysis_data
) %>% 
  list_modify(!!!map(models, pluck, "model_type") %>% set_names(~ str_c("MODEL_TYPE_", .))) %>% 
  list_modify(!!!wtp_stan_data)

# Stan Run ----------------------------------------------------------------

## Takeup ------------------------------------------------------------------
if (script_options$takeup) {
  models <- if (!is_null(script_options$models)) {
    models_to_run <- script_options$models %>% 
      map_if(str_detect(., r"{\d+}"), as.integer, .else = ~ str_which(names(models), .x)) %>% 
      unlist()
    
    models[models_to_run]
  } else {
    models
  }
  
  if (script_options$fit || script_options$prior) {
    dist_fit <- models %>% 
      stan_list(stan_data, script_options, use_cmdstanr = script_options$cmdstanr, include_paths = script_options$include_paths)
    
    if (script_options$cmdstanr) {
      dist_fit_obj <- dist_fit
      
      dist_fit %<>%
        imap(~ file.path(script_options$output_path, str_c(output_name, "_", .y, ".rds")))
     
      # BUG spaces in paths causing problems. Wait till it is fixed.
      try(iwalk(dist_fit_obj, ~ .x$save_output_files(dir = script_options$output_path, basename = str_c(output_name, .y, sep = "_"), timestamp = FALSE, random = FALSE)))
      
      dist_fit_obj %>% 
        iwalk(~ {
          cat(.y, "Diagnosis ---------------------------------\n")
          .x$cmdstan_diagnose()
        })
    }
    
    if (script_options$update_output) {
      new_dist_fit <- dist_fit
    
      dist_fit <- tryCatch({
        old_data_env <- new.env()
        load(output_file_name, envir = old_data_env)
        
        list_modify(old_data_env$dist_fit, !!!new_dist_fit)
      }, error = function(e) dist_fit)  
    }
  } else if (script_options$cv) {
    dist_kfold <- models %>% stan_list(stan_data, use_cmdstanr = script_options$cmdstanr, include_paths = script_options$include_paths)
    
    if (script_options$cmdstanr) {
      dist_kfold_obj <- dist_kfold
      
      dist_kfold %<>%
        imap(~ file.path(script_options$output_path, str_c(output_name, "_", .y, ".rds")))
    }
    
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
    if (script_options$fit || script_options$prior) {
      save(dist_fit, models, grid_dist, stan_data, file = output_file_name)
      
      if (script_options$cmdstanr) {
        walk2(dist_fit, dist_fit_obj, ~ .y$save_object(.x))
      }  
    } else if (script_options$cv) {
      save(dist_kfold, models, file = output_file_name)
      
      if (script_options$cmdstanr) {
        walk2(dist_kfold, dist_kfold_obj, ~ .y$save_object(.x))
      }  
    }
  }
## Beliefs -----------------------------------------------------------------
} else if (script_options$beliefs) {
  beliefs_model <- cmdstan_model(file.path("stan_models", "secobeliefs.stan"), include_paths = script_options$include_paths)
  
  beliefs_fit <- beliefs_model$sample(
    data = stan_data %>% discard(~ any(is.na(.x))),
    chains = script_options$chains,
    parallel_chains = script_options$chains,
    iter_warmup = script_options$iter %/% 2,
    iter_sampling = script_options$iter %/% 2,
    adapt_delta = 0.9
  )
  
  beliefs_fit$save_output_files(
    dir = script_options$output_path,
    basename = str_c(output_name, "_fit"), 
    timestamp = FALSE, random = FALSE
  )

  beliefs_fit$save_object(file.path(script_options$output_path, str_c(output_name, "_fit.rds")))

### Results -----------------------------------------------------------------

  beliefs_results <- lst(
    prob_knows = beliefs_fit$draws(c("prob_1ord", "prob_2ord")) %>% 
      posterior::as_draws_df() %>% 
      mutate(iter_id = .draw) %>% 
      pivot_longer(-c(iter_id, .draw, .chain, .iteration)) %>% 
      tidyr::extract(name, c("ord", "treatment_id"), r"{([12])ord\[(\d+)\]}", convert = TRUE) %>% 
      group_by(ord, treatment_id) %>% 
      summarize(
        per = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
        per_val = quantile(value, per),
        .groups = "drop"
      ) %>% 
      ungroup() %>% 
      pivot_wider(c(treatment_id, ord), names_from = per, values_from = per_val, names_prefix = "per_") %>% 
      right_join(mutate(cluster_treatment_map, treatment_id = seq(n())), ., by = "treatment_id") %>% 
      arrange(ord),
    
    ate_knows = beliefs_fit$draws(c("ate_1ord", "ate_2ord")) %>% 
      posterior::as_draws_df() %>% 
      mutate(iter_id = .draw) %>% 
      pivot_longer(-c(iter_id, .draw, .chain, .iteration)) %>% 
      tidyr::extract(name, c("ord", "ate_pair_index"), r"{([12])ord\[(\d+)\]}", convert = TRUE) %>% 
      group_by(ord, ate_pair_index) %>% 
      summarize(
        per = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
        per_val = quantile(value, per),
        .groups = "drop"
      ) %>% 
      ungroup() %>% 
      pivot_wider(c(ord, ate_pair_index), names_from = per, values_from = per_val, names_prefix = "per_") %>% 
      right_join(mutate(beliefs_ate_pairs, ate_pair_index = seq(n())), ., by = "ate_pair_index") %>% 
      right_join(mutate(cluster_treatment_map, treatment_id = seq(n())), ., by = c("treatment_id")) %>%
      right_join(mutate(cluster_treatment_map, treatment_id = seq(n())), ., by = c("treatment_id" = "treatment_id_control"), suffix = c("_right", "_left")) %>%
      select(-starts_with("treatment_id"), -ate_pair_index) %>% 
      arrange(ord)
  )
  
  write_rds(beliefs_results, file.path(script_options$output_path, str_c(output_name, "_results.rds")))
}
  
cat(str_glue("All done. Saved results to output ID '{output_name}'\n\n"))
