#!/usr/bin/Rscript

script_options <- docopt::docopt(
"Usage:
  run_stan_dist_sim [options]
  
Options:
  --num-sim=<num>, -n <num>  Number of simulations to test [default: 100]
  --num-cores=<num-cores>  Number of cores to use [default: 12]
  --chains=<chains>  Number of Stan chains [default: 4]
  --cmdstanr  Use cmdstan.
  --include-paths=<paths>  Includes path for cmdstanr [default: .]
  --sim-iter=<iter>  Number of iterations for simulation [default: 4000]
  --rf-dgp  Reduced form dgp
  --outputname=<name>  Name to use for output
  --output-path=<path>  Where to save output files [default: temp-data]
  --load-fit  Don't run simulation but instead load it and start postprocessing.
  --no-progress-bar
",

  # args = if (interactive()) "--outputname=test --cmdstanr --include-paths=~/Code/takeup/stan_models -n 6 --num-cores=12 --rf-dgp" else commandArgs(trailingOnly = TRUE)
  args = if (interactive()) "--outputname=test --cmdstanr --include-paths=~/Code/takeup/stan_models -n 3 --num-cores=12 --sim-iter=400" else commandArgs(trailingOnly = TRUE)
) 

library(magrittr)
library(tidyverse)
library(parallel)
library(pbmcapply)
library(HRW)
library(splines2)
library(furrr)

script_options %<>% 
  modify_at(c("chains", "num_cores", "num_sim", "sim_iter"), as.integer)

options(mc.cores = script_options$num_cores)

num_workers <- script_options$num_cores %/% script_options$chains

if (num_workers > 1) {
  if (interactive()) {
    plan(multisession, workers = num_workers)
  } else {
    plan(multicore, workers = num_workers)
  }
}

if (script_options$cmdstanr) {
  library(cmdstanr)
} else {
  library(rstan)
  
  rstan_options(auto_write = TRUE)
}


output_name <- script_options$outputname
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
  arrange(cluster_id)

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

# Run ------------------------------------------------------------------

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
    use_shifting_v_dist = FALSE,
    use_u_in_delta = TRUE,
    suppress_reputation = FALSE,
    suppress_shocks = FALSE,
    generate_sim = FALSE,
    thin = 1,
    alg_sol_f_tol = 0.001,
    alg_sol_max_steps = 1e9L,
    alg_sol_rel_tol = 0.0000001,

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
    use_cluster_effects = FALSE,
    use_county_effects = FALSE,
    use_param_dist_cluster_effects = FALSE,
    use_param_dist_county_effects = FALSE,
    use_mu_cluster_effects = FALSE,
    use_mu_county_effects = FALSE,
    use_shifting_v_dist = FALSE,
    suppress_reputation = TRUE,
    suppress_shocks = TRUE,
    
    structural_beta_county_sd_sd = 1,
    structural_beta_cluster_sd_sd = 1,
    
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
)

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
  use_u_in_delta = FALSE,
  cluster_log_lik = TRUE,
  generate_rep = TRUE,
  generate_sim = FALSE,
  fit_model_to_data = FALSE,
  cross_validate = FALSE,
  
  thin = thin_by,
  
  alg_sol_f_tol = 1e-5,
  alg_sol_rel_tol = 1e-5, 
  alg_sol_max_steps = 1e6L,
  
  # Priors
  
  u_splines_v_sigma_sd = 1,
  mu_rep_sd = 1,
  
  beta_control_sd = 0.5,
  beta_far_effect_sd = 0.25,
  beta_ink_effect_sd = 0.25,
  beta_calendar_effect_sd = 0.25,
  beta_bracelet_effect_sd = 0.25,
  
  dist_beta_v_sd = 0.25,
  
  structural_beta_county_sd_sd = 0.125,
  structural_beta_cluster_sd_sd = 0.125,
  
  analysis_data
) %>% 
  list_modify(!!!map(models, pluck, "model_type") %>% set_names(~ str_c("MODEL_TYPE_", .)))

prior_stan_data <- stan_data %>% 
  list_modify(!!!models$STRUCTURAL_LINEAR_U_SHOCKS) 

rf_prior_stan_data <- stan_data %>% 
  list_modify(!!!models$REDUCED_FORM_NO_RESTRICT) 

dist_sim <- if (script_options$rf_dgp) { 
  rf_prior_stan_data %>%  
    map_if(is.factor, as.integer) %>% 
    fit_model(script_options$chains, iter = 400, script_options$cmdstanr, script_options$include_paths)
} else {
  prior_stan_data %>%  
    map_if(is.factor, as.integer) %>% 
    fit_model(script_options$chains, iter = 400, script_options$cmdstanr, script_options$include_paths)
}
  # fit_model(1, iter = 10, script_options$cmdstanr, script_options$include_paths)

# Extracting simulation data ----------------------------------------------

ate <- tribble(
  ~ assigned_treatment_left, ~ mu_assigned_treatment_left, ~ assigned_treatment_right, ~ mu_assigned_treatment_right,
  
  "bracelet",                "bracelet",                   "control",                  "control",
  "bracelet",                "control",                    "control",                  "control",
)

if (!script_options$rf_dgp) {
  total_error_sd <- extract_obs_fit_level(dist_sim, par = "total_error_sd", stan_data = stan_data, iter_level = "none") %>%
    unnest(iter_data) %>%
    rename(total_error_sd = iter_est)
  
  mu_rep <- extract_obs_fit_level(dist_sim, par = "mu_rep", stan_data = stan_data, iter_level = "none", by_treatment = TRUE, mix = FALSE, summarize_est = FALSE) %>%
    unnest(iter_data) %>%
    rename(true_mu_rep = iter_est)
}

cutoff_cf <- extract_obs_cluster_cutoff_cf(dist_sim, stan_data = stan_data) %>% 
  unnest(iter_data) %>% 
  rename(cluster_cutoff = iter_est) %>% 
  mutate(
    mu_assigned_treatment = if (script_options$rf_dgp) coalesce(mu_assigned_treatment, assigned_treatment) else mu_assigned_treatment,
    observed = treatment_index == treatment_index_obs & mu_assigned_treatment == assigned_treatment
  )

sim_data <- cutoff_cf %>% 
  filter(!is.nan(cluster_cutoff)) %>% {
    if (script_options$rf_dgp) {
      mutate(., total_error_sd = 1)
    } else {
      left_join(., total_error_sd, by = "iter_id") 
    }
  } %>%
  # left_join( # Hardcoded prob
  #   tibble(
  #     treatment_index = 1:8,
  #     prob_takeup = c(0.75, 0.45, rep(0.75, 6)) # seq(0.4, by = 0.05, length.out = 8)
  #     # prob_takeup = 0.75 
  #   ),
  #   by = "treatment_index"
  # ) %>% 
  mutate(
    prob_takeup = 1 - pnorm(cluster_cutoff, sd = total_error_sd),
    num_takeup = rbinom(n(), cluster_size, prob_takeup),
    obs_takeup = map2(num_takeup, cluster_size, ~ sample(c(rep(1, .x), rep(0, .y - .x))))
  ) %>%
  group_nest(iter_id, .key = "cluster_data") %>% 
  sample_n(script_options$num_sim) %>% 
  mutate(
    takeup_levels = map(cluster_data, select, cluster_id, assigned_treatment, mu_assigned_treatment, assigned_dist_group, cluster_size, prob_takeup) %>% 
      map(~ {
        group_by(.x, across(c(contains("assigned_treatment"), assigned_dist_group))) %>% 
          summarize(true_takeup_level = weighted.mean(prob_takeup, cluster_size)) %>% 
          ungroup()
      }),
    
    treatment_effects = map(cluster_data, select, cluster_id, assigned_treatment, mu_assigned_treatment, assigned_dist_group, cluster_size, prob_takeup) %>%  
      map(
        ~ inner_join(ate, .x, by = c("assigned_treatment_left" = "assigned_treatment", 
                                     "mu_assigned_treatment_left" = "mu_assigned_treatment")) %>%  
          left_join(.x, by = c("cluster_id", "assigned_dist_group", "cluster_size", 
                               "assigned_treatment_right" = "assigned_treatment", 
                               "mu_assigned_treatment_right" = "mu_assigned_treatment"), suffix = c("_left", "_right")) %>% 
          mutate(treatment_effect = prob_takeup_left - prob_takeup_right) %>% 
          group_by(across(c(contains("assigned_treatment"), assigned_dist_group))) %>% 
          summarize(treatment_effect = weighted.mean(treatment_effect, cluster_size)) %>% 
          ungroup()
      )
    )

if (!script_options$rf_dgp) { 
  sim_data %<>% 
    nest_join(select(mu_rep, iter_id, treatment_index, true_mu_rep), name = "true_mu_rep", by = c("iter_id"))
}

# Run simulation ----------------------------------------------------------

sim_stan_data <- prior_stan_data %>% 
  list_modify(
    beta_control_sd = 1,
    beta_far_effect_sd = 0.25,
    beta_ink_effect_sd = 0.25,
    beta_calendar_effect_sd = 0.25,
    beta_bracelet_effect_sd = 0.25,
  )

rf_sim_stan_data <- rf_prior_stan_data %>%  
  list_modify(
    beta_control_sd = 1,
    beta_far_effect_sd = 1,
    beta_ink_effect_sd = 1,
    beta_calendar_effect_sd = 1,
    beta_bracelet_effect_sd = 1,
  )
  
sim_data %<>% 
  mutate(
    sim_id = seq(n()),
    sim_fit_file = str_glue("{script_options$output_path}/{str_c('dist_sim', script_options$outputname, sep = '_')}_fit_{sim_id}.rds"),
    sim_rf_fit_file = str_glue("{script_options$output_path}/{str_c('dist_sim', script_options$outputname, sep = '_')}_rf_fit_{sim_id}.rds")
  )

write_rds(sim_data, file = str_glue("{script_options$output_path}/{str_c('dist_sim_pre', script_options$outputname, sep = '_')}.rds"))

sim_fit_model <- function(data, fit_file, sim_stan_data, iter, save = TRUE) {
  browser()
  fit <- map_if(sim_stan_data, is.factor, as.integer) %>% 
    list_modify( 
      fit_model_to_data = TRUE,
      
      takeup = data %>% 
        filter(observed) %>%
        arrange(cluster_id) %>% 
        select(obs_takeup) %>% 
        unnest(obs_takeup) %>% 
        pull(obs_takeup)
    ) %>%
    fit_model(script_options$chains, iter = iter, script_options$cmdstanr, script_options$include_paths)
  
  if (save) {
    fit$save_object(fit_file)
  }
  
  return(fit)
}
  
if (script_options$load_fit) { 
  cat("\n\nLoad fit ...")
  
  sim_data %<>% 
    mutate(
      sim_fit = map(sim_fit_file, read_rds)
    )
  
  cat("done.\n")
} else {
  cat("\n\nStarting simulations...")
  
  sim_data %<>% 
    mutate(
      # sim_rf_fit = map2(cluster_data, sim_rf_fit_file, sim_fit_model, sim_stan_data = rf_sim_stan_data, iter = script_options$sim_iter),
      sim_rf_fit = future_map2(cluster_data, sim_rf_fit_file, sim_fit_model, sim_stan_data = rf_sim_stan_data, iter = script_options$sim_iter,
                               .options = furrr_options(seed = TRUE), .progress = !script_options$no_progress_bar),
    )
  
  if (!script_options$rf_dgp) {
    sim_data %<>% 
      mutate(
        # sim_fit = map2(cluster_data, sim_fit_file, sim_fit_model, sim_stan_data = sim_stan_data, iter = script_options$sim_iter) 
        sim_fit = future_map2(cluster_data, sim_fit_file, sim_fit_model, sim_stan_data = sim_stan_data, iter = script_options$sim_iter,
                              .options = furrr_options(seed = TRUE), .progress = !script_options$no_progress_bar),
      )
  }
  
  cat("done.\n") 
}

cat("Adding prior predicted fit...")

if (script_options$rf_dgp) {
  sim_data %<>% 
    add_row(
      iter_id = 0, 
      sim_rf_fit = list(fit_model(rf_sim_stan_data, script_options$chains, iter = 1000, script_options$cmdstanr, script_options$include_paths))
    )
} else {
  sim_data %<>% 
    add_row(
      iter_id = 0, 
      sim_fit = list(fit_model(sim_stan_data, script_options$chains, iter = 1000, script_options$cmdstanr, script_options$include_paths))
    )
}

cat("done.\n")

cat("Post-processing fit...")

calculate_sim_takeup_levels <- function(fit, true_levels, stan_data) {
  total_error_sd <- extract_obs_fit_level(fit, par = "total_error_sd", stan_data = stan_data, iter_level = "none")
  
  if (!is_null(total_error_sd)) {
    total_error_sd %<>% 
      unnest(iter_data) %>% 
      rename(total_error_sd = iter_est)
  }

  curr_levels <- extract_obs_cluster_cutoff_cf(fit, stan_data) %>% 
    mutate(observed = treatment_index == treatment_index_obs & mu_assigned_treatment == assigned_treatment) %>%  
    unnest(iter_data) %>% 
    rename(cluster_cutoff = iter_est) %>% { 
      if (!is_null(total_error_sd)) {
        left_join(., total_error_sd, by = "iter_id") 
      } else {
        mutate(., total_error_sd = 1)
      }
    } %>% 
    mutate(
      prob_takeup = 1 - pnorm(cluster_cutoff, sd = total_error_sd),
      num_takeup = if_else(observed, obs_num_takeup, rbinom(n(), cluster_size, prob_takeup)),
    ) %>% 
    group_by(iter_id, assigned_treatment, mu_assigned_treatment, assigned_dist_group) %>% 
    summarize(sim_takeup_level = weighted.mean(prob_takeup, cluster_size, na.rm = TRUE)) %>% 
    nest(iter_data = c(iter_id, sim_takeup_level)) %>% 
    mutate(mean_est = map_dbl(iter_data, ~ mean(.$sim_takeup_level, na.rm = TRUE)),
           quantiles_est = map(iter_data, quantilize_est, 
                               sim_takeup_level,
                               quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95),
                               na.rm = TRUE)) %>%
    unnest(quantiles_est)
  
  if (!is_null(true_levels)) {
    curr_levels %>% 
      mutate(mu_assigned_treatment = coalesce(mu_assigned_treatment, assigned_treatment)) %>% 
      left_join(true_levels, by = c("assigned_treatment", "mu_assigned_treatment", "assigned_dist_group")) %>% 
      mutate(in_0.8_ci = true_takeup_level <= per_0.9 & true_takeup_level >= per_0.1) 
  } else {
    return(curr_levels)
  }
}

sim_data %<>% 
  mutate(
    rf_sim_takeup_levels = future_map2(sim_rf_fit, takeup_levels, calculate_sim_takeup_levels, rf_prior_stan_data, .options = furrr_options(seed = TRUE), .progress = !script_options$no_progress_bar),
  )

if (!script_options$rf_dgp) {
  sim_data %<>% 
    mutate(
      sim_takeup_levels = future_map2(sim_fit, takeup_levels, calculate_sim_takeup_levels, prior_stan_data, .options = furrr_options(seed = TRUE), .progress = !script_options$no_progress_bar),

      sim_mu_rep = future_map2(sim_fit, true_mu_rep, ~ {
        curr_mu_rep <- extract_obs_fit_level(.x, par = "mu_rep", stan_data = stan_data, iter_level = "none", by_treatment = TRUE, mix = FALSE, summarize_est = TRUE) %>%
          unnest(quantiles_est)

        if (!is_null(.y)) {
          curr_mu_rep %>%
            left_join(.y, by = c("treatment_index")) %>%
            mutate(in_0.8_ci = true_mu_rep <= per_0.9 & true_mu_rep >= per_0.1)
        } else {
          return(curr_mu_rep)
        }
      }, .progress = !script_options$no_progress_bar),
    )
}

cat("done.\n")

plan(sequential)

cat("Saving...")

write_rds(sim_data, file = str_glue("{script_options$output_path}/{str_c('dist_sim', script_options$outputname, sep = '_')}.rds"))

cat("done.\n")

sim_data$sim_fit %>% 
  walk(~ .x$cmdstan_diagnose())
