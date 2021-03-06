#!u/sr/bin/Rscript
#
# This script is used to postprocess Stan fit for the various models, reduced form and structural. In addition to putting our analysis in a format that 
# allows for easy extraction of all levels and treatment effects, it allows handles imputing take-up levels for counterfactuals using Stan-generated cost-benefits 
# and reputational returns parameters: using these we can calculate the probability of take-up after calculating the v^* fixed point solution.
#

script_options <- docopt::docopt(
  stringr::str_glue(
"Usage:
  postprocess_dist_fit.R <fit-version> [--full-outputname --cores=<num-cores> --output-path=<path> --input-path=<path> --load-from-csv --no-rate-of-change --keep-fit]
  
Options:
  --cores=<num-cores>  Number of cores to use [default: 12]
  --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
  --output-path=<path>  Path to find results [default: temp-data]
  --keep-fit "), 
  # args = if (interactive()) "29" else commandArgs(trailingOnly = TRUE)
  # args = if (interactive()) "30" else commandArgs(trailingOnly = TRUE)
  # args = if (interactive()) "test3 --full-outputname" else commandArgs(trailingOnly = TRUE)
  # args = if (interactive()) "31 --cores=6" else commandArgs(trailingOnly = TRUE) 
  # args = if (interactive()) "test --full-outputname --cores=4 --input-path=/tigress/kn6838/takeup --output-path=/tigress/kn6838/takeup" else commandargs(trailingonly = true) 
  args = if (interactive()) "51 --cores=4 --load-from-csv" else commandArgs(trailingOnly = TRUE) 
)

library(magrittr)
library(tidyverse)
library(rlang)
library(cmdstanr)
library(furrr)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

fit_version <- script_options$fit_version
postprocess_cores <- as.integer(script_options$cores)

if (interactive()) {
  plan(multisession, workers = postprocess_cores)
} else {
  plan(multicore, workers = postprocess_cores)
}

# Analysis Data --------------------------------------------------------------------

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

analysis_data <- monitored_nosms_data

# Load Data ---------------------------------------------------------------

roc_param <- c("cluster_roc_diff", 
               str_c(rep(c("cluster_roc", "cluster_rep_return", "cluster_social_multiplier", "cluster_w_cutoff", "cluster_takeup_prop"), each = 2), c("_left", "_right")))

param_used <- c(
  "total_error_sd", "u_sd", "cluster_cf_cutoff", roc_param, "sim_delta", "obs_cluster_mu_rep", # "mu_beliefs_effect", 
  "group_dist_mean", "group_dist_sd", "group_dist_mix", "missing_cluster_standard_dist", 
  "prob_prefer_calendar", "strata_wtp_mu", "hyper_wtp_mu",
  "prob_1ord", "prob_2ord", "ate_1ord", "ate_2ord"
)

# param_used %<>% c("structural_cluster_takeup_prob", "missing_cluster_standard_dist") 

load_from_csv <- function(fit_file, input_path, param) {
  dir(input_path, pattern = str_c(str_remove(basename(fit_file), fixed(".rds")), r"{-\d+.csv}"), full.names = TRUE) %>% 
      read_cmdstan_csv(variables = param) 
}

load_fit <- function(fit_file, input_path = script_options$input_path, load_from_csv = script_options$load_from_csv, param = param_used) {
  if (load_from_csv) {
    load_from_csv(fit_file, input_path, param) %>% 
      pluck("post_warmup_draws") %>% 
      posterior::as_draws_df() %>% 
      mutate(iter_id = .draw) %>% 
      pivot_longer(!c(iter_id, .draw, .iteration, .chain), names_to = "variable", values_to = "iter_est") %>% 
      nest(iter_data = !variable)
  } else {
    read_rds(file.path(input_path, basename(fit_file)))
  }
}

# Metadata on the models fit
model_info <- tribble(
  ~ model,                          ~ model_name,                                ~ model_type,
  
  "REDUCED_FORM_NO_RESTRICT",       "Reduced Form",                              "reduced form",
  "STRUCTURAL_LINEAR_U_SHOCKS",     "Structural",                                "structural",
  # "STRUCTURAL_LINEAR",              "Structural",                                "structural",
  # "STRUCTURAL_QUADRATIC",           "Structural Quadratic Cost",                 "structural",
  # "STRUCTURAL_QUADRATIC_NO_SHOCKS", "Structural Quadratic Cost (No Shocks)",     "structural",
  # "STRUCTURAL_LINEAR_NO_SHOCKS",    "Structural Linear Cost (No Shocks)",        "structural",
  # "STRUCTURAL_QUADRATIC_SALIENCE",  "Structural Quadratic Cost With Salience",   "structural",
  # "STRUCTURAL_LINEAR_SALIENCE",     "Structural With Salience",                 "structural",
  
  "STACKED",                        "Stacked Model",                            "combined", # Includes both reduced form and structural models
  "STRUCTURAL_STACKED",             "Structural Stacked Model",                 "structural", 
) %>% 
  mutate(model_type = factor(model_type, levels = c("reduced form", "structural", "combined")))

# Stan fit 
dist_fit_data <- tryCatch({
  if (script_options$full_outputname) {
    load(file.path(script_options$input_path, str_interp("${fit_version}.RData")))
  } else {
    load(file.path(script_options$input_path, str_interp("dist_fit${fit_version}.RData")))
  }
  
  dist_fit %<>% 
    map_if(is.character, ~ tryCatch(load_fit(.x), error = function(err) load_fit(.x, param = NULL))) 

  if (has_name(dist_fit, "value")) {
    dist_fit_warnings <- dist_fit$warning
    dist_fit %<>% 
      list_modify(!!!.$value, value = NULL, warning = NULL)
  }

  # Convert to tibble with list column of fit objects
  dist_fit_data <- enframe(dist_fit, name = "model", value = "fit") %>% 
    mutate(stan_data = list(stan_data))
  
  rm(dist_fit, stan_data)
  
  dist_fit_data
}, error = function(err) NULL)

# Prior predicted fit for same models
dist_fit_data <- tryCatch({
  load(file.path(script_options$input_path, str_interp("dist_prior${fit_version}.RData")))
  
  dist_fit %<>% 
    map_if(is.character, ~ tryCatch(load_fit(.x), error = function(err) load_fit(.x, param = NULL))) 
  
  dist_fit_data %<>%
    bind_rows("fit" = .,
              "prior-predict" = enframe(dist_fit, name = "model", value = "fit") %>% 
                mutate(stan_data = list(stan_data)),
              .id = "fit_type") 
  
  rm(dist_fit, stan_data)
  
  dist_fit_data
},
error = function(err) {
  dist_fit_data %>% 
    mutate(fit_type = "fit")
})

dist_fit_data %<>% 
  mutate(fit_type = factor(fit_type, levels = c("fit", "prior-predict")))
    
dist_fit_data <- tryCatch({
  # Load cross-validation data
  load(file.path(script_options$input_path, str_interp("dist_kfold${fit_version}.RData")))
  
  dist_fit_data <- enframe(dist_kfold, name = "model", value = "kfold") %>% 
    inner_join(select(model_info, model, model_type), by = "model") %>% 
    mutate(
      fit_type = factor("fit", levels = c("fit", "prior-predict")),
      stacking_weight = map(kfold, pluck, "pointwise") %>% 
        do.call(rbind, .) %>% 
        t() %>% 
        loo::stacking_weights(),
      
      kfold = set_names(kfold, model)
    ) %>% 
    group_by(model_type) %>% 
    mutate(
      stacking_weight_by_type = map(kfold, pluck, "pointwise") %>% 
        do.call(rbind, .) %>% 
        t() %>% { 
        tryCatch(loo::stacking_weights(.), error = function(err) NA)
        }
    ) %>% 
    ungroup() %>% 
    left_join(kfold_compare(x = discard(.$kfold, is_null)), by = "model") %>% 
    left_join(dist_fit_data, ., by = c("model", "fit_type")) %>% 
    select(-one_of("model_type")) %>% 
    mutate(across(where(~ is(.x, "stacking_weights")), as.numeric))

  rm(dist_kfold)
  
  return(dist_fit_data)
},
error = function(error) dist_fit_data)

dist_fit_data %<>% 
  inner_join(select(model_info, model, model_type), by = "model") 

observed_takeup <- monitored_nosms_data %>% 
  group_by(cluster_id, assigned_treatment = assigned.treatment, assigned_dist = cluster.dist.to.pot) %>% 
  summarize(prop_takeup = mean(dewormed), 
            se = sqrt(prop_takeup * (1 - prop_takeup) / n()),
            prop_takeup_ub = prop_takeup + se,
            prop_takeup_lb = prop_takeup - se, 
            .groups = "drop") 

# Functions ---------------------------------------------------------------

quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

calculate_prob_and_num_takeup <- function(cutoffs, total_error_sd, force_draw_takeup = FALSE) {
  prob_iter_data <- if (!is_null(total_error_sd)) {
    mutate(cutoffs, iter_data = future_map(iter_data, 
                               ~ left_join(.x, .y, by = "iter_id", suffix = c("", "_total_error_sd")) %>% 
                                 mutate(prob = pnorm(- iter_est, sd = iter_est_total_error_sd)), 
                              total_error_sd$iter_data[[1]]))
  } else {
    mutate(cutoffs, iter_data = map(iter_data, ~ mutate(.x, prob = pnorm(- iter_est)))) 
  }
  
  prob_iter_data %>% 
    mutate(
      iter_data = future_pmap(lst(iter_data, cluster_size, obs_num_takeup), function(iter_data, cluster_size, obs_num_takeup, force_draw_takeup) {
        if (!is.na(obs_num_takeup) && !force_draw_takeup) {
          mutate(iter_data, iter_num_takeup = obs_num_takeup)
        } else {
          mutate(iter_data, iter_num_takeup = rbinom(n(), cluster_size, prob)) # If not observed treatment, draw random number of takers for cluster
        }
      }, force_draw_takeup = force_draw_takeup, .options = furrr_options(seed = TRUE))
    )
}

# Combine cluster level nested data to produce treatment level nested data
organize_by_treatment <- function(.data, condition_on_dist, ...) {
  .data %>% {
    if (condition_on_dist) { # Only use observed distance
      filter(., assigned_dist_group_obs == assigned_dist_group)
    } else .
  } %>% 
    select(..., cluster_size, iter_data) %>% 
    mutate(iter_data = map(iter_data, select, iter_id, prob, iter_num_takeup)) %>%
    unnest(iter_data) %>% 
    group_by(iter_id, ...) %>% 
    summarize(iter_prop_takeup = sum(iter_num_takeup) / sum(cluster_size), .groups = "drop") %>% 
    nest(iter_data = c(iter_id, iter_prop_takeup)) %>%
    mutate(
      mean_est = map_dbl(iter_data, ~ mean(.$iter_prop_takeup, na.rm = TRUE)),
      takeup_quantiles = map(iter_data, quantilize_est, iter_prop_takeup, wide = TRUE, quant_probs = c(quant_probs), na.rm = TRUE)
    ) %>% 
    unnest(takeup_quantiles)
}

level_stack_reducer <- function(accum, est_takeup_level, stacking_weight) {
  if (is.na(stacking_weight)) {
    return(accum)
  } else {
    weighed_level <- est_takeup_level %>% 
      select(-mean_est, -starts_with("per_")) %>% 
      mutate(iter_data = map(iter_data, mutate, iter_prop_takeup = iter_prop_takeup * stacking_weight),
             mu_assigned_treatment = if (!is.factor(mu_assigned_treatment)) factor(mu_assigned_treatment) else mu_assigned_treatment,
             mu_assigned_treatment = if_else(is.na(mu_assigned_treatment), assigned_treatment, mu_assigned_treatment))
     
    if (is_empty(accum)) {
      return(weighed_level)
    } else {
      inner_join(accum, weighed_level, by = setdiff(intersect(names(accum), names(weighed_level)), "iter_data"), suffix = c("_accum", "_weighed")) %>% 
        mutate(iter_data = map2(iter_data_accum, iter_data_weighed, inner_join, by = "iter_id", suffix = c("_accum", "_weighed")) %>% 
                 map(mutate, iter_prop_takeup = iter_prop_takeup_accum + iter_prop_takeup_weighed) %>% 
                 map(select, -ends_with("_accum"), -ends_with("_weighed"))) %>% 
        select(-ends_with("_accum"), -ends_with("_weighed"))
    }
  }
}

rate_of_change_stack_reducer <- function(accum, rate_of_change, stacking_weight) {
  if (is.na(stacking_weight)) {
    return(accum)
  } else {
    weighed_level <- rate_of_change %>% 
      select(-c(prob_takeup, social_multiplier, partial_bbar)) %>% 
      mutate(iter_data = map(iter_data, mutate_at, vars(iter_prob_takeup, iter_social_multiplier, iter_partial_bbar), ~ . * stacking_weight))
     
    if (is_empty(accum)) {
      return(weighed_level)
    } else {
      inner_join(accum, weighed_level, by = setdiff(intersect(names(accum), names(weighed_level)), "iter_data"), suffix = c("_accum", "_weighed")) %>% 
        mutate(iter_data = map2(iter_data_accum, iter_data_weighed, inner_join, by = "iter_id", suffix = c("_accum", "_weighed")) %>% 
                 map(mutate, 
                     iter_prob_takeup = iter_prob_takeup_accum + iter_prob_takeup_weighed,
                     iter_social_multiplier = iter_social_multiplier_accum + iter_social_multiplier_weighed,
                     iter_partial_bbar = iter_partial_bbar_accum + iter_partial_bbar_weighed) %>% 
                 map(select, -ends_with("_accum"), -ends_with("_weighed"))) %>% 
        select(-ends_with("_accum"), -ends_with("_weighed"))
    }
  }
}

summarize_roc <- function(param_data) {
  unnest(param_data, iter_data) %>%
    group_by(roc_distance_index, roc_distance, iter_id) %>% 
    summarize(iter_est = mean(iter_est), .groups = "drop") %>% 
    ungroup() %>% 
    nest(iter_data = c(iter_id, iter_est)) %>% 
    mutate(
      mean_est = map_dbl(iter_data, ~ mean(.x$iter_est, na.rm = TRUE)),
      quants = map(iter_data, quantilize_est, iter_est, quant_probs = quant_probs, na.rm = TRUE)
    ) %>% 
    unnest(quants)
}

extract_sim_delta <- function(fit, stan_data) {
  temp <- fit %>% 
    filter(str_detect(variable, "sim_delta"))
  
  if (nrow(temp) > 0) {
    temp %>% 
      transmute(w = stan_data$sim_delta_w, iter_data) %>% 
      mutate(
        mean_est = map_dbl(iter_data, ~ mean(.x$iter_est, na.rm = TRUE)),
        quants = map(iter_data, quantilize_est, iter_est, quant_probs = quant_probs, na.rm = TRUE)
      ) %>% 
      unnest(quants)
  } else return(NULL)
}

# Postprocessing ----------------------------------------------------------

dist_fit_data %<>% 
  mutate(
    total_error_sd = map2(fit, stan_data, ~ extract_obs_fit_level(.x, par = "total_error_sd", stan_data = .y, iter_level = "none", quant_probs = quant_probs)),
    u_sd = map2(fit, stan_data, ~ extract_obs_fit_level(.x, par = "u_sd", stan_data = .y, iter_level = "none", quant_probs = quant_probs)),
    
    obs_cluster_mu_rep = map(fit, filter, str_detect(variable, "obs_cluster_mu_rep")),
    
    cluster_cf_cutoff = pmap(lst(fit, stan_data, model_type), extract_obs_cluster_cutoff_cf, quant_probs = quant_probs) %>% 
      lst(cutoffs = ., total_error_sd, force_draw_takeup = fct_match(fit_type, "prior-predict")) %>%  
      pmap(calculate_prob_and_num_takeup), 
    
    est_takeup_level = list(cluster_cf_cutoff, FALSE) %>% # map_lgl(model_type, fct_match, "structural")) %>%  
      pmap(organize_by_treatment, mu_assigned_treatment, assigned_treatment, assigned_dist_group) %>% 
      map2(cluster_cf_cutoff,
           ~ filter(.y, assigned_dist_group_obs == assigned_dist_group) %>%
             organize_by_treatment(condition_on_dist = TRUE, mu_assigned_treatment, assigned_treatment) %>%
             bind_rows(.x, .)),
    
    obs_cluster_takeup_level = map(cluster_cf_cutoff, filter, !is.na(obs_prop_takeup)) %>% 
      map(mutate, quants = map(iter_data, quantilize_est, prob, quant_probs = quant_probs, na.rm = TRUE)) %>% 
      map(select, !iter_data) %>% 
      map(unnest, quants),
    
    group_dist_param = pmap(lst(fit, stan_data, model_type), get_dist_results),
    # imputed_dist = pmap(lst(fit, stan_data, model_type), get_imputed_dist),
    beliefs_results = map2(fit, stan_data, get_beliefs_results),
    wtp_results = map(fit, get_wtp_results),
  )

if (!script_options$no_rate_of_change) {
  dist_fit_data %<>% 
    mutate(
      map(roc_param, ~ map2(..2, ..3, extract_roc_param, ..1), fit, stan_data) %>% 
        set_names(roc_param) %>% 
        map(map_if, ~ !is_null(.x), summarize_roc) %>% 
        as_tibble() %>% 
        rename_with(str_replace_all, c("_left" = "_bracelet", "_right" = "_control"), .cols = everything()),
      
      sim_delta = map2(fit, stan_data, extract_sim_delta)
    )
}

dist_fit_data %<>% 
  select(!c(ends_with("_dist_cost"), cluster_cf_cutoff, any_of("structural_cluster_benefit")))

if (!script_options$keep_fit) {
  dist_fit_data %<>% 
    select(!fit)
}

## Combine models using stacking -------------------------------------------

if (has_name(dist_fit_data, "stacking_weight")) {
  # The two functions below are used to combined multiple takeup levels, from different models, into a single (stacked) takeup level.
  
  dist_fit_data %<>% 
    add_row(
      model = factor("STACKED"), 
      model_type = factor("combined"),
      fit_type = factor("fit", levels = c("fit", "prior-predict")),
      est_takeup_level = list(
        reduce2(.$est_takeup_level, .$stacking_weight, level_stack_reducer, .init = tibble()) %>% 
          mutate(
            mean_est = map_dbl(iter_data, ~ mean(.$iter_prop_takeup)),
            takeup_quantiles = map(iter_data, quantilize_est, iter_prop_takeup, wide = TRUE, quant_probs = c(quant_probs))
          ) %>% 
          unnest(takeup_quantiles)
      ),
    ) %>%  
    add_row(
      model = factor("STRUCTURAL_STACKED"), 
      model_type = factor("structural"),
      fit_type = factor("fit", levels = c("fit", "prior-predict")),
      est_takeup_level = list(
        reduce2(.$est_takeup_level, .$stacking_weight_by_type, level_stack_reducer, .init = tibble()) %>% 
          mutate(
            mean_est = map_dbl(iter_data, ~ mean(.$iter_prop_takeup)),
            takeup_quantiles = map(iter_data, quantilize_est, iter_prop_takeup, wide = TRUE, quant_probs = c(quant_probs))
          ) %>% 
          unnest(takeup_quantiles)
      ),
    ) 
}

dist_fit_data %<>%
  inner_join(select(model_info, model, model_name), by = "model") %>% 
  mutate(model = factor(model, levels = model_info$model)) %>% 
  arrange(model)


## ATE ---------------------------------------------------------------------

# Select all the estimate differences that need to be calculated
ate_combo <- dist_fit_data %>% 
  # filter(fct_match(model_type, "structural")) %>% 
  slice(1) %$% 
  select(est_takeup_level[[1]], mu_assigned_treatment:assigned_dist_group) %>% {
    bind_cols(rename_all(., str_c, "_left"), rename_all(., str_c, "_right"))
  } %>% {
    bind_rows(
      expand(., crossing(!!!syms(names(.)))) %>%
        filter(mu_assigned_treatment_left != mu_assigned_treatment_right |
                 assigned_treatment_left != assigned_treatment_right |
                 (fct_match(assigned_dist_group_left, "far") & fct_match(assigned_dist_group_right, "close"))),
      select(., -contains("dist_group")) %>% 
        expand(crossing(!!!syms(names(.)))) %>%
        filter(mu_assigned_treatment_left != mu_assigned_treatment_right |
                 assigned_treatment_left != assigned_treatment_right) 
    )
  }

dist_fit_data %<>% 
  mutate(
    # Calculate treatment effects
    est_takeup_te = est_takeup_level %>% 
      map(select, mu_assigned_treatment:assigned_dist_group, iter_data) %>%
      map(select_if, ~ !all(is.na(.))) %>% # Get rid of the NA mu_assigned_treatment in the reduced form results 
      map(function(level_data, ate_combo) {
        present_col <- intersect(names(level_data), c("mu_assigned_treatment", "assigned_treatment", "assigned_dist_group"))
        
        left_data <- inner_join(select(ate_combo, str_c(rep(present_col, each = 2), c("_left", "_right"))) %>% 
                                  distinct_all(.keep_all = TRUE),
                                level_data,
                                by = present_col %>% set_names(str_c(.,"_left")))
        
        inner_join(left_data,
                   level_data,
                   by = present_col %>% set_names(str_c(., "_right")),
                   suffix = c("_left", "_right"))
      },
      ate_combo = ate_combo) %>% 
    map_if(
      fct_match(model_type, "structural"),
      ~ filter(., 
        mu_assigned_treatment_left != mu_assigned_treatment_right |
        assigned_treatment_left != assigned_treatment_right |
        (is.na(assigned_dist_group_left) & is.na(assigned_dist_group_right)) |
        assigned_dist_group_left != assigned_dist_group_right),
      .else = ~ filter(., 
        assigned_treatment_left != assigned_treatment_right |
        (is.na(assigned_dist_group_left) & is.na(assigned_dist_group_right)) |
        assigned_dist_group_left != assigned_dist_group_right)) %>% 
    map(mutate, iter_data = map2(iter_data_left, iter_data_right, inner_join, by = "iter_id", suffix = c("_left", "_right")) %>% 
          map(mutate, iter_takeup_te = iter_prop_takeup_left - iter_prop_takeup_right)) %>% 
    map(select, -iter_data_left, -iter_data_right) %>% 
    map(mutate, 
        mean_est = map_dbl(iter_data, ~ mean(.$iter_takeup_te, na.rm = TRUE)),
        takeup_te_quantiles = map(iter_data, quantilize_est, iter_takeup_te, wide = TRUE, quant_probs = c(quant_probs), na.rm = TRUE)) %>% 
    map(unnest, takeup_te_quantiles),
   
    # Calculate treatment effects based on distance 
    est_takeup_dist_te = est_takeup_te %>%
      map(filter, !is.na(assigned_dist_group_left), !is.na(assigned_dist_group_right)) %>% 
      map_if(fct_match(model_type, "structural"), filter, mu_assigned_treatment_left == assigned_treatment_left, mu_assigned_treatment_right == assigned_treatment_right) %>% 
      map_if(fct_match(model_type, "structural"), select, -starts_with("mu_assigned_treatment")) %>% 
      map(filter, assigned_treatment_left == assigned_treatment_right, fct_match(assigned_dist_group_left, "far"), fct_match(assigned_dist_group_right, "close")) %>% 
      map(mutate, assigned_treatment_right = "control") %>% 
      map(select, -mean_est, -starts_with("per_"), -starts_with("assigned_dist_group")) %>% 
      map(~ inner_join(., select(., -assigned_treatment_right), by = c("assigned_treatment_right" = "assigned_treatment_left"), suffix = c("_left", "_right"))) %>% 
      map(filter, assigned_treatment_left != assigned_treatment_right) %>% 
      map(mutate, iter_data = map2(iter_data_left, iter_data_right, inner_join, by = "iter_id", suffix = c("_left", "_right")) %>% 
            map(transmute, iter_id, iter_takeup_dist_te = iter_takeup_te_left - iter_takeup_te_right)) %>% 
      map(select, -iter_data_left, -iter_data_right) %>% 
      map(mutate, 
          mean_est = map_dbl(iter_data, ~ mean(.$iter_takeup_dist_te)),
          takeup_te_quantiles = map(iter_data, quantilize_est, iter_takeup_dist_te, wide = TRUE, quant_probs = quant_probs, na.rm = TRUE)) %>% 
      map(unnest, takeup_te_quantiles),
    
      est_takeup_level = map(est_takeup_level, mutate, iter_data = map(iter_data, as_tibble))  
  )


save(dist_fit_data, file = file.path(script_options$output_path, str_interp("processed_dist_fit${fit_version}.RData")))

dist_fit_data %<>% 
  mutate(across(where(is_list), map_if, ~ is_tibble(.x) && has_name(.x, "iter_data"), ~ select(.x, !iter_data)))

save(dist_fit_data, file = file.path(script_options$output_path, str_interp("processed_dist_fit${fit_version}_lite.RData")))

plan(sequential)

cat(str_glue("Post processing completed [version {fit_version}]\n\n"))
