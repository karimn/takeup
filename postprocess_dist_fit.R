#!/usr/bin/Rscript
#
# This script is used to postprocess Stan fit for the various models, reduced form and structural. In addition to putting our analysis in a format that 
# allows for easy extraction of all levels and treatment effects, it allows handles imputing take-up levels for counterfactuals using Stan-generated cost-benefits 
# and reputational returns parameters: using these we can calculate the probability of take-up after calculating the v^* fixed point solution.
#

script_options <- docopt::docopt(
"Usage:
  postprocess_dist_fit <fit-version>",

  args = if (interactive()) "28" else commandArgs(trailingOnly = TRUE) 
) 

library(magrittr)
library(tidyverse)
library(rlang)
library(rstan)
library(pbmcapply)
library(nleqslv)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

fit_version <- script_options$`fit-version`


# Load Data ---------------------------------------------------------------

# Stan fit 
load(file.path("data", "stan_analysis_data", str_interp("dist_fit${fit_version}.RData")))

# Convert to tibble with list column of fit objects
dist_fit_data <- enframe(dist_fit, name = "model", value = "fit")

rm(dist_fit)

# Metadata on the models fit
model_info <- tribble(
  ~ model,                         ~ model_name,                                ~ model_type,
  
  "REDUCED_FORM_NO_RESTRICT",      "Reduced Form Discrete Cost",                "reduced form",
  "STRUCTURAL_QUADRATIC",          "Structural Quadratic Cost",                 "structural",
  "STRUCTURAL_QUADRATIC_SALIENCE", "Structural Quadratic Cost With Salience",   "structural",
  
  "STACKED",                        "Stacked Model",                            "combined", # Includes both reduced form and structural models
  "STRUCTURAL_STACKED",             "Structural Stacked Model",                 "structural", 
) %>% 
  mutate(model_type = factor(model_type, levels = c("reduced form", "structural", "combined")))

# How much to thin the posterior samples, mostly so the postprocessing fits in memory
thin_model <- c(
  "REDUCED_FORM_NO_RESTRICT" = 1L,
  "REDUCED_FORM_SEMIPARAM" = 1L,
  "STRUCTURAL_LINEAR" = 1L,
  "STRUCTURAL_QUADRATIC" = 1L,
  "STRUCTURAL_QUADRATIC_SALIENCE" = 4L
)

dist_fit_data %<>% 
  left_join(enframe(thin_model, name = "model", value = "thin"), by = "model") 

# Prior predicted fit for same models
dist_fit_data <- tryCatch({
  load(file.path("data", "stan_analysis_data", str_interp("dist_fit_prior${fit_version}.RData")))
  
  dist_fit_data %<>%
    bind_rows("fit" = .,
              "prior-predict" = enframe(dist_fit, name = "model", value = "fit") %>% 
                mutate(thin = 1L),
              .id = "fit_type") 
  
  rm(dist_fit)
  
  dist_fit_data
},
error = function(err) {
  dist_fit_data %>% 
    mutate(fit_type = "fit")
})

dist_fit_data %<>% 
  mutate(fit_type = factor(fit_type, levels = c("fit", "prior-predict")))

# Load cross-validation data
load(file.path("data", "stan_analysis_data", str_interp("dist_kfold${fit_version}.RData")))
    
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
  select(-one_of("model_type"))

rm(dist_kfold)

dist_fit_data %<>% 
  inner_join(select(model_info, model, model_type), by = "model") %>% 
  mutate(thin = coalesce(thin, 1L))

observed_takeup <- monitored_nosms_data %>% 
  group_by(cluster_id, assigned_treatment = assigned.treatment, assigned_dist = cluster.dist.to.pot) %>% 
  summarize(prop_takeup = mean(dewormed), 
            se = sqrt(prop_takeup * (1 - prop_takeup) / n()),
            prop_takeup_ub = prop_takeup + se,
            prop_takeup_lb = prop_takeup - se) %>% 
  ungroup()


# Functions ---------------------------------------------------------------

quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
postprocess_cores <- 12

# Join benefit-cost, reputation mu, and total error sd
join_benefit_cost_mu_error_sd <- function(benefit_cost, mu_rep, total_error_sd) {
  if (!is_null(mu_rep)) {
    mutate(benefit_cost, 
           iter_data = pbmclapply(mc.cores = postprocess_cores,
                                  iter_data, 
                                  left_join, 
                                  unnest(mu_rep, iter_data) %>% 
                                    select(mu_assigned_treatment = assigned_treatment, iter_id, iter_est) %>%  
                                    inner_join(total_error_sd, by = c("iter_id"), suffix = c("_mu", "_error_sd")),
                                  by = c("iter_id")) %>% 
             map(mutate, obs_num_takeup = if_else(assigned_treatment == mu_assigned_treatment, obs_num_takeup, NA_integer_)))
  } else if (!is_null(total_error_sd)) { 
    mutate(benefit_cost, 
           iter_data = pbmclapply(mc.cores = postprocess_cores,
                                  iter_data, 
                                  inner_join, 
                                  total_error_sd, 
                                  by = c("iter_id"), 
                                  suffix = c("", "_error_sd")) %>% 
             map(mutate, mu_assigned_treatment = NA_character_, iter_est_mu = 0))
  } else {
    mutate(benefit_cost,
           iter_data = map(iter_data, mutate, iter_est_error_sd = 1, iter_est_mu = 0) %>% 
             map(mutate, mu_assigned_treatment = NA_character_, iter_est_mu = 0))
  } 
}

# Given multiple variables, generate quantiles and pack them
prep_multiple_est <- function(accum_data, next_col) {
  group_name <- str_remove(as.character(next_col), "^iter_")
                                           
  mutate(accum_data,
         mean_est = map_dbl(iter_data, ~ mean(pull(., !!next_col))), 
         quants = map(iter_data, quantilize_est, !!next_col, wide = TRUE, quant_probs = c(quant_probs))) %>% 
    unnest(quants) %>% 
    pack({{ group_name }} := c(mean_est, starts_with("per_"))) 
}

# For a range of v^* calculate the social multiplier effect as well as the partial differentiation of E[Y] wrt to \bar{B} 
simulate_social_multiplier <- function(mu_rep, total_error_sd, fit_type) {
  multiplier_v_range <- seq(-2, 2, 0.25) # Simulated values for v^*, used calculate reputational multiplier effect 
  
  if (fct_match(fit_type, "fit") && !is_null(mu_rep) && !is_null(total_error_sd)) {
    mu_sd <- mu_rep %>% 
      select(assigned_treatment, iter_data) %>% 
      unnest(iter_data) %>% 
      select(assigned_treatment, iter_id, iter_est) %>% 
      left_join(total_error_sd %>% 
                  select(iter_data) %>% 
                  unnest(iter_data),
                by = "iter_id", suffix = c("_mu", "_sd"))
    
    multiplier_v_range %>% 
      map_df(~ mu_sd %>% 
               mutate(
                 v = .x,
                 iter_prob_takeup = pnorm(v, sd = iter_est_sd, lower.tail = FALSE),  
                 iter_social_multiplier = social_multiplier(v, iter_est_mu),
                 iter_partial_bbar = expect_y_partial_bbar(v, iter_est_mu, iter_est_sd)
               ) %>% 
      nest(iter_data = c(iter_id, iter_est_mu, iter_est_sd, iter_prob_takeup, iter_social_multiplier, iter_partial_bbar))) %>% 
      reduce(exprs(iter_prob_takeup, iter_social_multiplier, iter_partial_bbar), prep_multiple_est, .init = .)
  } 
}

# Combine cluster level nested data to produce treatment level nested data
organize_by_treatment <- function(.data, ...) {
  .data %>%
    select(
      cluster_id, 
      assigned_treatment, assigned_dist_group, mu_assigned_treatment, 
      assigned_dist_obs, assigned_treatment_obs, assigned_dist_group_obs, cluster_size, 
      iter_data
    ) %>% 
    mutate(iter_data = map(iter_data, select, iter_id, prob, iter_num_takeup)) %>%
    unnest(iter_data) %>% 
    group_by(iter_id, !!!exprs(...)) %>% 
    summarize(iter_prop_takeup = sum(iter_num_takeup) / sum(cluster_size)) %>% 
    ungroup() %>% 
    nest(iter_data = c(iter_id, iter_prop_takeup)) %>%
    mutate(
      mean_est = map_dbl(iter_data, ~ mean(.$iter_prop_takeup)),
      takeup_quantiles = map(iter_data, quantilize_est, iter_prop_takeup, wide = TRUE, quant_probs = c(quant_probs))
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

# Postprocessing ----------------------------------------------------------

dist_fit_data %<>% 
  mutate(
    total_error_sd = map(fit, extract_obs_fit_level, par = "total_error_sd", stan_data = stan_data, iter_level = "none", quant_probs = quant_probs),
    
    mu_rep = map(fit, extract_obs_fit_level, par = "mu_rep", stan_data = stan_data, iter_level = "treatment", mix = FALSE, quant_probs = quant_probs),
    
    cluster_cf_benefit_cost = map2(fit, thin, ~ extract_obs_cf(.x, par = "cluster_cf_benefit_cost", stan_data = stan_data, iter_level = "cluster", quant_probs = quant_probs, thin = .y)) %>% 
      map(nest, iter_data = c(treatment_index, obs_num_takeup, matches("^assigned_(treatment|dist_group)$"), iter_data)) %>% 
      map(mutate, iter_data = map(iter_data, unnest, iter_data)) %>% 
      lst(benefit_cost = ., 
          mu_rep, 
          total_error_sd = map_if(total_error_sd, ~ !is_null(.x), select, -rhat, -starts_with("ess")) %>% 
            map_if(~ !is_null(.x), unnest, iter_data)) %>%
      pmap(join_benefit_cost_mu_error_sd) %>% 
      map(
        mutate,
        iter_data = pbmclapply(mc.cores = postprocess_cores,
                               iter_data,
                               # In this function solve for the fixed point equilibrium level v^*, given current benefit-cost and mu
                               function(curr_iter_data) {
                                 mutate(
                                   curr_iter_data,
                                   no_rep_prob = pnorm(iter_est, sd = iter_est_error_sd), # Probability of take-up in cluster without reputation
                                   v_cutoff = map2_dbl(iter_est, iter_est_mu, ~ nleqslv(x = - ..1, fn = generate_v_cutoff_fixedpoint(..1, ..2)) %>% pluck("x")), # v^* solution
                                   delta_rep = rep_normal(v_cutoff), 
                                   prob = pnorm(- v_cutoff, sd = iter_est_error_sd) # Probability of take-up in cluster
                                 )
                               }) %>% 
          lst(iter_data = ., cluster_size) %>%
          pmap(~ mutate(..1, iter_num_takeup = if_else(!is.na(obs_num_takeup), obs_num_takeup, rbinom(n(), ..2, prob)))) # If not observed treatment, draw random number of takers for cluster
      ) %>% 
      # Change nesting to be at the treatment level
      map(mutate, iter_data = map(iter_data, group_by_at, vars(treatment_index, matches("^(mu_)?assigned_(treatment|dist_group)$"), obs_num_takeup)) %>%  
            map(group_nest,  .key = "iter_data")) %>% 
      map(unnest, iter_data),
    
    y_rate_of_change = pmap(lst(mu_rep, total_error_sd, fit_type), simulate_social_multiplier),
    
    est_takeup_level = map(cluster_cf_benefit_cost, organize_by_treatment, mu_assigned_treatment, assigned_treatment, assigned_dist_group) %>% 
      map2(cluster_cf_benefit_cost, 
           ~ filter(.y, assigned_dist_group_obs == assigned_dist_group) %>% 
             organize_by_treatment(mu_assigned_treatment, assigned_treatment) %>% 
             bind_rows(.x, .))
  ) 

# Combine models using stacking
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
      
      y_rate_of_change = list(
        reduce2(.$y_rate_of_change, .$stacking_weight_by_type, rate_of_change_stack_reducer, .init = tibble()) %>% 
          reduce(exprs(iter_prob_takeup, iter_social_multiplier, iter_partial_bbar), prep_multiple_est, .init = .)
      ) 
    ) 
}

dist_fit_data %<>%
  right_join(select(model_info, model, model_name), by = "model") %>% 
  mutate(model = factor(model, levels = model_info$model)) %>% 
  arrange(model)

# Select all the estimate differences that need to be calculated
ate_combo <- dist_fit_data %>% 
  filter(fct_match(model, "STRUCTURAL_QUADRATIC")) %$%
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
        mean_est = map_dbl(iter_data, ~ mean(.$iter_takeup_te)),
        takeup_te_quantiles = map(iter_data, quantilize_est, iter_takeup_te, wide = TRUE, quant_probs = c(quant_probs))) %>% 
    map(unnest, takeup_te_quantiles),
    
    # Calculate differences in the rate of change of E[Y] wrt to benefit-cost
    diff_y_rate_of_change = y_rate_of_change %>% 
      map_if(~ !is_null(.x), select, assigned_treatment, v, iter_data) %>% 
      map_if(~ !is_null(.x),
           function(level_data, ate_combo) {
             ate_combo %>% 
               select(str_c("assigned_treatment", c("_left", "_right"))) %>% 
               distinct_all(.keep_all = TRUE) %>% 
               inner_join(level_data, by = c("assigned_treatment_left" = "assigned_treatment")) %>% 
               inner_join(level_data, by = c("assigned_treatment_right" = "assigned_treatment", "v"), suffix = c("_left", "_right"))
          },
         ate_combo = ate_combo) %>% 
      map_if(~ !is_null(.x), mutate, iter_data = map2(iter_data_left, iter_data_right, inner_join, by = "iter_id", suffix = c("_left", "_right")) %>% 
            map(mutate, 
                iter_diff_prob_takeup = iter_prob_takeup_left - iter_prob_takeup_right,
                iter_diff_social_multiplier = iter_social_multiplier_left - iter_social_multiplier_right,
                iter_diff_partial_bbar = iter_partial_bbar_left - iter_partial_bbar_right)) %>% 
      map_if(~ !is_null(.x), select, -iter_data_left, -iter_data_right) %>% 
      map_if(~ !is_null(.x), ~ reduce(exprs(iter_diff_prob_takeup, iter_diff_social_multiplier, iter_diff_partial_bbar), prep_multiple_est, .init = .)), 
   
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
          takeup_te_quantiles = map(iter_data, quantilize_est, iter_takeup_dist_te, wide = TRUE, quant_probs = quant_probs)) %>% 
      map(unnest, takeup_te_quantiles),
    
      est_takeup_level = map(est_takeup_level, mutate, iter_data = map(iter_data, as_tibble))  
  )

# Posterior samples for assigned distance model
group_dist_param <- dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), !is_null(fit), fct_match(model_type, "structural")) %>% 
  pull(fit) %>% 
  first() %>% 
  as.data.frame(pars = c("group_dist_mean", "group_dist_sd", "group_dist_mix")) %>% 
  sample_n(1500) %>% 
  mutate(iter_id = seq(n())) %>% 
  pivot_longer(names_to = c(".value", "assigned_dist_group", "mix_index"), names_pattern = "([^\\[]+)\\[(\\d),(\\d)", cols = -iter_id) %>% 
  mutate(assigned_dist_group = factor(assigned_dist_group, levels = 1:2, labels = stan_data$cluster_treatment_map[, 2, drop = TRUE] %>% levels()))

save(dist_fit_data, stan_data, group_dist_param, file = file.path("temp-data", str_interp("processed_dist_fit${fit_version}.RData")))
