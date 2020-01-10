library(magrittr)
library(tidyverse)
library(rstan)
library(pbmcapply)
library(nleqslv)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

load(file.path("stan_analysis_data", "dist_fit21.RData"))

dist_fit_data <- enframe(dist_fit, name = "model", value = "fit")

rm(dist_fit)

load(file.path("stan_analysis_data", "dist_kfold_16.RData"))

dist_fit_data %<>%
  left_join(enframe(dist_kfold, name = "model", value = "kfold") %>% 
              mutate(stacking_weight = map(kfold, pluck, "pointwise") %>% 
                       do.call(rbind, .) %>% 
                       t() %>% 
                       loo::stacking_weights()),
            by = "model") %>% 
  mutate(kfold = set_names(kfold, model)) %>% 
  left_join(kfold_compare(x = discard(.$kfold, is_null)), by = "model") 

rm(dist_kfold)

# load(file.path("stan_analysis_data", "dist_fit_prior18.RData"))
# 
# dist_fit_data %<>%
#   bind_rows("fit" = .,
#             "prior-predict" = enframe(dist_fit, name = "model", value = "fit"))


model_names <- c(
  "REDUCED_FORM_NO_RESTRICT" = "Reduced Form Discrete Cost",
  "REDUCED_FORM_SEMIPARAM" = "Reduced Form Semiparametric Cost",
  "STRUCTURAL_LINEAR" = "Structural Linear Cost",
  "STRUCTURAL_LINEAR_SALIENCE" = "Structural Linear Cost With Salience"
)

thin_model <- c(
  "REDUCED_FORM_NO_RESTRICT" = 1L,
  "REDUCED_FORM_SEMIPARAM" = 1L,
  "STRUCTURAL_LINEAR" = 1L
)

dist_fit_data %<>% 
  left_join(enframe(model_names, name = "model", value = "model_name"), by = "model") %>% 
  left_join(enframe(thin_model, name = "model", value = "thin"), by = "model") %>% 
  mutate(
    model = factor(model, levels = names(model_names)),
    thin = coalesce(thin, 1L)
  )

observed_takeup <- monitored_nosms_data %>% 
  group_by(cluster_id, assigned_treatment = assigned.treatment, assigned_dist = cluster.dist.to.pot) %>% 
  summarize(prop_takeup = mean(dewormed), 
            se = sqrt(prop_takeup * (1 - prop_takeup) / n()),
            prop_takeup_ub = prop_takeup + se,
            prop_takeup_lb = prop_takeup - se) %>% 
  ungroup()

quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
postprocess_cores <- 12


dist_fit_data %<>% 
  mutate(
    total_error_sd = map(fit, extract_obs_fit_level, stan_data = stan_data, par = "total_error_sd", iter_level = "none"),
    # cluster_dist_cost = map(fit, extract_obs_fit_level, par = "cluster_dist_cost", stan_data = stan_data, iter_level = "cluster", quant_probs = quant_probs),
    # net_benefit = map(fit, extract_obs_fit_level, par = "structural_cluster_benefit_cost", stan_data = stan_data, iter_level = "cluster", quant_probs = quant_probs),
    # v_cutoff = map(fit, extract_obs_fit_level, par = "structural_cluster_obs_v", stan_data = stan_data, iter_level = "cluster", quant_probs = quant_probs),
    # takeup_prob = map(fit, extract_obs_fit_level, par = "structural_cluster_takeup_prob", stan_data = stan_data, iter_level = "cluster", mix = TRUE, quant_probs = quant_probs),
    # structural_beta = map(fit, extract_obs_fit_level, par = "beta", stan_data = stan_data, iter_level = "treatment", mix = FALSE, quant_probs = quant_probs),
    # dist_beta_v = map(fit, extract_obs_fit_level, par = "dist_beta_v", stan_data = stan_data, iter_level = "treatment", mix = FALSE, quant_probs = quant_probs),
    mu_rep = map(fit, extract_obs_fit_level, par = "mu_rep", stan_data = stan_data, iter_level = "treatment", mix = FALSE, quant_probs = quant_probs) %>% 
        map_if(~ !is_null(.), mutate, compare_to = "control") %>% 
        map_if(~ !is_null(.), ~ left_join(., select(., assigned_treatment, iter_data), by = c("compare_to" = "assigned_treatment"), suffix = c("", "_compare"))) %>% 
        map_if(~ !is_null(.),
            mutate, 
            iter_data = map2(iter_data, iter_data_compare, left_join, by = "iter_id", suffix = c("", "_compare")) %>% 
              map(mutate, 
                  iter_est_vs_compare = iter_est - iter_est_compare,
                  iter_est_vs_compare_ratio = iter_est / iter_est_compare,
                  per = if_else(rank(iter_est) %in% c((n() + 1) %/% 2, ((n() + 1) %/% 2) + ((n() + 1) %% 2)), 0.5, NA_real_)
              ) %>% 
              map2(quantiles_est, ~ mutate(.x, 
                                           per_group = unlist(.y) %>%
                                             set_names(str_extract(names(.), "0\\.\\d+")) %>% { 
                                               cut(iter_est, breaks = c(0, ., Inf), labels = c("0.0", names(.))) } %>% 
                                             as.numeric())),
            mean_est_vs_compare = map_dbl(iter_data, ~ mean(.$iter_est_vs_compare)),
            mean_est_vs_compare_ratio = map_dbl(iter_data, ~ mean(.$iter_est_vs_compare_ratio)),
            quantiles_est_vs_compare = map(iter_data, quantilize_est, iter_est_vs_compare, quant_probs = quant_probs), 
            quantiles_est_vs_compare_ratio = map(iter_data, quantilize_est, iter_est_vs_compare_ratio, quant_probs = quant_probs)
        ),
    reduced_mu_rep = map_if(mu_rep, ~ !is_null(.), mutate, iter_data = map(iter_data, arrange, iter_est) %>% map(filter, (row_number() %% 100) == 0)),
    
    # structural_param = pmap(lst(cluster_dist_cost, net_benefit, v_cutoff, takeup_prob, mu_rep),
    #                  function(cluster_dist_cost, net_benefit, v_cutoff, takeup_prob, mu_rep) {
    #                    inner_join(net_benefit, v_cutoff, by = c("obs_index", "assigned_treatment", "assigned_dist", "assigned_dist_standard"), suffix = c("", "_v_cutoff")) %>% 
    #                      inner_join(cluster_dist_cost, 
    #                                 by = c("obs_index", "assigned_treatment", "assigned_dist", "assigned_dist_standard"), 
    #                                 suffix = c("", "_cluster_dist_cost")) %>% { 
    #                        if (!is_null(mu_rep)) {
    #                          inner_join(., takeup_prob, by = c("obs_index", "assigned_treatment", "assigned_dist", "assigned_dist_standard"), suffix = c("", "_takeup_prob")) %>% 
    #                            left_join(., mu_rep, by = c("assigned_treatment"), suffix = c("_net_benefit", "_mu_rep")) 
    #                        } else {
    #                          inner_join(., 
    #                                     takeup_prob, 
    #                                     by = c("obs_index", "assigned_treatment", "assigned_dist", "assigned_dist_standard"), 
    #                                     suffix = c("_net_benefit", "_takeup_prob")) 
    #                        }
    #                      } %>% 
    #                       mutate(
    #                         iter_data = map2(iter_data_net_benefit, iter_data_v_cutoff, inner_join, by = "iter_id", suffix = c("", "_v_cutoff")) %>% {
    #                           if (!is_null(mu_rep)) {
    #                             map2(., iter_data_takeup_prob, inner_join, by = "iter_id", suffix = c("", "_takeup_prob")) %>% 
    #                               map2(iter_data_mu_rep, left_join, by = "iter_id", suffix = c("_net_benefit", "_mu_rep")) 
    #                           } else {
    #                             map2(., iter_data_takeup_prob, inner_join, by = "iter_id", suffix = c("_net_benefit", "_takeup_prob"))  
    #                           }
    #                         } %>%  
    #                           map(mutate, 
    #                               iter_est_rep = - iter_est_v_cutoff - iter_est_net_benefit)
    #                       ) %>% 
    #                       left_join(select(observed_takeup, -c(assigned_treatment, assigned_dist)), by = c("obs_index" = "cluster_id")
    #                       )
    #                  }),
    
    cluster_cf_benefit_cost = map2(fit, thin, ~ extract_obs_cf(.x, par = "cluster_cf_benefit_cost", stan_data = stan_data, iter_level = "cluster", quant_probs = quant_probs, thin = .y)) %>% 
      # map(select, -starts_with("ess"), -one_of("rhat")) %>% 
      map(nest, iter_data = c(treatment_index, obs_num_takeup, matches("^assigned_(treatment|dist_group)$"), iter_data)) %>% 
      map(mutate, iter_data = map(iter_data, unnest, iter_data)) %>% 
      lst(benefit_cost = ., 
          mu_rep, 
          total_error_sd = map(total_error_sd, select, -rhat, -starts_with("ess")) %>% 
            map(unnest, iter_data)) %>%
      pmap(function(benefit_cost, mu_rep, total_error_sd) {
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
        } else {
          mutate(benefit_cost, 
                 iter_data = pbmclapply(mc.cores = postprocess_cores,
                                        iter_data, 
                                        inner_join, 
                                        total_error_sd, 
                                        by = c("iter_id"), 
                                        suffix = c("", "_error_sd")) %>% 
                   map(mutate, mu_assigned_treatment = NA_character_, iter_est_mu = 0))
        } 
      }) %>% 
      map(
        mutate,
        iter_data = pbmclapply(mc.cores = postprocess_cores,
                               iter_data,
                               function(curr_iter_data) {
                                 mutate(
                                   curr_iter_data,
                                   no_rep_prob = pnorm(iter_est, sd = iter_est_error_sd),
                                   v_cutoff = map2_dbl(iter_est, iter_est_mu, ~ nleqslv(x = - ..1, fn = generate_v_cutoff_fixedpoint(..1, ..2)) %>% pluck("x")),
                                   delta_rep = rep_normal(v_cutoff), 
                                   prob = pnorm(- v_cutoff, sd = iter_est_error_sd)
                                 )
                               }) %>% 
          lst(iter_data = ., cluster_size) %>%
          pmap(~ mutate(..1, iter_num_takeup = if_else(!is.na(obs_num_takeup), obs_num_takeup, rbinom(n(), ..2, prob))))
      ) %>% 
      map(mutate, iter_data = map(iter_data, group_by_at, vars(treatment_index, matches("^(mu_)?assigned_(treatment|dist_group)$"), obs_num_takeup)) %>% 
            map(group_nest,  .key = "iter_data")) %>% 
      map(unnest, iter_data),
    
    est_takeup_level = map(cluster_cf_benefit_cost,
      function(.data) {
        .data %>%
          select(
            cluster_id, 
            assigned_treatment, assigned_dist_group, mu_assigned_treatment, 
            assigned_dist_obs, assigned_treatment_obs, assigned_dist_group_obs, cluster_size, 
            iter_data
          ) %>% 
          mutate(iter_data = map(iter_data, select, iter_id, prob, iter_num_takeup)) %>%
          unnest(iter_data) %>% 
          group_by(iter_id, mu_assigned_treatment, assigned_treatment, assigned_dist_group) %>% 
          summarize(iter_prop_takeup = sum(iter_num_takeup) / sum(cluster_size)) %>% 
          ungroup() %>% 
          nest(iter_data = c(iter_id, iter_prop_takeup)) %>%
          mutate(
            mean_est = map_dbl(iter_data, ~ mean(.$iter_prop_takeup)),
            takeup_quantiles = map(iter_data, quantilize_est, iter_prop_takeup, wide = TRUE, quant_probs = c(quant_probs))
          ) %>% 
          unnest(takeup_quantiles)
      }),
    
    # stacked_est_takeup_level = map2(stacking_weight, est_takeup_level, {
    #   if (!is.na(.x)) {
    #     .y %>% 
    #       select(-mean_est, -starts_with("per_")) %>% 
    #       mutate(iter_data = map(iter_data, mutate, iter_prop_takeup = iter_prop_takeup * .x))
    #   } else {
    #     return(NULL)
    #   }
    # }),
  ) 

dist_fit_data %<>% 
  add_row(
    model = "STACKED", 
    model_name = "Stacked Model", 
    est_takeup_level = list(
      reduce2(.$est_takeup_level, .$stacking_weight, 
              function(accum, est_takeup_level, stacking_weight) {
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
              },
              .init = tibble()) %>% 
        mutate(
          mean_est = map_dbl(iter_data, ~ mean(.$iter_prop_takeup)),
          takeup_quantiles = map(iter_data, quantilize_est, iter_prop_takeup, wide = TRUE, quant_probs = c(quant_probs))
        ) %>% 
        unnest(takeup_quantiles)
    ))

ate_combo <- dist_fit_data %>% 
  filter(fct_match(model, "STRUCTURAL_LINEAR")) %$%
  select(est_takeup_level[[1]], mu_assigned_treatment:assigned_dist_group) %>% {
    bind_cols(rename_all(., str_c, "_left"), rename_all(., str_c, "_right"))
  } %>%
  expand(crossing(!!!syms(names(.)))) %>%
  filter(mu_assigned_treatment_left != mu_assigned_treatment_right |
           assigned_treatment_left != assigned_treatment_right |
           (fct_match(assigned_dist_group_left, "far") & fct_match(assigned_dist_group_right, "close")))

dist_fit_data %<>% 
  mutate(
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
      ~ any(str_detect(names(.), "mu_assigned_treatment")),
      ~ filter(., 
        mu_assigned_treatment_left != mu_assigned_treatment_right |
        assigned_treatment_left != assigned_treatment_right |
        assigned_dist_group_left != assigned_dist_group_right),
      .else = ~ filter(., 
        assigned_treatment_left != assigned_treatment_right |
        assigned_dist_group_left != assigned_dist_group_right)) %>% 
    map(mutate, iter_data = map2(iter_data_left, iter_data_right, inner_join, by = "iter_id", suffix = c("_left", "_right")) %>% 
          map(mutate, iter_takeup_te = iter_prop_takeup_left - iter_prop_takeup_right)) %>% 
    map(select, -iter_data_left, -iter_data_right) %>% 
    map(mutate, 
        mean_est = map_dbl(iter_data, ~ mean(.$iter_takeup_te)),
        takeup_te_quantiles = map(iter_data, quantilize_est, iter_takeup_te, wide = TRUE, quant_probs = c(quant_probs))) %>% 
    map(unnest, takeup_te_quantiles),
    
    est_takeup_dist_te = est_takeup_te %>% 
      map_if(~ any(str_detect(names(.), "mu_assigned_treatment")), filter, mu_assigned_treatment_left == assigned_treatment_left, mu_assigned_treatment_right == assigned_treatment_right) %>% 
      map_if(~ any(str_detect(names(.), "mu_assigned_treatment")), select, -starts_with("mu_assigned_treatment")) %>% 
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
      map(unnest, takeup_te_quantiles)
  )

dist_fit_data %<>% 
  mutate(est_takeup_level = map(est_takeup_level, mutate, iter_data = map(iter_data, as_tibble)))  

group_dist_param <- dist_fit[[1]] %>% 
  as.data.frame(pars = c("group_dist_mean", "group_dist_sd", "group_dist_mix")) %>% 
  sample_n(1500) %>% 
  mutate(iter_id = seq(n())) %>% 
  pivot_longer(names_to = c(".value", "assigned_dist_group", "mix_index"), names_pattern = "([^\\[]+)\\[(\\d),(\\d)", cols = -iter_id) %>% 
  mutate(assigned_dist_group = factor(assigned_dist_group, levels = 1:2, labels = stan_data$cluster_treatment_map[, 2, drop = TRUE] %>% levels()))

save(dist_fit_data, stan_data, group_dist_param, file = file.path("temp-data", "processed_dist_fit21.RData")) 
