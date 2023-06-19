#!/usr/bin/Rscript
script_options <- docopt::docopt(
  stringr::str_glue(
"Usage:
  quick_ate_postprocess.R <fit-version> [options] [<chain>...]
  
Options:
  --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
  --output-path=<path>  Path to find results [default: temp-data]
  --model=<model>  Which model to postprocess
  --prior  Postprocess the prior predictive
  
  "), 
  args = if (interactive()) "
  95
  --output-path=temp-data
  --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_STRATA_SOB
  3 4
  " else commandArgs(trailingOnly = TRUE)
)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(tidybayes)



## Load analysis data
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
analysis_data <- monitored_nosms_data %>% 
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group)

## Load Stan output
load_param_draws = function(fit_version, model, chain, prior_predictive = FALSE, ...) {
  if (prior_predictive == TRUE) {
    fit_str = str_glue(
      "data/stan_analysis_data/dist_prior{fit_version}_{model}-{chain}.csv"
    )
  } else {
    fit_str = str_glue(
      "data/stan_analysis_data/dist_fit{fit_version}_{model}-{chain}.csv"
    )
  }

  fit_obj = as_cmdstan_fit(fit_str)
  draws = spread_rvars(
    fit_obj,
    ...
  ) %>%
    mutate(model = model, fit_version = fit_version, fit_type = if_else(prior_predictive, "prior-predict", "fit"))
  return(draws)
}

# N.B. treat_idx (the second idx, is the mu (signalling) idx)
mu_idx_mapper = tibble(
  treat_idx = 1:4,
  mu_treatment = c("control", "ink", "calendar", "bracelet")
) %>%
  mutate(mu_treatment = factor(mu_treatment, levels = c("bracelet", "calendar", "ink", "control")))
dist_idx_mapper = tibble(
  dist_treat_idx = 1:8,
  dist_treatment = rep(c("control", "ink", "calendar", "bracelet"), 2),
  dist_group = rep(c("close", "far"), each = 4)
) %>%
  mutate(dist_treatment = factor(dist_treatment, levels = c("bracelet", "calendar", "ink", "control")))
cluster_mapper = analysis_data %>%
  select(
    cluster_id,
    assigned_treatment,
    assigned_dist_group
  ) %>% unique()


cluster_error_draws_raw = load_param_draws(
  fit_version = script_options$fit_version,
  model = script_options$model,
  chain = script_options$chain,
  prior_predictive = script_options$prior,
  cluster_cf_cutoff[dist_treat_idx, treat_idx, cluster_idx],
  total_error_sd[treat_idx]
)

rvar_pnorm = rfun(pnorm)
cluster_error_draws = cluster_error_draws_raw %>%
  mutate(
    pr_takeup = rvar_pnorm(-cluster_cf_cutoff, sd = total_error_sd)
  ) %>%
  left_join(
    dist_idx_mapper,
    by = c("dist_treat_idx")
  ) %>%
  left_join(
    mu_idx_mapper,
    by = "treat_idx"
 ) %>%
  left_join(
    cluster_mapper,
    by = c("cluster_idx" = "cluster_id")
  )

fit_type_str = if_else(script_options$prior, "prior", "fit")
if (length(script_options$chain) > 1) {
  chain_str = str_glue("{min(script_options$chain)}-{max(script_options$chain)}")
} else {
  chain_str = script_options$chain
}

cluster_error_draws %>%
  saveRDS(
    file.path(
      script_options$output_path,
      str_glue(
        "rvar_processed_dist_{fit_type_str}{script_options$fit_version}_cluster_{script_options$model}_{chain_str}.rds"
      )
    ) 
  )


incentive_tes =  bind_rows(
  cluster_error_draws %>%
    filter(dist_group == assigned_dist_group) %>%
    filter(mu_treatment == dist_treatment) %>%
    group_by(fit_type, model, fit_version, dist_treatment, dist_group) %>%
    summarise(
      pr_takeup = rvar_mean(pr_takeup)
    )   %>%
    group_by(dist_group) %>%
    mutate(
      pr_takeup = if_else(dist_treatment == "control", pr_takeup, pr_takeup - pr_takeup[dist_treatment == "control"])
    ),
  cluster_error_draws %>%
    filter(dist_group == assigned_dist_group) %>%
    filter(mu_treatment == dist_treatment) %>%
    group_by(fit_type, model, fit_version, dist_treatment) %>%
    summarise(
      pr_takeup = rvar_mean(pr_takeup)
    ) %>%
    mutate(
      pr_takeup = if_else(dist_treatment == "control", pr_takeup, pr_takeup - pr_takeup[dist_treatment == "control"])
    )  %>%
    mutate(dist_group = "combined")
) %>%
  mutate(
    dist_group = factor(dist_group, levels = c("far", "close", "combined"))
  ) 


create_tes = function(.data, group_var) {
  tes =  bind_rows(
    .data %>%
      group_by(fit_type, model, fit_version, {{ group_var }}, dist_group) %>%
      summarise(
        pr_takeup = rvar_mean(pr_takeup),
        .groups = "drop"
      )   %>%
      group_by(dist_group) %>%
      mutate(
        pr_takeup = if_else({{ group_var }} == "control", pr_takeup, pr_takeup - pr_takeup[{{ group_var }} == "control"])
      ),
    .data %>%
      group_by(fit_type, model, fit_version, {{ group_var }}) %>%
      summarise(
        pr_takeup = rvar_mean(pr_takeup),
        .groups = "drop"
      ) %>%
      mutate(
        pr_takeup = if_else({{ group_var }} == "control", pr_takeup, pr_takeup - pr_takeup[{{ group_var }} == "control"])
      )  %>%
      mutate(dist_group = "combined")
  ) %>%
    mutate(
      dist_group = factor(dist_group, levels = c("far", "close", "combined"))
    ) 
  return(tes)
}
  
incentive_tes = cluster_error_draws %>%
      filter(dist_group == assigned_dist_group) %>%
      filter(mu_treatment == dist_treatment) %>%
      create_tes(group_var = dist_treatment)

signal_tes = cluster_error_draws %>%
  filter(dist_treatment == "control") %>%
  create_tes(group_var = mu_treatment)

private_tes = cluster_error_draws %>%
  filter(dist_group == assigned_dist_group) %>%
  filter(mu_treatment == "control") %>%
  create_tes(group_var = dist_treatment) %>%
  filter(dist_group == "combined")


all_tes = bind_rows(
  incentive_tes %>% mutate(estimand = "overall"),
  signal_tes %>% mutate(estimand = "signal"),
  private_tes %>% mutate(estimand = "private")
) %>%
  ungroup()

all_tes %>%
  saveRDS(
    file.path(
      script_options$output_path,
      str_glue(
        "rvar_processed_dist_{fit_type_str}{script_options$fit_version}_ates_{script_options$model}_{chain_str}.rds"
      )
    ) 
  )


## Munging stuff

# create_cis = function(.data) {
#   med_fun = function(x) {
#     mean_x = mean(x) %>% round(3)
#     conf.low = quantile(x, 0.05) %>% round(3)
#     conf.high = quantile(x, 0.95) %>% round(3)
#     return(paste0(
#       mean_x, " (", conf.low, ", ", conf.high, ")"
#     ))
#   }
#   .data %>%
#     mutate(across(where(is_rvar), med_fun))
# }

# wide_struct_tes = all_tes %>% 
#     filter(estimand == "overall") %>%
#     select(
#       dist_treatment,
#       dist_group,
#       pr_takeup
#     ) %>%
#     pivot_wider(
#       names_from = dist_group,
#       values_from = pr_takeup
#     ) %>%
#     select(dist_treatment, combined, close, far) %>%
#     arrange(dist_treatment) %>%
#     bind_rows(
#       # bracelet minus calendar row
#       all_tes %>%
#         filter(dist_treatment %in% c("bracelet", "calendar")) %>%
#         filter(estimand == "overall") %>%
#         pivot_wider(names_from = dist_treatment, values_from = pr_takeup) %>%
#         mutate(
#           bracelet_minus_calendar = bracelet - calendar
#         ) %>%
#         select(dist_group, bracelet_minus_calendar) %>%
#         pivot_wider(
#           names_from = dist_group,
#           values_from = bracelet_minus_calendar
#         ) %>%
#         mutate(dist_treatment = "bracelet_minus_calendar")
#     )

# wide_struct_tes %>%
#   mutate(dist_treatment = factor(dist_treatment, levels = c(
#     "bracelet",
#     "calendar",
#     "ink",
#     "bracelet_minus_calendar",
#     "control"
#   ))) %>%
#   mutate(far_minus_close = far - close) %>%
#   arrange(dist_treatment)   %>%
#   create_cis() %>%
#   View()
