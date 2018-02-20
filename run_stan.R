#!/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript

script_options <- docopt::docopt(
"Usage:
  run_stan dynamic [--analysis-data-only |[--gumbel --num-chains=<num-chains> --num-iterations=<iterations> --adapt-delta=<adapt-delta>]] [--output-name=<output-name>]
  run_stan static [--analysis-data-only |[--num-chains=<num-chains> --num-iterations=<iterations> --adapt-delta=<adapt-delta>]] [--output-name=<output-name>]

 Options:
  --num-chains=<num-chains>, -c <num-chains>  Number of Stan chains [default: 1]
  --num-iterations=<iterations>, -i <iterations>  Number of sampling iterations [default: 300]
  --adapt-delta=<adapt-delta>, -d <adapt-delta>  Stan control adapt_delta [default: 0.8]
  --output-name=<output-name>, -o <output-name>  Name to use in stanfit .csv files and analysis data .RData file [default: param]
  --analysis-data-only  Don't run sampling, just produce analysis data
  --gumbel, -g  Gumbel link"
)

library(magrittr)
library(plyr)
library(forcats)
library(broom)
library(tidyverse)
library(modelr)
library(rstan)

source("analysis_util.R")

load(file.path("data", "analysis.RData"))

fit_version <- script_options$`output-name`

# Analysis Data -----------------------------------------------------------

if (script_options$dynamic) {
  stan_analysis_data <- analysis.data %>%
    mutate(private_value = fct_collapse(assigned.treatment, 
                                        control = c("control", "ink"), 
                                        "calendar" = c("calendar", "bracelet")), 
           social_value = fct_collapse(assigned.treatment, control = c("control", "calendar")),
           sms.treatment.2 = fct_recode(sms.treatment.2, control = "sms.control")) %>% 
    filter(!name_matched, sms.treatment.2 == "control") #, sms.treatment.2 == "control") #, !hh.baseline.sample)
  
  static_treatment_map <- param_dyn_analysis_data %>% 
    #data_grid(private_value, social_value, sms.treatment.2, dist.pot.group, phone_owner) %>% #, name_matched) %>% 
    data_grid(private_value, social_value, dist.pot.group, phone_owner) %>% #, name_matched) %>% 
    #filter(sms.treatment.2 == "control" | phone_owner,
           # !name_matched | sms.treatment.2 == "control",
           # sms.treatment.2 != "reminder.only" | (private_value == "control" & social_value == "control"),
    filter(private_value == "control" | social_value != "ink") %>%
    prepare_treatment_map()
  
  all_ate <- get_dyn_ate() %>% 
    filter(!name_matched) %>% 
    select(-name_matched) %>%
    filter(sms.treatment.2_left == "control",
           sms.treatment.2_right == "control") %>% 
    select(-starts_with("sms.treatment"), -starts_with("reminder_info_stock")) %>% 
    distinct()
} else {
  stan_analysis_data <- analysis.data %>%
    mutate(private_value = fct_collapse(assigned.treatment, 
                                        control = c("control", "ink"), 
                                        "calendar" = c("calendar", "bracelet")), 
           social_value = fct_collapse(assigned.treatment, control = c("control", "calendar")),
           sms.treatment.2 = fct_recode(sms.treatment.2, control = "sms.control")) %>% 
    filter(!name_matched, sms.treatment.2 == "control") 
  
  static_treatment_map <- param_dyn_analysis_data %>% 
    #data_grid(private_value, social_value, sms.treatment.2, dist.pot.group, phone_owner) %>% #, name_matched) %>% 
    data_grid(private_value, social_value, dist.pot.group, phone_owner) %>% #, name_matched) %>% 
    #filter(sms.treatment.2 == "control" | phone_owner,
           # !name_matched | sms.treatment.2 == "control",
           # sms.treatment.2 != "reminder.only" | (private_value == "control" & social_value == "control"),
    filter(private_value == "control" | social_value != "ink") %>%
    prepare_treatment_map()
  
  all_ate <- get_dyn_ate() %>% 
    filter(!name_matched) %>% 
    select(-name_matched) %>%
    filter(sms.treatment.2_left == "control",
           sms.treatment.2_right == "control") %>% 
    select(-starts_with("sms.treatment"), -starts_with("reminder_info_stock")) %>% 
    select(-starts_with("signal_observed"), -starts_with("incentive_shift"), -starts_with("dyn_dist_pot")) %>% 
    distinct()
}

param_stan_data <- prepare_bayesian_analysis_data(
  stan_analysis_data,
  wtp.data, 
   
  prepared_treatment_maps = TRUE, 
  treatment_map = static_treatment_map,
    
  treatment_formula = ~ (private_value + social_value * dist.pot.group) * phone_owner,
  #treatment_formula = ~ (private_value + social_value * dist.pot.group) * phone_owner + (social_value * dist.pot.group) : sms.treatment.2 + sms.treatment.2,
  # treatment_formula = ~ (private_value + social_value * dist.pot.group) * phone_owner * name_matched + (social_value * dist.pot.group) : sms.treatment.2 + sms.treatment.2,
  subgroup_col = "phone_owner",
  drop_intercept_from_dm = FALSE, 
  nonparam_dynamics = script_options$dynamic,
  param_poly_order = 2,
  
  all_ate = all_ate,
  
  scale_sigma = 1,
  cluster_scale_sigma = 0.01,
  cluster_intercept_scale_sigma = 1,
  hyper_coef_sigma = 1,
  hyper_intercept_sigma = 5,
  
  lkj_df = 2,
  
  use_logit = !(script_options$dynamic && script_options$gumbel),
  
  estimate_ate = 1
)

save(param_stan_data, file = file.path("stan_analysis_data", str_interp("model_${fit_version}.RData")))

if (script_options$`analysis-data-only`) quit()

# Initializer Factory -------------------------------------------------------------

gen_initializer <- function(stan_data_list, dynamic = TRUE) {
  if (dynamic) {
    function() {
      lst(
        strata_beta_day1_corr_mat_non_phone = with(stan_data_list, rethinking::rlkjcorr(1, subgroup_treatment_col_sizes[1], lkj_df)),
        strata_beta_day1_corr_mat_phone = with(stan_data_list, rethinking::rlkjcorr(1, subgroup_treatment_col_sizes[2], lkj_df)),
        strata_beta_day1_L_corr_mat_non_phone = t(chol(strata_beta_day1_corr_mat_non_phone)),
        strata_beta_day1_L_corr_mat_phone = t(chol(strata_beta_day1_corr_mat_phone)),
        
        # hyper_beta_day1 = rep(0, stan_data_list$num_all_treatment_coef),
        # hyper_beta_day1 = rnorm(stan_data_list$num_all_treatment_coef),
        # hyper_baseline_dyn_effect = rep(0, stan_data_list$num_deworming_days - 1),
        # hyper_baseline_dyn_effect = rnorm(stan_data_list$num_deworming_days - 1),
        # hyper_treat_beta_dyn_effect = rep(0, stan_data_list$num_param_dyn_coef),
        strata_baseline_dyn_effect_raw = with(stan_data_list, matrix(rnorm((num_deworming_days - 1) * num_strata, sd = 0.5), nrow = num_deworming_days - 1, ncol = num_strata)),
        QR_strata_beta_day1 = with(stan_data_list, matrix(rnorm(num_all_treatment_coef * num_strata, sd = 0.005), nrow = num_all_treatment_coef, ncol = num_strata)),
        QR_strata_beta_dyn_effect = with(stan_data_list, matrix(rnorm(num_param_dyn_coef * num_strata, sd = 0.005), nrow = num_param_dyn_coef, ncol = num_strata)),
        # cluster_effect = rep(0, stan_data_list$num_clusters)
      )
    } 
  } else {
    function() {}
  }
}

# Run ---------------------------------------------------------------------

num_chains <- as.integer(script_options$`num-chains`)

if (num_chains > parallel::detectCores()) stop("Not enough cores.")

options(mc.cores = num_chains)
rstan_options(auto_write = TRUE)

if (script_options$dynamic) {
  model_param <- stan_model(file = file.path("stan_models", "takeup_model_3_param.stan"), model_name = "model_3_param")
} else {
  model_param <- stan_model(file = file.path("stan_models", "takeup_model_4_static_param.stan"), model_name = "model_4_static_param")
}

cat(str_interp("Output name: ${fit_version}\n"))

model_fit <- param_stan_data %>% 
  sampling(model_param, data = ., 
           chains = num_chains,
           iter = as.integer(script_options$`num-iterations`),
           control = lst(max_treedepth = 15, adapt_delta = as.numeric(script_options$`adapt-delta`)), 
           init = if (script_options$dynamic && script_options$gumbel) gen_initializer(.) else "random",
           sample_file = file.path("stanfit", str_interp("model_3_${fit_version}.csv")))

