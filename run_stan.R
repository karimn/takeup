#!/opt/microsoft/ropen/3.4.3/lib64/R/bin/Rscript

script_options <- docopt::docopt(sprintf(
"Usage:
  run_stan dynamic [--analysis-data-only |[--gumbel --num-chains=<num-chains> --num-iterations=<iterations> --adapt-delta=<adapt-delta> --max-treedepth=<max-treedepth> --include-latent-var-data]] [--separate-private-value --include-name-matched --no-private-value-interact --output-name=<output-name> --output-dir=<output-dir>]
  run_stan static [--analysis-data-only |[--num-chains=<num-chains> --num-iterations=<iterations> --adapt-delta=<adapt-delta> --max-treedepth=<max-treedepth>]] [--separate-private-value --sms-control-only --include-name-matched --no-private-value-interact --model-levels=<model-levels> --output-name=<output-name> --output-dir=<output-dir>]

 Options:
  --num-chains=<num-chains>, -c <num-chains>  Number of Stan chains [default: 1]
  --num-iterations=<iterations>, -i <iterations>  Number of sampling iterations [default: 300]
  --adapt-delta=<adapt-delta>, -d <adapt-delta>  Stan control adapt_delta [default: 0.8]
  --max-treedepth=<max-treedepth>, -t <max-treedepth>  Stan control max_treedepth [default: 15]
  --output-name=<output-name>, -o <output-name>  Name to use in stanfit .csv files and analysis data .RData file [default: param]
  --output-dir=<output-dir>, -p <output-dir>  Directory analysis data and stanfit output will be stored [default: %s]
  --analysis-data-only  Don't run sampling, just produce analysis data
  --sms-control-only  Exclude SMS treatment data
  --separate-private-value  Use separate private values for calendars and bracelets
  --include-name-matched  Include unmonitored sample (name matched against census)
  --no-private-value-interact  Allow private value to interact with distance
  --model-levels=<model-levels>  How deep is the multilevel model [default: 3]
  --gumbel, -g  Gumbel link
  --include-latent-var-data, -l  Save latent variable while sampling", getwd())) 

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

#   Analysis Data -----------------------------------------------------------

subgroups <- "phone_owner"

stan_analysis_data <- analysis.data %>%
  mutate(private_value = fct_collapse(assigned.treatment, 
                                      control = if (script_options$`separate-private-value`) c("control", "ink", "bracelet") else c("control", "ink"), 
                                      "calendar" = if (script_options$`separate-private-value`) "calendar" else c("calendar", "bracelet")), 
         social_value = fct_collapse(assigned.treatment, control = c("control", "calendar")),
         sms.treatment.2 = fct_recode(sms.treatment.2, control = "sms.control")) 

all_ate <- get_dyn_ate() 
  # bind_rows(Busia = ., Kakamega = ., Siaya = ., .id = "stratum")

if (script_options$`no-private-value-interact`) {
  treatment_formula <- ~ (private_value + social_value * dist.pot.group) * phone_owner
} else {
  treatment_formula <- ~ (private_value + social_value) * dist.pot.group * phone_owner
}

if (script_options$`separate-private-value`) {
  all_ate %<>% 
    mutate(
      private_value_left = if_else(social_value_left == "bracelet", "control", private_value_left),
      private_value_right = if_else(social_value_right == "bracelet", "control", private_value_right),
      incentive_shift_left = if_else(signal_observed_left == "bracelet", "control", incentive_shift_left),
      incentive_shift_right = if_else(signal_observed_right == "bracelet", "control", incentive_shift_right))
}

if (script_options$dynamic || script_options$`sms-control-only`) {
  stan_analysis_data %<>% 
    filter(!name_matched | script_options$`include-name-matched`, 
           sms.treatment.2 == "control") 
  
  if (script_options$`include-name-matched`) {
    subgroups %<>% c("name_matched")
    
    static_treatment_map <- stan_analysis_data %>% 
      data_grid(private_value, social_value, dist.pot.group, phone_owner, name_matched) %>% 
      filter(private_value == "control" | social_value != "ink") %>%
      prepare_treatment_map()
    
    treatment_formula %<>% 
      update.formula(~ . * name_matched)
  } else {
    static_treatment_map <- stan_analysis_data %>% 
      data_grid(private_value, social_value, dist.pot.group, phone_owner) %>% 
      filter(private_value == "control" | social_value != "ink") %>%
      prepare_treatment_map()
    
     all_ate %<>% 
      filter(!name_matched) %>% 
      select(-name_matched) 
  }
  
  all_ate %<>%
    filter(sms.treatment.2_left == "control", sms.treatment.2_right == "control") %>% 
    select(-starts_with("reminder_info_stock"), -starts_with("sms.treatment")) %>% 
    distinct()
  
  if (!script_options$dynamic) {
    all_ate %<>% 
      select(-starts_with("reminder_info_stock"), -starts_with("signal_observed"), -starts_with("incentive_shift"), -starts_with("dyn_dist_pot")) %>% 
      distinct()
  }
} else {
  stan_analysis_data %<>% 
    filter(!name_matched | (script_options$`include-name-matched` & sms.treatment.2 == "control")) 
  
  if (script_options$`include-name-matched`) {
    subgroups %<>% c("name_matched")
    
    static_treatment_map <- stan_analysis_data %>% 
      data_grid(private_value, social_value, sms.treatment.2, dist.pot.group, phone_owner, name_matched) %>%
      filter(sms.treatment.2 == "control" | phone_owner,
             !name_matched | sms.treatment.2 == "control",
             sms.treatment.2 != "reminder.only" | (private_value == "control" & social_value == "control"),
             private_value == "control" | social_value != "ink") %>%
      prepare_treatment_map()
    
    treatment_formula %<>% 
      update.formula(~ . * name_matched)
  } else {
    static_treatment_map <- stan_analysis_data %>% 
      data_grid(private_value, social_value, sms.treatment.2, dist.pot.group, phone_owner) %>% #, name_matched) %>%
      filter(sms.treatment.2 == "control" | phone_owner,
             sms.treatment.2 != "reminder.only" | (private_value == "control" & social_value == "control"),
             private_value == "control" | social_value != "ink") %>%
      prepare_treatment_map()
    
     all_ate %<>% 
      filter(!name_matched) %>% 
      select(-name_matched) 
  }
  
  treatment_formula %<>% 
    update.formula(~ . + (social_value * dist.pot.group) : sms.treatment.2 + sms.treatment.2)
  
  all_ate %<>% 
    select(-starts_with("reminder_info_stock"), -starts_with("signal_observed"), -starts_with("incentive_shift"), -starts_with("dyn_dist_pot")) %>% 
    distinct()
  
}

param_stan_data <- prepare_bayesian_analysis_data(
  stan_analysis_data,
  wtp.data, 
   
  prepared_treatment_maps = TRUE, 
  treatment_map = static_treatment_map,
    
  treatment_formula = treatment_formula,
  subgroup_col = subgroups,
  drop_intercept_from_dm = FALSE, 
  nonparam_dynamics = script_options$dynamic,
  param_poly_order = 2,
  
  all_ate = all_ate,
 
  scale_sigma = 1,
  cluster_scale_sigma = 0.25,
  hyper_coef_sigma = 1,
  hyper_intercept_sigma = 5,
  
  lkj_df = 2,
  
  use_logit = !(script_options$dynamic && script_options$gumbel),
  
  estimate_ate = 1,
  
  model_levels = as.integer(script_options$`model-levels`)
)

save(param_stan_data, file = file.path(script_options$`output-dir`, "stan_analysis_data", str_interp("model_${fit_version}.RData")))

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
           # include = script_options$`include-latent-var-data`, pars = if (script_options$`include-latent-var-data`) NA else c("cluster_latent_var_map"),           
           control = lst(max_treedepth = as.integer(script_options$`max-treedepth`), 
                         adapt_delta = as.numeric(script_options$`adapt-delta`)), 
           init = if (script_options$dynamic && script_options$gumbel) gen_initializer(.) else "random",
           sample_file = file.path(script_options$`output-dir`, "stanfit", str_interp("model_${fit_version}.csv")))

