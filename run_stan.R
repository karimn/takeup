library(magrittr)
library(plyr)
library(forcats)
library(broom)
library(tidyverse)
library(modelr)
library(rstan)
library(docopt)

script_options <- "Usage:
  run_stan [--analysis-data-only | --num-chains=<num-chains> --num-iterations=<iterations> --adapt-delta=<adapt-delta>] [--output-name=<output-name>]

 Options:
  --num-chains=<num-chains>, -c <num-chains>  Number of Stan chains [default: 1]
  --num-iterations=<iterations>, -i <iterations>  Number of sampling iterations [default: 300]
  --adapt-delta=<adapt-delta>, -d <adapt-delta>  Stan control adapt_delta [default: 0.8]
  --output-name=<output-name>, -o <output-name>  Name to use in stanfit .csv files and analysis data .RData file [default: param]
  --analysis-data-only  Don't run sampling, just produce analysis data
" %>% 
  docopt()

source("analysis_util.R")

load(file.path("data", "analysis.RData"))

dyn_fit_version <- script_options$`output-name`

# Analysis Data -----------------------------------------------------------

param_dyn_analysis_data <- analysis.data %>%
  mutate(private_value = fct_collapse(assigned.treatment, 
                                      control = c("control", "ink"), 
                                      "calendar" = c("calendar", "bracelet")), 
         social_value = fct_collapse(assigned.treatment, control = c("control", "calendar")),
         sms.treatment.2 = fct_recode(sms.treatment.2, control = "sms.control")) %>% 
  filter(!name_matched, sms.treatment.2 == "control") #, sms.treatment.2 == "control") #, !hh.baseline.sample)

dyn_static_treatment_map <- param_dyn_analysis_data %>% 
  #data_grid(private_value, social_value, sms.treatment.2, dist.pot.group, phone_owner) %>% #, name_matched) %>% 
  data_grid(private_value, social_value, dist.pot.group, phone_owner) %>% #, name_matched) %>% 
  #filter(sms.treatment.2 == "control" | phone_owner,
         # !name_matched | sms.treatment.2 == "control",
         # sms.treatment.2 != "reminder.only" | (private_value == "control" & social_value == "control"),
  filter(private_value == "control" | social_value != "ink") %>%
  prepare_treatment_map()

param_dyn_stan_data <- prepare_bayesian_analysis_data(
  param_dyn_analysis_data,
  wtp.data, 
   
  prepared_treatment_maps = TRUE, 
  treatment_map = dyn_static_treatment_map,
    
  treatment_formula = ~ (private_value + social_value * dist.pot.group) * phone_owner,
  #treatment_formula = ~ (private_value + social_value * dist.pot.group) * phone_owner + (social_value * dist.pot.group) : sms.treatment.2 + sms.treatment.2,
  # treatment_formula = ~ (private_value + social_value * dist.pot.group) * phone_owner * name_matched + (social_value * dist.pot.group) : sms.treatment.2 + sms.treatment.2,
  subgroup_col = "phone_owner",
  drop_intercept_from_dm = FALSE, 
  nonparam_dynamics = TRUE,
  param_poly_order = 2,
  
  all_ate = get_dyn_ate() %>% 
    filter(!name_matched) %>% 
    select(-name_matched) %>%
    filter(sms.treatment.2_left == "control",
           sms.treatment.2_right == "control") %>% 
    select(-starts_with("sms.treatment"), -starts_with("reminder_info_stock")) %>% 
    distinct(),
  
  scale_sigma = 1,
  cluster_scale_sigma = 0.5,
  hyper_coef_sigma = 1,
  hyper_intercept_sigma = 5,
  
  lkj_df = 2,
  
  use_logit = 1,
  
  estimate_ate = 1
)

save(param_dyn_stan_data, file = file.path("stan_analysis_data", str_interp("model_3_${dyn_fit_version}.RData")))

if (script_options$`analysis-data-only`) quit()

# Run ---------------------------------------------------------------------

num_chains <- as.integer(script_options$`num-chains`)
# num_chains <- commandArgs(trailingOnly = TRUE)[1]

# if (is.na(num_chains)) num_chains <- max(1, parallel::detectCores())
if (num_chains > parallel::detectCores()) stop("Not enough cores.")

options(mc.cores = num_chains)
rstan_options(auto_write = TRUE)

model_3_param <- stan_model(file = file.path("stan_models", "takeup_model_3_param.stan"), model_name = "model_3_param")

cat(str_interp("Detected cores ${parallel::detectCores()}\n"))

model_3_fit <- param_dyn_stan_data %>% 
  sampling(model_3_param, data = ., 
           chains = num_chains, 
           iter = as.integer(script_options$`num-iterations`),
           control = lst(max_treedepth = 15, adapt_delta = as.numeric(script_options$`adapt-delta`)), 
           sample_file = file.path("stanfit", str_interp("model_3_${dyn_fit_version}.csv")))

