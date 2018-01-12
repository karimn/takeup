library(magrittr)
library(plyr)
library(forcats)
library(broom)
library(tidyverse)
library(modelr)

library(rstan)

source("analysis_util.R")

options(mc.cores = max(1, parallel::detectCores()))
rstan_options(auto_write = TRUE)

load(file.path("data", "analysis.RData"))

# ATE ---------------------------------------------------------------------

dyn_sms_control_ate <- tribble(
  ~ private_value_left, ~ social_value_left, ~ private_value_right, ~ social_value_right,
  "control",            "ink",               "control",             "control",
  "calendar",           "control",           "control",             "control",
  "calendar",           "bracelet",          "control",             "control",
  "control",            "bracelet",           "control",            "control" 
) %>% 
  mutate(signal_observed_left = social_value_left,
         signal_observed_right = social_value_right) %>% 
  bind_rows(
    tribble(
      ~ private_value_left, ~ social_value_left, ~ signal_observed_left, ~ private_value_right, ~ social_value_right, ~ signal_observed_right,
      "control",            "bracelet",          "control",              "control",             "control",            "control",
      "control",            "ink",               "control",              "control",             "control",            "control",
      "control",            "bracelet",          "bracelet",             "control",             "bracelet",           "control",
      "control",            "ink",               "ink",                  "control",             "ink",                "control"
  )) %>%
  mutate(incentive_shift_left = private_value_left,
         incentive_shift_right = private_value_right) %>% 
  bind_rows(
    tribble(
      ~ private_value_left, ~ incentive_shift_left,  ~ private_value_right, ~ incentive_shift_right,
      "calendar",           "calendar",              "calendar",            "control" 
    ) %>% 
    mutate(social_value_left = "control",
           social_value_right = "control",
           signal_observed_left = "control",
           signal_observed_right = "control")) %>%
  mutate(dyn_dist_pot_left = NA,
         dyn_dist_pot_right = NA) %>% 
  bind_rows(
    tribble(
      ~ social_value_left, ~ dyn_dist_pot_left, ~ social_value_right, ~ dyn_dist_pot_right,
      "control",              "close",             "control",               "far",
      "bracelet",             "close",             "bracelet",              "far",
      "ink",                  "close",             "ink",                   "far"
    ) %>% 
      mutate(private_value_left = "control",
             private_value_right = "control",
             signal_observed_left = social_value_left,
             signal_observed_right = social_value_right,
             incentive_shift_left = private_value_left,
             incentive_shift_right = private_value_right) 
  ) %>% 
  bind_rows("close" = ., "far" = ., .id = "dist.pot.group") %>% 
  bind_rows(`TRUE` = ., `FALSE` = ., .id = "phone_owner") %>%
  bind_rows(`TRUE` = ., `FALSE` = ., .id = "name_matched") %>%
  mutate(sms.treatment.2_left = "control", 
         sms.treatment.2_right = "control",
         name_matched = FALSE,
         # hh.baseline.sample = FALSE,
         reminder_info_stock_left =  sms.treatment.2_left,
         reminder_info_stock_right = sms.treatment.2_right) %>% 
  mutate_at(vars(phone_owner, name_matched), as.logical) %>%
  mutate_at(vars(starts_with("dyn_dist_pot")), funs(coalesce(., dist.pot.group))) %>% 
  filter(!name_matched) %>% select(-name_matched)

dyn_phone_owners_ate <- tribble(
  ~ private_value_left, ~ social_value_left, ~ sms.treatment.2_left, ~ private_value_right, ~ social_value_right, ~ sms.treatment.2_right,
  "control",            "control",           "reminder.only",        "control",             "control",            "control",
  "control",            "control",           "social.info",          "control",             "control",            "control",
  "control",            "control",           "social.info",          "control",             "control",            "reminder.only",

  "control",            "ink",               "social.info",          "control",             "control",            "control",
  "control",            "ink",               "social.info",          "control",             "control",            "social.info",
  "control",            "ink",               "social.info",          "control",             "ink",                "control",

  "calendar",           "control",           "social.info",          "control",             "control",            "control",
  "calendar",           "control",           "social.info",          "control",             "control",            "social.info",
  "calendar",           "control",           "social.info",          "calendar",            "control",            "control",

  "calendar",           "bracelet",          "social.info",          "control",             "control",            "control",
  "calendar",           "bracelet",          "social.info",          "control",             "control",            "social.info",
  "calendar",           "bracelet",          "social.info",          "calendar",            "control",            "social.info",
  "calendar",           "bracelet",          "social.info",          "calendar",            "bracelet",           "control",
  
  "control",           "bracelet",           "social.info",          "control",             "control",            "control",
  "control",           "bracelet",           "social.info",          "control",             "control",            "social.info",
  "control",           "bracelet",           "social.info",          "calendar",            "control",            "social.info",
  "control",           "bracelet",           "social.info",          "control",             "bracelet",           "control"
) %>%
  mutate(signal_observed_left = social_value_left,
         signal_observed_right = social_value_right,
         reminder_info_stock_left = sms.treatment.2_left,
         reminder_info_stock_right = sms.treatment.2_right) %>% 
  bind_rows(
    tribble(
      ~ social_value_left, ~ social_value_right, ~ sms.treatment.2_right,
      "bracelet",          "bracelet",           "control",
      "bracelet",          "control",            "social.info",
      "ink",               "ink",                "control",    
      "ink",               "control",            "social.info"
    ) %>% 
    mutate(private_value_left = "control",
           private_value_right = "control",
           signal_observed_left = "control",
           signal_observed_right = "control",
           reminder_info_stock_left = "control",
           reminder_info_stock_right = "control",
           sms.treatment.2_left = "social.info")) %>% 
  bind_rows(
    tribble(
      ~ private_value_left, ~ private_value_right, ~ sms.treatment.2_right,
      "calendar",           "calendar",            "control",
      "calendar",           "control",             "social.info"
    ) %>% 
    mutate(social_value_left = "control",
           social_value_right = "control",
           signal_observed_left = "control",
           signal_observed_right = "control",
           reminder_info_stock_left = "control",
           reminder_info_stock_right = "control",
           sms.treatment.2_left = "social.info")) %>% 
  bind_rows(
    tribble(
      ~ social_value_left, ~ reminder_info_stock_left, ~ social_value_right,  ~ sms.treatment.2_right, 
      "bracelet",          "social.info",              "bracelet",            "social.info",   
      "bracelet",          "control",                  "bracelet",            "control",       
      "ink",               "social.info",              "ink",                 "social.info",   
      "ink",               "control",                  "ink",                 "control"       
  ) %>% 
    mutate(private_value_left = "control",
           private_value_right = "control",
           signal_observed_left = social_value_left,
           signal_observed_right = social_value_right,
           sms.treatment.2_left = "social.info",
           reminder_info_stock_right = "control")) %>% 
  bind_rows(
    tribble(
       ~ social_value_left, ~ private_value_left, ~ dyn_dist_pot_left, ~ dyn_dist_pot_right,
       "control",           "control",            "close",             "far",
       "control",           "calendar",           "close",             "far",
       "ink",               "control",            "close",             "far",
       "bracelet",          "control",            "close",             "far" 
    ) %>%
    mutate(
      private_value_right = private_value_left,
      social_value_right = social_value_left,
      signal_observed_left = "control",
      signal_observed_right = "control",
      sms.treatment.2_left = "social.info",
      sms.treatment.2_right = "social.info",
      reminder_info_stock_left = sms.treatment.2_left,
      reminder_info_stock_right = sms.treatment.2_right,
      incentive_shift_left = private_value_left,
      incentive_shift_right = private_value_right
    )
  ) %>% 
  bind_rows("close" = ., "far" = ., .id = "dist.pot.group") %>%
  mutate(
    phone_owner = TRUE,
    incentive_shift_left = "control",
    incentive_shift_right = "control",
    # hh.baseline.sample = FALSE, 
    name_matched = FALSE
  ) %>%
  mutate_at(vars(starts_with("dyn_dist_pot")), funs(coalesce(., dist.pot.group))) %>% 
  select(-name_matched)
  # bind_rows(filter(., sms.treatment.2_left == "control" & sms.treatment.2_right == "control") %>% mutate(name_matched = TRUE))

dyn_all_ate <- bind_rows(dyn_sms_control_ate, dyn_phone_owners_ate) 

assertthat::assert_that(all(map_lgl(dyn_all_ate, ~ !any(is.na(.)))))

# Analysis Data -----------------------------------------------------------


dyn_analysis_data <- analysis.data %>%
  mutate(private_value = fct_collapse(assigned.treatment, 
                                      control = c("control", "ink"), 
                                      "calendar" = c("calendar", "bracelet")), 
         social_value = fct_collapse(assigned.treatment, control = c("control", "calendar")),
         sms.treatment.2 = fct_recode(sms.treatment.2, control = "sms.control")) %>% 
  filter(!name_matched) #, sms.treatment.2 == "control") #, !hh.baseline.sample)

dyn_static_treatment_map <- dyn_analysis_data %>% 
  data_grid(private_value, social_value, sms.treatment.2, dist.pot.group, phone_owner) %>% 
  filter(sms.treatment.2 == "control" | phone_owner,
         sms.treatment.2 != "reminder.only" | (private_value == "control" & social_value == "control"),
         private_value == "control" | social_value != "ink") %>%
  prepare_treatment_map()

dyn_stan_data <- prepare_bayesian_analysis_data(
  dyn_analysis_data,
  wtp.data, 
   
  prepared_treatment_maps = TRUE, 
  treatment_map = dyn_static_treatment_map,
    
  # treatment_formula = ~ (private_value + social_value * dist.pot.group) * phone_owner,
  treatment_formula = ~ (private_value + social_value * dist.pot.group) * phone_owner + (social_value * dist.pot.group) : sms.treatment.2 + sms.treatment.2,
  subgroup_col = "phone_owner",
  drop_intercept_from_dm = FALSE, 
  nonparam_dynamics = TRUE,
  param_poly_order = 2,
  
  all_ate = dyn_all_ate,
    # filter(sms.treatment.2_left == "control",
    #        sms.treatment.2_right == "control") %>% 
    # select(-starts_with("sms.treatment"), -starts_with("reminder_info_stock")) %>% 
    # distinct(),
  
  scale_sigma = 1,
  hyper_coef_sigma = 1,
  hyper_intercept_sigma = 5,
  
  lkj_df = 2,
  
  use_logit = 1,
  
  estimate_ate = 1
)

# Run ---------------------------------------------------------------------

dyn_fit_version <- "param_sms_ctrl_6"

model_3_param <- stan_model(file = file.path("stan_models", "takeup_model_3_param.stan"), model_name = "model_3_param")

model_3_fit <- dyn_stan_data %>% 
  sampling(model_3_param, data = ., 
           chains = 4, iter = 600,
           control = lst(max_treedepth = 15), #adapt_delta = 0.99), 
           sample_file = file.path("stanfit", str_interp("model_3_${dyn_fit_version}.csv")))

save(param_dyn_stan_data, file = file.path("stan_analysis_data", str_interp("model_3_${dyn_fit_version}.RData")))