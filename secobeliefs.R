library(magrittr)
library(plyr)
library(forcats)
library(broom)
library(tidyverse)
library(modelr)
library(rstan)

load(file.path("data", "analysis.RData"))

options(mc.cores = 4)
rstan_options(auto_write = TRUE)

# Prepare Data ------------------------------------------------------------

beliefs_data <- analysis.data %>% 
  add_count(cluster.id, village, name = "village_size") %>% 
  select(KEY.individ, county, village, KEY.individ, phone_owner, village_size) %>% 
  inner_join(endline.know.table.data, "KEY.individ") %>% 
  filter(know.table.type == "table.A") %>% 
  mutate(
    sms.treatment = factor(sms.treatment, levels = c("sms.control", "reminder.only", "social.info")) %>% 
      fct_recode(control = "sms.control"),
    num.recognized = coalesce(num.recognized, 0L)
  )

social_links <- beliefs_data %>% 
  group_by(KEY.individ, assigned.treatment, dist.pot.group, sms.treatment, phone_owner, village_size) %>% 
  summarize(
    obs_know_person = sum(num.recognized),
    obs_know_person_prop = mean(num.recognized),
    thinks_other_knows = sum(second.order %in% c("yes", "no"), na.rm = TRUE),
    thinks_other_knows_yes = sum(second.order == "yes", na.rm = TRUE),
    
    .groups = "drop"
  ) %>% 
  mutate(obs_index = seq_len(n())) 

treatment_map <- social_links %>% 
  distinct(assigned.treatment, dist.pot.group, phone_owner, sms.treatment) %>% 
  mutate(treatment_id = seq_len(n()))

control_treatment_map <- treatment_map %>% 
  filter(assigned.treatment == "control")

treatment_map %<>% 
  left_join(control_treatment_map, c("phone_owner", "sms.treatment", "dist.pot.group"), suffix = c("", "_control"))

treatment_map_design_matrix <- treatment_map %>% 
  model_matrix(~ assigned.treatment * dist.pot.group * (phone_owner + sms.treatment)) %>% 
  magrittr::extract(, -1)

social_links %<>% 
  left_join(treatment_map, by = c("assigned.treatment", "dist.pot.group", "sms.treatment", "phone_owner")) 

obs_missing_data <- treatment_map %>% 
  group_by(treatment_id) %>% 
  do(treatment_obs_ids = semi_join(social_links, ., "treatment_id"),
     treatment_missing_ids = anti_join(semi_join(social_links, ., "phone_owner"), ., "treatment_id")) %>% 
  ungroup()

# Sample ------------------------------------------------------------------

secobeliefs_model <- stan_model(file = file.path("stan_models", "secobeliefs.stan"))

secobeliefs_fit <- 
  sampling(secobeliefs_model, 
           data = lst(
             num_obs = nrow(social_links),
             num_treatments = nrow(treatment_map),
             num_with_links_obs = social_links %>% filter(obs_know_person > 0) %>% nrow(),
             num_thinks_other_knows = social_links %>% filter(thinks_other_knows > 0) %>% nrow(),
             
             treatment_ids = social_links$treatment_id,
             control_treatment_ids = social_links$treatment_id_control,
             
             control_treatment_map = treatment_map$treatment_id_control,
             
             treatment_obs_ids = obs_missing_data %>% unnest(treatment_obs_ids) %>% pull(obs_index),
             treatment_missing_ids = obs_missing_data %>% unnest(treatment_missing_ids) %>% pull(obs_index),
             
             treatment_obs_sizes = obs_missing_data %$% map_int(treatment_obs_ids, nrow),
             treatment_missing_sizes = obs_missing_data %$% map_int(treatment_missing_ids, nrow),
             
             know_treatment_obs_ids_data = obs_missing_data %>% unnest(treatment_obs_ids) %>% filter(obs_know_person > 0) %>% select(treatment_id, obs_index),
             know_treatment_missing_ids_data = obs_missing_data %>% unnest(treatment_missing_ids) %>% filter(obs_know_person > 0) %>% select(treatment_id, obs_index),
             
             know_treatment_obs_sizes = know_treatment_obs_ids_data %>% count(treatment_id) %>% arrange(treatment_id) %>% pull(n),
             know_treatment_missing_sizes = know_treatment_missing_ids_data %>% count(treatment_id) %>% arrange(treatment_id) %>% pull(n),
             
             know_treatment_obs_ids = know_treatment_obs_ids_data %>% pull(obs_index),
             know_treatment_missing_ids = know_treatment_missing_ids_data %>% pull(obs_index),
             
             other_know_treatment_obs_ids_data = obs_missing_data %>% unnest(treatment_obs_ids) %>% filter(thinks_other_knows > 0) %>% select(treatment_id, obs_index),
             other_know_treatment_missing_ids_data = obs_missing_data %>% unnest(treatment_missing_ids) %>% filter(thinks_other_knows > 0) %>% select(treatment_id, obs_index),
             other_know_treatment_obs_sizes = other_know_treatment_obs_ids_data %>% count(treatment_id) %>% arrange(treatment_id) %>% pull(n),
             other_know_treatment_missing_sizes = other_know_treatment_missing_ids_data %>% count(treatment_id) %>% arrange(treatment_id) %>% pull(n),
             
             other_know_treatment_obs_ids = other_know_treatment_obs_ids_data %>% pull(obs_index),
             other_know_treatment_missing_ids = other_know_treatment_missing_ids_data %>% pull(obs_index),
             
             with_links_ids = social_links %>% filter(obs_know_person > 0) %>% pull(obs_index),
             thinks_other_knows_ids = social_links %>% filter(thinks_other_knows > 0) %>% pull(obs_index),
             
             obs_know_person = social_links$obs_know_person,
             
             village_size = social_links$village_size,
             
             num_treatment_coef = ncol(treatment_map_design_matrix),
             treatment_map_design_matrix,
               
             thinks_other_knows = social_links %>% pull(thinks_other_knows),
             thinks_other_knows_yes = social_links %>% pull(thinks_other_knows_yes)
           ), 
           chains = 4,
           iter = 600, 
           # control = lst(max_treedepth = , 
           #               adapt_delta = ), 
           sample_file = file.path("data", "stanfit", "secobeliefs_fit.csv"))

