library(magrittr)
library(tidyverse)
library(modelr)
library(rstan)

load(file.path("data", "analysis.RData"))

chains <- 8

options(mc.cores = chains)
rstan_options(auto_write = TRUE)

# Prepare Data ------------------------------------------------------------

beliefs_data <- analysis.data %>%
  add_count(cluster.id, village, name = "village_size") %>%
  select(KEY.individ, county, village, KEY.individ, phone_owner, village_size) %>%
  inner_join(endline.know.table.data, by = "KEY.individ") %>%
  mutate(
    sms.treatment = factor(sms.treatment, levels = c("sms.control", "reminder.only", "social.info")) %>%
      fct_recode(control = "sms.control"),
    num.recognized = coalesce(num.recognized, 0L)
  ) %>% 
  filter(
    fct_match(know.table.type, "table.A"),
    # fct_match(sms.treatment, "control")
  ) 

social_links <- beliefs_data %>%
  group_by(KEY.individ, assigned.treatment, dist.pot.group, sms.treatment, phone_owner, county, cluster.id, village_size) %>%
  summarize(
    obs_know_person = sum(num.recognized),
    obs_know_person_prop = mean(num.recognized),
    knows_other_dewormed = sum(fct_match(dewormed, c("yes", "no")), na.rm = TRUE),
    knows_other_dewormed_yes = sum(fct_match(dewormed, "yes"), na.rm = TRUE),
    thinks_other_knows = sum(fct_match(second.order, c("yes", "no")), na.rm = TRUE),
    thinks_other_knows_yes = sum(fct_match(second.order, "yes"), na.rm = TRUE),

    .groups = "drop"
  ) %>%
  filter(obs_know_person > 0) %>% 
  mutate(obs_index = seq_len(n()))

treatment_formula <- ~ assigned.treatment  

treatment_map <- social_links %>%
  distinct(!!!syms(all.vars(treatment_formula))) %>% 
  arrange(!!!syms(all.vars(treatment_formula))) %>% 
  mutate(treatment_id = seq_len(n()))

control_treatment_map <- treatment_map %>%
  filter(fct_match(assigned.treatment, "control")) %>%
  select(-assigned.treatment) %>% 
  rename(treatment_id_control = treatment_id)

if (ncol(control_treatment_map) > 1) {
  treatment_map %<>%
    left_join(control_treatment_map)
} else {
  treatment_map %<>%
    bind_cols(control_treatment_map)
}

treatment_map_design_matrix <- treatment_map %>%
  model_matrix(treatment_formula)

ate_pairs <- treatment_map %>% 
  filter(treatment_id != treatment_id_control) %>% 
  select(treatment_id, treatment_id_control) %>% 
  arrange(treatment_id, treatment_id_control) 

social_links %<>%
  left_join(treatment_map)

# obs_missing_data <- treatment_map %>%
#   group_by(treatment_id) %>%
#   do(treatment_obs_ids = semi_join(social_links, ., "treatment_id"),
#      treatment_missing_ids = anti_join(social_links, ., "treatment_id")) %>%
#   ungroup() %>%
#   mutate(across(ends_with("_ids"), map, select, obs_index, obs_know_person, thinks_other_knows))

cluster_idx <- social_links %>% 
  group_by(cluster.id) %>% 
  summarize(cluster_index = cur_group_id())

social_links %<>% 
  left_join(cluster_idx, by = "cluster.id")

# Sample ------------------------------------------------------------------

secobeliefs_model <- stan_model(file = file.path("stan_models", "secobeliefs.stan"))

secobeliefs_stan_data <- lst(
  know_table_A_sample_size = 10,
  
  num_obs = nrow(social_links),
  num_clusters = n_distinct(social_links$cluster.id),
  num_strata = n_distinct(social_links$county),
  num_treatments = nrow(treatment_map),
  
  treatment_map_design_matrix,
  
  cluster_index = social_links$cluster_index,
  stratum_index = social_links$county %>% as.integer(),
  
  treatment_id = social_links$treatment_id,
  
  num_recognized = social_links$obs_know_person,
  num_knows_1ord = social_links$knows_other_dewormed,
  num_knows_2ord = social_links$thinks_other_knows,
  
  ate_pairs,
  num_ate_pairs = nrow(ate_pairs),
)

secobeliefs_fit <- sampling(
  secobeliefs_model, 
  data = secobeliefs_stan_data,
  chains = chains,
  iter = 1000,
  pars = c("prob_1ord", "ate_1ord", "prob_2ord", "ate_2ord"),
  control = lst(adapt_delta = 0.9),
)

# Results -----------------------------------------------------------------

secobeliefs_results <- lst(
  prob_2ord = secobeliefs_fit %>% 
    as.data.frame(pars = c("prob_2ord")) %>% 
    mutate(iter_id = seq(n())) %>% 
    pivot_longer(-iter_id) %>% 
    tidyr::extract(name, "treatment_id", r"{\[(\d+)\]}", convert = TRUE) %>% 
    group_by(treatment_id) %>% 
    summarize(
      per = c(0.1, 0.25, 0.5, 0.75, 0.9),
      per_val = quantile(value, per)
    ) %>% 
    ungroup() %>% 
    pivot_wider(treatment_id, names_from = per, values_from = per_val, names_prefix = "per_") %>% 
    right_join(treatment_map, ., by = "treatment_id"),
  
  ate_2ord = secobeliefs_fit %>% 
    as.data.frame(pars = c("ate_2ord")) %>% 
    mutate(iter_id = seq(n())) %>% 
    pivot_longer(-iter_id) %>% 
    tidyr::extract(name, "ate_pair_index", r"{\[(\d+)\]}", convert = TRUE) %>% 
    group_by(ate_pair_index) %>% 
    summarize(
      per = c(0.1, 0.25, 0.5, 0.75, 0.9),
      per_val = quantile(value, per)
    ) %>% 
    ungroup() %>% 
    pivot_wider(ate_pair_index, names_from = per, values_from = per_val, names_prefix = "per_") %>% 
    right_join(mutate(ate_pairs, ate_pair_index = seq(n())), ., by = "ate_pair_index") %>% 
    right_join(treatment_map, ., by = c("treatment_id", "treatment_id_control")) %>% 
    right_join(treatment_map, ., by = c("treatment_id" = "treatment_id_control"), suffix = c("_left", "_right")) %>% 
    select(-starts_with("treatment_id"), -ate_pair_index)
)

write_rds(secobeliefs_results, file.path("data", "stan_analysis_data", "secobeliefs_results.rds"))
  

