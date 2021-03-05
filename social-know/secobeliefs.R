library(magrittr)
library(tidyverse)
library(modelr)
library(cmdstanr)

load(file.path("data", "analysis.RData"))

chains <- 8

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

treatment_formula <- ~ assigned.treatment * dist.pot.group  

treatment_map <- social_links %>%
  distinct(!!!syms(all.vars(treatment_formula))) %>% 
  arrange(!!!syms(all.vars(treatment_formula))) %>% 
  mutate(treatment_id = seq_len(n()))

treatment_map_design_matrix <- treatment_map %>%
  model_matrix(treatment_formula)

ate_pairs <- treatment_map %>% {
  bind_rows(
    left_join(., filter(., fct_match(assigned.treatment, "control")), by = c("dist.pot.group"), suffix = c("", "_control")) %>% 
      filter(assigned.treatment != assigned.treatment_control) %>% 
      select(treatment_id, treatment_id_control),
    
    left_join(., filter(., fct_match(dist.pot.group, "close")), by = c("assigned.treatment"), suffix = c("", "_control")) %>% 
      filter(dist.pot.group != dist.pot.group_control) %>% 
      select(treatment_id, treatment_id_control),
  )
} %>%
  arrange(treatment_id, treatment_id_control) 
    

social_links %<>%
  left_join(treatment_map)

cluster_idx <- social_links %>% 
  group_by(cluster.id) %>% 
  summarize(cluster_index = cur_group_id())

social_links %<>% 
  left_join(cluster_idx, by = "cluster.id")

# Sample ------------------------------------------------------------------

secobeliefs_model <- cmdstan_model(file.path("stan_models", "secobeliefs.stan"))

secobeliefs_stan_data <- lst(
  use_stratum_level = TRUE,
  use_cluster_level = TRUE,
  use_obs_level = TRUE,
  
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

secobeliefs_fit <- secobeliefs_model$sample( 
  data = secobeliefs_stan_data,
  chains = chains,
  parallel_chains = chains,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.9
)

secobeliefs_fit$save_output_files(
  dir = file.path("~/Code/takeup", "data", "stan_analysis_data"), # Temp until cmdstanr is fixed and can handle paths with spaces
  basename = "secobeliefs_fit", 
  timestamp = FALSE, random = FALSE
)

secobeliefs_fit$save_object(file.path("data", "stan_analysis_data", "secobeliefs_fit.rds"))

# Results -----------------------------------------------------------------

secobeliefs_results <- lst(
  prob_knows = secobeliefs_fit$draws(c("prob_1ord", "prob_2ord")) %>% 
    posterior::as_draws_df() %>% 
    mutate(iter_id = .draw) %>% 
    pivot_longer(-c(iter_id, .draw, .chain, .iteration)) %>% 
    tidyr::extract(name, c("ord", "treatment_id"), r"{([12])ord\[(\d+)\]}", convert = TRUE) %>% 
    group_by(ord, treatment_id) %>% 
    summarize(
      per = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
      per_val = quantile(value, per),
      .groups = "drop"
    ) %>% 
    ungroup() %>% 
    pivot_wider(c(treatment_id, ord), names_from = per, values_from = per_val, names_prefix = "per_") %>% 
    right_join(treatment_map, ., by = "treatment_id") %>% 
    arrange(ord),
  
  ate_knows = secobeliefs_fit$draws(c("ate_1ord", "ate_2ord")) %>% 
    posterior::as_draws_df() %>% 
    mutate(iter_id = .draw) %>% 
    pivot_longer(-c(iter_id, .draw, .chain, .iteration)) %>% 
    tidyr::extract(name, c("ord", "ate_pair_index"), r"{([12])ord\[(\d+)\]}", convert = TRUE) %>% 
    group_by(ord, ate_pair_index) %>% 
    summarize(
      per = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
      per_val = quantile(value, per),
      .groups = "drop"
    ) %>% 
    ungroup() %>% 
    pivot_wider(c(ord, ate_pair_index), names_from = per, values_from = per_val, names_prefix = "per_") %>% 
    right_join(mutate(ate_pairs, ate_pair_index = seq(n())), ., by = "ate_pair_index") %>% 
    right_join(treatment_map, ., by = c("treatment_id")) %>%
    right_join(treatment_map, ., by = c("treatment_id" = "treatment_id_control"), suffix = c("_right", "_left")) %>%
    select(-starts_with("treatment_id"), -ate_pair_index) %>% 
    arrange(ord)
)

write_rds(secobeliefs_results, file.path("data", "stan_analysis_data", "secobeliefs_results.rds"))
  

