library(tidyverse)
library(cmdstanr)

params = list(structural_takeup_version = 66)



dist_fit_data <- str_split(params$structural_takeup_version, ",") %>% 
  pluck(1) %>% 
  as.integer() %>% 
  map_dfr(~ {
    temp_env <- new.env()
    load(file.path("temp-data", str_glue("processed_dist_fit{.x}.RData")), envir = temp_env)
    
    temp_env$dist_fit_data %>% 
      filter(fct_match(model_type, "structural")) %>% 
      mutate(version = .x)
  }) %>% 
  group_by(model, fit_type, model_type) %>% 
  filter(min_rank(version) == n()) %>% 
  ungroup()



collapse_estim = dist_fit_data %>% 
  filter(fct_match(model_type, "structural")) %>% 
  select(fit_type, obs_cluster_mu_rep)  %>%
  unnest(obs_cluster_mu_rep) %>%
  unnest(iter_data) %>%
  mutate(cluster_id = str_extract(variable, "\\d+") %>% as.numeric) %>%
  arrange(cluster_id)  %>%
  left_join(
    stan_data[[1]]$analysis_data %>%
        select(cluster_id,assigned_dist_group, assigned.treatment, cluster.dist.to.pot), 
        by = "cluster_id"
  )  %>%
  group_by(assigned.treatment, assigned_dist_group) %>%
  summarise(
    median = median(iter_est), 
    mean = mean(iter_est), 
    conf.low = quantile(iter_est, 0.25), 
    conf.high = quantile(iter_est, 0.75)) 

collapse_estim %>%
    ggplot(aes(
        x = median, 
        xmin = conf.low, 
        xmax = conf.high, 
        colour = assigned.treatment, 
        y = assigned_dist_group
    )) +
    geom_pointrange(position = position_dodge(0.5))  +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(title = "Estimated Mu Rep, collapsed by close-far", 
    x = "Estimate", 
    caption = "Fifty percent credible intervals")
ggsave("temp-data/mu-by-d-collapsed.png", width = 10, height = 10, dpi = 500)

dist_fit_data %>% 
  filter(fct_match(model_type, "structural")) %>% 
  select(fit_type, obs_cluster_mu_rep)  %>%
  unnest(obs_cluster_mu_rep) %>%
  unnest(iter_data) %>%
  group_by(variable) %>%
  summarise(
    median = median(iter_est), 
    mean = mean(iter_est), 
    conf.low = quantile(iter_est, 0.25), 
    conf.high = quantile(iter_est, 0.75)) %>%
  mutate(cluster_id = str_extract(variable, "\\d+") %>% as.numeric) %>%
  arrange(cluster_id) %>%
  left_join(
    stan_data[[1]]$analysis_data %>%
        select(cluster_id,assigned_dist_group, assigned.treatment, cluster.dist.to.pot), 
        by = "cluster_id"
  ) %>%
#   filter(assigned_dist_group == "close") %>%
  ggplot(aes(
    x = cluster.dist.to.pot, 
    y = median, 
    ymin = conf.low, 
    ymax = conf.high,
    colour  = assigned.treatment 
  )) +
  geom_pointrange() +
  facet_wrap(~assigned_dist_group, ncol = 1) +
  theme_bw() +
  labs(title = "Estimated Mu Rep, by cluster over distance", 
  caption = "fifty percent credible intervals (v skewed upwards)") +
  theme( 
    legend.position = "bottom"
  ) 


ggsave("temp-data/mu-by-d.png", width = 10, height = 10, dpi = 500)