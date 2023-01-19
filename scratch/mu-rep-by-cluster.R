library(tidyverse)
library(cmdstanr)

params = list(structural_takeup_version = 71)



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


load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)
# stick to monitored sms.treatment group
# remove sms.treatment.2
monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_nosms_data %>% 
mutate(
    assigned_treatment = assigned.treatment, 
    assigned_dist_group = dist.pot.group,
    sms_treatment = factor(sms.treatment.2))

collapse_estim = dist_fit_data %>% 
  filter(fct_match(model_type, "structural")) %>% 
  select(fit_type, obs_cluster_mu_rep)  %>%
  unnest(obs_cluster_mu_rep) %>%
  unnest(iter_data) %>%
  mutate(cluster_id = str_extract(variable, "\\d+") %>% as.numeric) %>%
  arrange(cluster_id)  %>%
  left_join(
    analysis_data %>%
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
    analysis_data %>%
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
dist_fit_data %>% 
  filter(fct_match(model_type, "structural"))  %>%
  filter(fit_type == "fit") %>%
  select(fit_type, model, obs_cluster_mu_rep)   %>%
  unnest(obs_cluster_mu_rep) %>%
  unnest(iter_data)  %>%
  group_by(variable, model) %>%
  summarise(
    median = median(iter_est), 
    mean = mean(iter_est), 
    conf.low = quantile(iter_est, 0.25), 
    conf.high = quantile(iter_est, 0.75)) %>%
  mutate(cluster_id = str_extract(variable, "\\d+") %>% as.numeric) %>%
  arrange(cluster_id) %>%
  left_join(
    analysis_data %>%
        select(cluster_id,assigned_dist_group, assigned.treatment, cluster.dist.to.pot), 
        by = "cluster_id"
  ) %>%
#   filter(assigned_dist_group == "close") %>%
  # filter(model == "STRUCTURAL_LINEAR_U_SHOCKS") %>%
  ggplot(aes(
    x = cluster.dist.to.pot, 
    y = median, 
    ymin = conf.low, 
    ymax = conf.high,
    colour  = assigned.treatment 
  )) +
  geom_point() +
  geom_line() +
  facet_wrap(~model, scales = "free", ncol = 1) +
  theme_bw() +
  labs(title = "Estimated Mu Rep, by cluster over distance", 
  caption = "fifty percent credible intervals (v skewed upwards)") +
  theme( 
    legend.position = "bottom"
  ) 

ggsave(
  "temp-data/mu-rep.png", width = 10, height = 10, dpi = 500)
