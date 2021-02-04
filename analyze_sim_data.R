library(magrittr)
library(tidyverse)
library(cmdstanr)

# sim_data <- read_rds("temp-data/dist_sim.rds")
sim_data <- read_rds("/tigress/kn6838/takeup/dist_sim.rds")

sim_data %>% 
  select(iter_id, sim_takeup_levels) %>%
  mutate(
    fit_type = if_else(iter_id == 0, "prior", "sim") %>% factor()
  ) %>% 
  unnest(sim_takeup_levels) %>% 
  filter(
    (fct_match(assigned_treatment, "control") & fct_match(mu_assigned_treatment, "control")) | 
      (fct_match(assigned_treatment, "bracelet") & fct_match(mu_assigned_treatment, "bracelet")) |
      (fct_match(assigned_treatment, "control") & fct_match(mu_assigned_treatment, "bracelet"))) %>% 
  ggplot(aes(factor(iter_id))) +
  geom_linerange(aes(ymin = per_0.1, ymax = per_0.9, color = fit_type)) +
  geom_linerange(aes(ymin = per_0.05, ymax = per_0.95, color = fit_type), alpha = 0.5, size = 2) + 
                 # data = . %>% filter(fct_match(fit_type, "prior"))) +
  geom_point(aes(y = true_takeup_level, color = in_0.8_ci)) +
  scale_color_manual("", values = c(`FALSE` = "red", `TRUE` = "black", prior = "orange", sim = "black")) +
  labs(x = "", y = "") +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(rows = vars(assigned_dist_group), cols = vars(assigned_treatment, mu_assigned_treatment)) +
  theme(axis.text.x.bottom = element_blank(), axis.ticks = element_blank(), legend.position = "none") +
  NULL 

sim_data %>% 
  select(iter_id, sim_takeup_levels) %>% 
  unnest(sim_takeup_levels) %>% 
  group_by(assigned_treatment, mu_assigned_treatment) %>% 
  summarize(coverage = mean(in_0.8_ci))
