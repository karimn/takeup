library(tidyverse)
library(tidybayes)
library(cmdstanr)


beliefs_results = read_rds("data/stan_analysis_data/beliefs_results.rds")


beliefs_results$ate %>%
    filter(assigned_treatment_right == "control") %>%
    filter(assigned_dist_group_right == assigned_dist_group_left) %>%
    janitor::clean_names() %>%
    ggplot(aes(
        x = per_0_5,
        xmin = per_0_1,
        xmax = per_0_9,
        y = assigned_treatment_left, 
        colour = assigned_dist_group_left
    )) +
    geom_pointrange(size = 2, position = position_dodge(0.5)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_wrap(~ord, ncol = 1) +
    scale_colour_discrete(direction = -1) +
    labs(
        title = "ATE Beliefs Results - Estimated Separately", 
        x = "ATE", 
        y = "Treatment"
    ) +
    geom_vline(xintercept = 0, linetype = "longdash")
ggsave("temp-data/ATE-beliefs-estimates.png", width = 10, height = 10, dpi = 500)

beliefs_results$prob_knows %>%
    janitor::clean_names() %>%
    ggplot(aes(
        x = per_0_5,
        xmin = per_0_1,
        xmax = per_0_9,
        y = assigned_treatment, 
        colour = assigned_dist_group
    )) +
    geom_pointrange(size = 2, position = position_dodge(0.5)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_wrap(~ord, ncol = 1) +
    scale_colour_discrete(direction = -1) +
    labs(
        title = "Proportion Beliefs Results - Estimated Separately", 
        x = "Proportion", 
        y = "Treatment"
    )

ggsave("temp-data/prop-beliefs-estimates.png", width = 10, height = 10, dpi = 500)

## Data Prep

# source("analysis_util.R")
# source(file.path("multilvlr", "multilvlr_util.R"))
# source("dist_structural_util.R")

# # Data --------------------------------------------------------------------

# load(file.path("data", "analysis.RData"))

# standardize <- as_mapper(~ (.) / sd(.))
# unstandardize <- function(standardized, original) standardized * sd(original)

# monitored_nosms_data <- analysis.data %>% 
#   filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
#   left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
#             by = "cluster.id") %>% 
#   mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
#   group_by(cluster.id) %>% 
#   mutate(cluster_id = cur_group_id()) %>% 
#   ungroup()

# nosms_data <- analysis.data %>% 
#   filter(sms.treatment.2 == "sms.control") %>% 
#   left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
#             by = "cluster.id") %>% 
#   mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
#   group_by(cluster.id) %>% 
#   mutate(cluster_id = cur_group_id()) %>% 
#   ungroup()

# analysis_data <- monitored_nosms_data %>% 
#   mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group)


# ##


# beliefs_results$ate_knows  %>%
#     filter(assigned_treatment_right == "control") %>%
#     ggplot(aes(
#         x = `per_0.5`,
#         xmin = `per_0.1`,
#         xmax = `per_0.95`, 
#         y = 
#     )) +
#     facet_wrap(~)



# beliefs_files = fs::dir_ls("data/stan_analysis_data", regex = "beliefs_fit.*\\.csv")
# beliefs_fit =  as_cmdstan_fit(beliefs_files)


# beliefs_fit$summary()
# centered_obs_beta_2ord, centered_obs_dist_beta_2ord

# beta_belief_draws = beliefs_fit %>%
#     gather_draws(
#         ate_1ord[i],
#         ate_2ord[i]
#         )

# beta_belief_draws %>%
#     median_qi() %>%
#     to_broom_names() %>%
#     ggplot(aes(
#         x = estimate, 
#         xmin = conf.low, 
#         xmax = conf.high, 
#         y = factor(i), 
#         colour = term
#     ))  +
#     geom_pointrange(position = position_dodge(0.5))

# cluster_treatment_map = distinct(analysis_data, assigned_treatment, assigned_dist_group) %>% 
#   arrange(assigned_dist_group, assigned_treatment) # We must arrange by distance first

# beliefs_ate_pairs <- cluster_treatment_map %>% 
#   #   if (script_options$no_dist) {
#   #     distinct(., assigned_treatment)
#   #   } else .
#   # } %>%  
#   # filter(fct_match(assigned_dist_group, "close")) %>% 
#   mutate(treatment_id = seq(n())) %>% {
#     # if (script_options$no_dist) {
#     #   mutate(., treatment_id_control = 1) %>% 
#     #     filter(treatment_id != treatment_id_control) %>% 
#     #     select(treatment_id, treatment_id_control)
#     # } else {
#       bind_rows(
#         left_join(., filter(., fct_match(assigned_treatment, "control")), by = c("assigned_dist_group"), suffix = c("", "_control")) %>% 
#           filter(assigned_treatment != assigned_treatment_control) %>% 
#           select(assigned_treatment, assigned_dist_group, treatment_id, treatment_id_control),
        
#         left_join(., filter(., fct_match(assigned_dist_group, "close")), by = c("assigned_treatment"), suffix = c("", "_control")) %>% 
#           filter(assigned_dist_group != assigned_dist_group_control) %>% 
#           select(assigned_treatment, assigned_dist_group, treatment_id, treatment_id_control),
#       )
#     # }
# } %>%
#   arrange(treatment_id, treatment_id_control) 

# beliefs_ate_pairs


# beta_belief_draws
