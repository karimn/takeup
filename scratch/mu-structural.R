library(tidyverse)
library(magrittr)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

struct_fit = read_rds("data/stan_analysis_data/dist_fit66_STRUCTURAL_LINEAR_U_SHOCKS.rds")


struct_fit


load("data/stan_analysis_data/dist_fit66.RData")

struct_beliefs_results <- struct_fit$draws(c("prob_1ord", "prob_2ord", "ate_1ord", "ate_2ord")) %>% 
    posterior::as_draws_df() %>% 
    mutate(iter_id = .draw) %>% 
    pivot_longer(!c(iter_id, .draw, .iteration, .chain), names_to = "variable", values_to = "iter_est") %>% 
    nest(iter_data = !variable) %>% 
    get_beliefs_results(stan_data)

beliefs_results = read_rds("data/stan_analysis_data/beliefs_results.rds")



beliefs_results$ate_knows = bind_rows(
    beliefs_results$ate_knows %>% mutate(model = "only beliefs model"),
    struct_beliefs_results$ate_knows %>% mutate(model = "structural")
)

beliefs_results$prob = bind_rows(
    beliefs_results$prob %>% mutate(model = "only beliefs model"),
    struct_beliefs_results$prob %>% mutate(model = "structural")
)






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
    geom_pointrange(position = position_dodge(0.5)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    facet_grid(model~ord) +
    scale_colour_discrete(direction = -1) +
    labs(
        title = "ATE Beliefs Results - Estimated Separately vs in Joint Model", 
        x = "ATE", 
        y = "Treatment",
        caption = "Panel 1 indicates first order, panel 2 indicates second order beliefs"
    ) +
    geom_vline(xintercept = 0, linetype = "longdash") 


ggsave("temp-data/comp-belief-ate-66.png", width = 10, height = 10)
