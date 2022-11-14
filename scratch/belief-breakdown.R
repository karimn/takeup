## 
## Misc plots for Karing
##
library(tidyverse)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

# Data --------------------------------------------------------------------

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_nosms_data %>% 
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group)


# WTP Stan Data -----------------------------------------------------------

wtp_stan_data <- analysis.data %>% 
  mutate(stratum = county) %>% 
  prepare_bayes_wtp_data(
    wtp.data,
    
    preference_value_diff = seq(-100, 100, 10), 
    num_preference_value_diff = length(preference_value_diff), 
    
    wtp_utility_df = 3,
    tau_mu_wtp_diff = 100,
    mu_wtp_df_student_t = 7,
    tau_sigma_wtp_diff = 50,
    sigma_wtp_df_student_t = 2.5
  )

# Treatment Details -------------------------------------------------------

treatment_formula <- ~ assigned_treatment * assigned_dist_group 

cluster_treatment_map = distinct(analysis_data, assigned_treatment, assigned_dist_group) %>% 
  arrange(assigned_dist_group, assigned_treatment) # We must arrange by distance first

treatment_map_design_matrix <- cluster_treatment_map %>%
  modelr::model_matrix(treatment_formula)

# Beliefs Data ------------------------------------------------------------

beliefs_treatment_formula <- ~ assigned_treatment 

beliefs_treatment_map_design_matrix <- cluster_treatment_map %>%
  modelr::model_matrix(beliefs_treatment_formula) %>% 
  distinct()

analysis_data %<>% 
  nest_join(
    endline.know.table.data %>% 
      filter(fct_match(know.table.type, "table.A")),
    by = "KEY.individ", 
    name = "knowledge_data"
  ) %>% 
  mutate(
    map_dfr(knowledge_data, ~ {
      tibble(
        obs_know_person = sum(.x$num.recognized),
        obs_know_person_prop = mean(.x$num.recognized),
        knows_other_dewormed = sum(fct_match(.x$dewormed, c("yes", "no")), na.rm = TRUE),
        knows_other_dewormed_yes = sum(fct_match(.x$dewormed, "yes"), na.rm = TRUE),
        knows_other_dewormed_no = sum(fct_match(.x$dewormed, "no"), na.rm = TRUE),
        thinks_other_knows = sum(fct_match(.x$second.order, c("yes", "no")), na.rm = TRUE),
        thinks_other_knows_yes = sum(fct_match(.x$second.order, "yes"), na.rm = TRUE),
        thinks_other_knows_no = sum(fct_match(.x$second.order, "no"), na.rm = TRUE),
      )
    }
  ))



analysis_data %>%
    filter(obs_know_person > 0)  %>%
    select(KEY.individ, contains("know"), assigned.treatment, dist.pot.group, assigned_dist_group) %>%
    mutate(
        doesnt_know_other_dewormed = obs_know_person - knows_other_dewormed, 
        doesnt_think_other_knows = obs_know_person - thinks_other_knows
    ) %>% 
    select(KEY.individ, 
           assigned.treatment,
           assigned_dist_group,
           obs_know_person,
           knows_other_dewormed_yes,
           knows_other_dewormed_no, 
           doesnt_know_other_dewormed
           ) %>%
    gather(variable, value, 
        knows_other_dewormed_yes:doesnt_know_other_dewormed) %>%
    mutate(prop = value/obs_know_person) %>%
    filter(prop > 1)
stop()
comp_belief_data = analysis_data %>%
    filter(obs_know_person > 0)  %>%
    select(KEY.individ, contains("know"), assigned.treatment, dist.pot.group, assigned_dist_group) %>%
    mutate(
        doesnt_know_other_dewormed = obs_know_person - knows_other_dewormed, 
        doesnt_think_other_knows = obs_know_person - thinks_other_knows
    ) %>% 
    select(KEY.individ, 
           assigned.treatment,
           assigned_dist_group,
           obs_know_person,
           knows_other_dewormed_yes,
           knows_other_dewormed_no,
           doesnt_know_other_dewormed, 
           thinks_other_knows_yes, 
           thinks_other_knows_no, 
           doesnt_think_other_knows
           ) %>%
    gather(variable, value, 
        knows_other_dewormed_yes:doesnt_think_other_knows)   %>%
    mutate(knowledge_type = case_when(
        str_detect(variable, "_yes") ~ "yes",
        str_detect(variable, "_no") ~ "no",
        str_detect(variable, "doesnt") ~ "doesn't know"
    )) %>%
    mutate(belief_type = if_else(str_detect(variable, "think"), "2ord", "1ord")) %>%
    mutate(prop = value/obs_know_person) %>%
    group_by(assigned.treatment, assigned_dist_group, knowledge_type, belief_type) %>% 
    summarise(
        mean_est = mean(prop), 
        std.error = sd(prop)/sqrt(n()),
        per_0.5 = median(prop), 
        per_0.05 = per_0.5 - cv(0.05)*std.error, 
        per_0.95 = per_0.5 + cv(0.05)*std.error,
        per_0.1 =  per_0.5 - cv(0.1)*std.error,
        per_0.9 =  per_0.5 + cv(0.1)*std.error,
        per_0.25 =  per_0.5 - cv(0.25)*std.error,
        per_0.75 =  per_0.5 + cv(0.25)*std.error
    )

cv = function(alpha){
  qnorm(1 - (alpha/2))
}
cv(0.05)
comp_belief_data %>%
    # filter(assigned_dist_group == "close") %>%
    ggplot(
        aes(
            x = estimate, 
            # xmin = conf.low, 
            # xmax = conf.high, 
            y = knowledge_type, 
            colour = assigned.treatment,
            fill = assigned.treatment
        )
    ) +
    geom_col(colour = "black", position = position_dodge(0.8)) +
    facet_grid(assigned_dist_group ~ belief_type) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "Proportion", title = "Beliefs")

ggsave("temp-data/prop-belief.png", width = 10, height = 10, dpi = 600)    



comp_belief_data %>%
  filter(belief_type == "1ord") %>%
  filter(assigned_dist_group == "close") %>%
  ggplot(aes(
    x = estimate, 
    xmin = conf.low,
    xmax = conf.high,
    y = assigned.treatment, 
    colour = knowledge_type
  )) +
  geom_pointrange(position = position_dodge(0.5)) +
  theme_minimal() 


pos_dodge = position_dodge(width = 0.3)
crossbar_width = 0.2

comp_belief_data %>%
  filter(belief_type == "1ord") %>%
  filter(assigned_dist_group == "close") %>%
  ggplot(aes(y = assigned.treatment, group = knowledge_type)) +
  geom_linerange(aes(xmin = per_0.05, xmax = per_0.95, color = knowledge_type), position = pos_dodge, size = 0.3) +
  geom_crossbar(aes(x = per_0.5, xmin = per_0.1, xmax = per_0.9, color = knowledge_type), position = pos_dodge, fatten = 2, size = 0.4, width = crossbar_width) +
  geom_linerange(aes(xmin = per_0.25, xmax = per_0.75, color = knowledge_type), position = pos_dodge, alpha = 0.4, size = 2.25) +
  geom_point(aes(x = per_0.5, color = knowledge_type), position = pos_dodge, size = 1.8) +
  geom_point(aes(x = per_0.5), position = pos_dodge, color = "white", size = 0.6) +
  scale_y_discrete(drop = FALSE) +
  scale_color_canva("", labels = "FOB", palette = canva_palette_vibrant) + 
  labs(
    title = "FOB",
    subtitle = "",
    x = "", y = "") +
  theme(
    legend.position = "bottom"
  ) + 
  NULL


plot_single_beliefs_est <- function(beliefs_results_type_df, 
                                    order, 
                                    top_title = NULL, 
                                    width = 0.3, 
                                    crossbar_width = 0.2, 
                                    vline = TRUE) {
  pos_dodge <- position_dodge(width = width)
  if (order == 1) {
    str_title = "First Order Beliefs"
  } else {
    str_title = "Second Order Beliefs"
  }
  belief_plot = 
      beliefs_results_type_df %>% 
        filter(ord == order, assigned_dist_group_left == assigned_dist_group_right) %>% 
        ggplot(aes(y = assigned_treatment_left, group = assigned_dist_group_left)) +
        geom_linerange(aes(xmin = per_0.05, xmax = per_0.95, color = assigned_dist_group_left), position = pos_dodge, size = 0.3) +
        geom_crossbar(aes(x = per_0.5, xmin = per_0.1, xmax = per_0.9, color = assigned_dist_group_left), position = pos_dodge, fatten = 2, size = 0.4, width = crossbar_width) +
        geom_linerange(aes(xmin = per_0.25, xmax = per_0.75, color = assigned_dist_group_left), position = pos_dodge, alpha = 0.4, size = 2.25) +
        geom_point(aes(x = per_0.5, color = assigned_dist_group_left), position = pos_dodge, size = 1.8) +
        geom_point(aes(x = per_0.5), position = pos_dodge, color = "white", size = 0.6) +
        scale_y_discrete(drop = FALSE) +
        scale_color_canva("", labels = str_to_title, palette = canva_palette_vibrant) + 
        labs(
          title = str_title,
          subtitle = "",
          x = "", y = "") +
        theme(
          legend.position = "bottom"
        ) + 
        NULL

  if (vline == TRUE) {
    belief_plot = belief_plot +
         geom_vline(xintercept = 0, linetype = "dotted") 
  }
        
  return(belief_plot)
} 