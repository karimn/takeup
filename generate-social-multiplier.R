library(tidyverse)
library(ggthemes)


canva_palette_vibrant = "Primary colors with a vibrant twist"

params = list(fit_version = 86)
params$models = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP"

sm_df = map_dfr(
  c("bracelet", "calendar", "ink", "control"), 
  ~read_csv(
    str_glue(
      "optim/data/pred-social-multiplier-fit{params$fit_version}-cutoff-b-{.x}-mu-{.x}-{params$models}.csv"
    )
  ) %>% mutate(treatment = .x)
)


sm_df = sm_df %>%
    mutate(
        sm_fix_mu = (-delta + 0) / (1 + delta_v_star_deriv)
    )

sm_df = sm_df %>%
    mutate(
        sm = (-delta + mu_rep_deriv*delta_v_star) / (1 + mu_rep*delta_v_star_deriv),
        sm_delta_part = -delta / (1 + mu_rep*delta_v_star_deriv),
        sm_mu_deriv_part = mu_rep_deriv*delta_v_star  / (1 + mu_rep*delta_v_star_deriv),
        sm_fix_mu = (-delta) / (1 + mu_rep*delta_v_star_deriv)
    ) 
# TODO: Create Delta thing
# sm = (-delta + mu_rep_deriv*delta_v_star) / (1 + mu_rep*delta_v_star_deriv),
# -delta / (1 + mu Delta') + mu' Delta / (1 + mu Delta')




summ_sm_df = sm_df %>%
  mutate(
    across(matches("^sm"), ~.x/delta )
    ) %>%
  select(matches("^sm"), treatment, dist) %>%
  gather(variable, value, -treatment, -dist) %>%
  group_by(treatment, dist, variable) %>%
  summarise(
    mean_est = mean(value), 
    per_0.05 = quantile(value, 0.05),
    per_0.5 = quantile(value, 0.5),
    per_0.95 = quantile(value, 0.95),
    per_0.25 = quantile(value, 0.25),
    per_0.75 = quantile(value, 0.75)
  ) %>%
  mutate(treatment = factor(treatment, levels = c('control', "ink", 'calendar', 'bracelet')))

summ_sm_df

summ_sm_df %>%
  filter(variable == "sm") %>%
  filter(dist > 0) %>%
  mutate(dist = dist/1000) %>%
  mutate(
    treatment = factor(treatment, levels = c("control","bracelet", "ink", "calendar")), 
    treatment = fct_relabel(treatment, str_to_title)
  ) %>%
  ggplot(aes(
    x = dist, y = per_0.5
  )) +
  geom_line(aes(colour = treatment), size = 1) +
  labs(
    x = "Distance to Treatment (d) [km]" ,
     y = "Social Multiplier", 
     title = "Social Multiplier",
   colour = ""
  ) +
  geom_hline(
    yintercept = -1, 
    linetype = "longdash"
  ) + 
    guides(fill = "none") +
  ggthemes::scale_color_canva(
    palette = canva_palette_vibrant
  ) +
  annotate(
    "text", 
    x = 0.0 + 0.2, 
    y = -1 - 0.02,
    label = "Amplification", 
    size = 3, 
    alpha = 0.7
  ) +
  annotate(
    "text", 
    x = 2.5 , 
    y = -1 + 0.02,
    label = "Mitigation", 
    size = 3,
    alpha = 0.7
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")



summ_sm_df %>%
  filter(variable %in% c("sm", "sm_delta_part", "sm_mu_deriv_part")) %>%
  filter(dist > 0) %>%
  mutate(dist = dist/1000) %>%
  mutate(
    treatment = factor(treatment, levels = c("control","bracelet", "ink", "calendar")), 
    treatment = fct_relabel(treatment, str_to_title),
    type = factor(if_else(str_detect(variable, "part"), "part", "full"), levels = c("part", "full"))
  ) %>%
  ggplot(aes(
    x = dist, y = per_0.5,
    colour = treatment
  )) +
  geom_line(aes(linetype = variable), size = 1) +
  labs(
    x = "Distance to Treatment (d) [km]" ,
   y = "Social Multiplier", 
   title = "Decomposing into Multiplier Components (kinda)",
   colour = ""
  ) +
  geom_hline(
    yintercept = -1, 
    linetype = "longdash"
  ) + 
    guides(fill = "none") +
  ggthemes::scale_color_canva(
    "",
    palette = canva_palette_vibrant
  ) +
  annotate(
    "text", 
    x = 0.0 + 0.2, 
    y = -1 - 0.02,
    label = "Amplification", 
    size = 3, 
    alpha = 0.7
  ) +
  annotate(
    "text", 
    x = 2.5 , 
    y = -1 + 0.02,
    label = "Mitigation", 
    size = 3,
    alpha = 0.7
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")  +
  facet_wrap(~type, ncol = 2)



summ_sm_df %>%
  filter(variable == "sm_fix_mu") %>%
  filter(dist > 0) %>%
  mutate(dist = dist/1000) %>%
  mutate(
    treatment = factor(treatment, levels = c("control","bracelet", "ink", "calendar")), 
    treatment = fct_relabel(treatment, str_to_title)
  ) %>%
  ggplot(aes(
    x = dist, y = per_0.5,
    linetype = variable
  )) +
  geom_line(aes(colour = treatment), size = 1) +
  labs(
    x = "Distance to Treatment (d) [km]" ,
   y = "Social Multiplier", 
   title = "Setting Mu = 1, Mu' = 0",
   colour = ""
  ) +
  geom_hline(
    yintercept = -1, 
    linetype = "longdash"
  ) + 
    guides(fill = "none") +
  ggthemes::scale_color_canva(
    palette = canva_palette_vibrant
  ) +
  annotate(
    "text", 
    x = 0.0 + 0.2, 
    y = -1 - 0.02,
    label = "Amplification", 
    size = 3, 
    alpha = 0.7
  ) +
  annotate(
    "text", 
    x = 2.5 , 
    y = -1 + 0.02,
    label = "Mitigation", 
    size = 3,
    alpha = 0.7
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


## Endog v adjust


sm_endog_fix_df = map_dfr(
  c("bracelet", "calendar", "ink", "control"), 
  ~read_csv(
    str_glue(
      "optim/data/pred-social-multiplier-fix-mu-fit{params$fit_version}-cutoff-b-{.x}-mu-{.x}-{params$models}.csv"
    )
  ) %>% mutate(treatment = .x)
)


sm_endog_fix_df

summ_sm_endog_fix_df = sm_endog_fix_df %>%
  mutate(
        sm = (-delta + mu_rep_deriv*delta_v_star) / (1 + mu_rep*delta_v_star_deriv),
        sm = sm/delta
    ) %>%
  select(-draw) %>%
  gather(variable, value, -treatment, -dist) %>%
  group_by(treatment, dist, variable) %>%
  summarise(
    mean_est = mean(value), 
    per_0.05 = quantile(value, 0.05),
    per_0.5 = quantile(value, 0.5),
    per_0.95 = quantile(value, 0.95),
    per_0.25 = quantile(value, 0.25),
    per_0.75 = quantile(value, 0.75)
  ) %>%
  mutate(treatment = factor(treatment, levels = c('control', "ink", 'calendar', 'bracelet')))



summ_sm_endog_fix_df %>%
  filter(variable == "sm") %>%
  filter(dist > 0) %>%
  mutate(dist = dist/1000) %>%
  mutate(
    treatment = factor(treatment, levels = c("control","bracelet", "ink", "calendar")), 
    treatment = fct_relabel(treatment, str_to_title)
  ) %>%
  ggplot(aes(
    x = dist, y = per_0.5
  )) +
  geom_line(aes(colour = treatment), size = 1) +
  labs(
    x = "Distance to Treatment (d) [km]" ,
     y = "Social Multiplier", 
     title = "Social Multiplier",
   colour = ""
  ) +
  geom_hline(
    yintercept = -1, 
    linetype = "longdash"
  ) + 
    guides(fill = "none") +
  ggthemes::scale_color_canva(
    palette = canva_palette_vibrant
  ) +
  annotate(
    "text", 
    x = 0.0 + 0.2, 
    y = -1 - 0.02,
    label = "Amplification", 
    size = 3, 
    alpha = 0.7
  ) +
  annotate(
    "text", 
    x = 2.5 , 
    y = -1 + 0.02,
    label = "Mitigation", 
    size = 3,
    alpha = 0.7
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~treatment)




summ_sm_endog_fix_df %>%
  filter(dist > 0) %>%
  mutate(dist = dist/1000) %>%
  mutate(
    treatment = factor(treatment, levels = c("control","bracelet", "ink", "calendar")), 
    treatment = fct_relabel(treatment, str_to_title)
  ) %>%
  ggplot(aes(
    x = dist, y = per_0.5,
    colour = treatment
  )) +
  geom_point() +
  facet_wrap(~variable, scales = "free")


comp_sm_df = sm_df %>%
  select(-draw) %>%
  gather(variable, value, -treatment, -dist) %>%
  group_by(treatment, dist, variable) %>%
  summarise(
    mean_est = mean(value), 
    per_0.05 = quantile(value, 0.05),
    per_0.5 = quantile(value, 0.5),
    per_0.95 = quantile(value, 0.95),
    per_0.25 = quantile(value, 0.25),
    per_0.75 = quantile(value, 0.75)
  ) %>%
  mutate(treatment = factor(treatment, levels = c('control', "ink", 'calendar', 'bracelet')))



comp_sm_df %>%
  filter(dist > 0) %>%
  mutate(dist = dist/1000) %>%
  mutate(
    treatment = factor(treatment, levels = c("control","bracelet", "ink", "calendar")), 
    treatment = fct_relabel(treatment, str_to_title)
  ) %>%
  ggplot(aes(
    x = dist, y = per_0.5,
    colour = treatment
  )) +
  geom_point() +
  facet_wrap(~variable, scales = "free")
