
library(magrittr)
library(tidyverse)
library(broom)
library(ggrepel)
library(ggmap)
library(ggstance)
library(gridExtra)
library(cowplot)
library(rgeos)
library(sp)
library(knitr)
library(modelr)
library(car)
library(rstan)
library(latex2exp)
library(ggthemes)

library(econometr)
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))


source(file.path("multilvlr", "multilvlr_util.R"))

fit_version <- 66
default_top_levels = c("Bracelet", "Combined")

# 66 ed fit
# 60 Karim fit
# 62 also Karim fit


model_fit_by = if_else(fit_version %in% c(60, 62), "Karim", "Ed")


models_we_want = c(
  "STRUCTURAL_LINEAR_U_SHOCKS"
)


quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

output_basepath = str_glue("temp-data/output_dist_fit{fit_version}")

canva_palette_vibrant <- "Primary colors with a vibrant twist"

theme_set(theme_minimal() +
            theme(legend.position = "bottom"))

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

rct.schools.data <- read_rds(file.path("data", "takeup_rct_schools.rds"))
rct.cluster.selection <- read_rds(file.path("data", "rct_cluster_selection_2.0.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))

load(file.path("data", "takeup_village_pot_dist.RData"))

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_nosms_data


## Fit Loading

load(file.path("temp-data", str_interp("processed_dist_fit${fit_version}_lite.RData")))

dist_fit_data %<>%
  left_join(tribble(
    ~ fit_type,        ~ model_color,
      "fit",           "black", 
      "prior-predict", "darkgrey",
  ), by = "fit_type")

rf_analysis_data <- dist_fit_data %>% 
  filter(
    fct_match(model_type, "reduced form"),
    fct_match(fit_type, "fit"),
  ) %$% 
  stan_data[[1]]$analysis_data 
delta <- function(v, ...) dnorm(v, ...) / ((pnorm(v, ...) * pnorm(v, ..., lower.tail = FALSE)))

belief_data = dist_fit_data %>%
  filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) %>%
  pull(beliefs_results) %>%
  first() 

#### Separate Beliefs ####
belief_data$ate_knows %>%
  filter(assigned_treatment_left != "control") %>%
  mutate(assigned_treatment_left = fct_drop(assigned_treatment_left)) %>%
  plot_single_beliefs_est(
    width = 0.7, 
    order = 1,
    crossbar_width = 0.4) +
    labs(subtitle = "Treatment Effects") + 
    NULL

ggsave(file.path(output_basepath, str_glue("belief-te-1ord-plots.png")), width = 7.5, height = 5.0, dpi = 500)

belief_data$ate_knows %>%
  filter(assigned_treatment_left != "control") %>%
  mutate(assigned_treatment_left = fct_drop(assigned_treatment_left)) %>%
  plot_single_beliefs_est(
    width = 0.7, 
    order = 2,
    crossbar_width = 0.4) +
    labs(subtitle = "Treatment Effects") +
    NULL

ggsave(file.path(output_basepath, str_glue("belief-te-2ord-plots.png")), width = 7.5, height = 5.0, dpi = 500)

## Props

belief_data$prob_knows %>%
    mutate(assigned_treatment = fct_drop(assigned_treatment))  %>%
    mutate(assigned_treatment_left = assigned_treatment, assigned_treatment_right = assigned_treatment) %>%
    mutate(assigned_dist_group_left = assigned_dist_group, assigned_dist_group_right = assigned_dist_group) %>%
  plot_single_beliefs_est(
    width = 0.7, 
    order = 1,
    crossbar_width = 0.4, 
    vline = FALSE) +
    labs(subtitle = "Proportion") +
    NULL

ggsave(file.path(output_basepath, str_glue("belief-prop-1ord-plots.png")), width = 7.5, height = 5.0, dpi = 500)

belief_data$prob_knows %>%
    mutate(assigned_treatment = fct_drop(assigned_treatment)) %>%
    mutate(assigned_treatment_left = assigned_treatment, assigned_treatment_right = assigned_treatment) %>%
    mutate(assigned_dist_group_left = assigned_dist_group, assigned_dist_group_right = assigned_dist_group) %>%
    plot_single_beliefs_est(
    width = 0.7, 
    order = 2,
    crossbar_width = 0.4, 
    vline = FALSE) +
    labs(subtitle = "Proportion") +
    NULL

ggsave(file.path(output_basepath, str_glue("belief-prop-2ord-plots.png")), width = 7.5, height = 5.0, dpi = 500)

#### Belief ATEs combined ####
belief_data_control = belief_data$prob_knows %>%
    filter(assigned_treatment == "control")  %>%
    select(assigned_treatment, assigned_dist_group, iter_data, ord) %>%
    unnest(iter_data)

combined_belief_ate = belief_data$prob_knows %>%
    filter(assigned_treatment != "control") %>%
    select(assigned_treatment, iter_data, ord) %>%
    unnest(iter_data) %>%
    rename(assigned_treatment_left = assigned_treatment) %>%
    left_join(
        belief_data_control %>% 
            rename(assigned_treatment_right = assigned_treatment) %>%
            select(iter_est_right = iter_est, iter_id, ord, assigned_treatment_right),
        by = c("iter_id", "ord")
    ) %>%
    mutate(
        iter_est_te = iter_est - iter_est_right
    )  %>%
    select(-.chain, -.iteration, -.draw) %>%
    nest(iter_data = c(iter_id, iter_est_te, iter_est, iter_est_right)) %>%
    mutate(
        mean_est = map_dbl(iter_data, ~mean(.x$iter_est_te)), 
        quants = map(iter_data, quantilize_est, iter_est_te, quant_probs = quant_probs, na.rm = TRUE)
    )  %>%
    unnest(quants) %>%
    mutate(assigned_dist_group_left = "Combined", assigned_dist_group_right = "Combined") %>%
    mutate(
        assigned_treatment_left = fct_drop(assigned_treatment_left), 
        assigned_treatment_right = fct_drop(assigned_treatment_right)
    )  %>%
    mutate(assigned_treatment_left = fct_relabel(assigned_treatment_left, str_to_title)) 
combined_belief_prop = belief_data$prob_knows %>%
    select(assigned_treatment, iter_data, ord) %>%
    unnest(iter_data) %>%
    rename(assigned_treatment_left = assigned_treatment) %>%
    select(-.chain, -.iteration, -.draw) %>%
    nest(iter_data = c(iter_id, iter_est)) %>%
    mutate(
        mean_est = map_dbl(iter_data, ~mean(.x$iter_est)), 
        quants = map(iter_data, quantilize_est, iter_est, quant_probs = quant_probs, na.rm = TRUE)
    )  %>%
   unnest(quants) %>%
    mutate(assigned_dist_group_left = "Combined", assigned_dist_group_right = "Combined") %>%
    mutate(
        assigned_treatment_right = fct_drop(assigned_treatment_left)
    )  %>%
    mutate(assigned_treatment_left = fct_relabel(assigned_treatment_left, str_to_title))

combined_belief_prop %>%
  plot_single_beliefs_est(
    width = 0.7, 
    order = 1,
    crossbar_width = 0.4, 
    vline = FALSE) +
    labs(subtitle = "Proportion") +
    guides(colour = "none") +
    NULL 

ggsave(
    file.path(output_basepath, str_glue("combined-prop-fob.png")), 
    width = 7.5, 
    height = 5.0, 
    dpi = 500
    )

combined_belief_prop %>%
  plot_single_beliefs_est(
    width = 0.7, 
    order = 2,
    crossbar_width = 0.4, 
    vline = FALSE) +
    labs(subtitle = "Proportion") +
    guides(colour = "none") +
    NULL

ggsave(
    file.path(output_basepath, str_glue("combined-prop-sob.png")), 
    width = 7.5, 
    height = 5.0, 
    dpi = 500
    )


combined_belief_ate %>%
    plot_single_beliefs_est(order = 1)  +
    labs(subtitle = "Treatment Effects") +
    guides(colour = "none") 
ggsave(
    file.path(output_basepath, str_glue("combined-ate-fob.png")), 
    width = 7.5, 
    height = 5.0, 
    dpi = 500
    )

combined_belief_ate %>%
    plot_single_beliefs_est(order = 2)  +
    labs(subtitle = "Treatment Effects") +
    guides(colour = "none")
ggsave(
    file.path(output_basepath, str_glue("combined-ate-sob.png")), 
    width = 7.5, 
    height = 5.0, 
    dpi = 500
    )

#### Disaggregated FOBs ####

cv = function(alpha){
  qnorm(1 - (alpha/2))
}

comp_belief_data = analysis_data %>%
  mutate(assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group) %>%
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
  )) %>%
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
        per_0.5 = mean(prop), 
        per_0.05 = per_0.5 - cv(0.05)*std.error, 
        per_0.95 = per_0.5 + cv(0.05)*std.error,
        per_0.1 =  per_0.5 - cv(0.1)*std.error,
        per_0.9 =  per_0.5 + cv(0.1)*std.error,
        per_0.25 =  per_0.5 - cv(0.5)*std.error,
        per_0.75 =  per_0.5 + cv(0.5)*std.error
    ) %>%
    mutate(
      knowledge_type = factor(knowledge_type, levels = c("yes", "no", "doesn't know")), 
      knowledge_type = fct_relabel(knowledge_type, str_to_title) %>% fct_rev
    ) %>%
    mutate(
      assigned.treatment = factor(assigned.treatment, 
                                  levels = c("bracelet",
                                             "calendar", 
                                             "ink",
                                             "control")) %>%
                          fct_relabel(str_to_title) %>%
                          fct_rev
    )



comp_belief_data %>%
  filter(belief_type == "1ord") %>%
  filter(assigned_dist_group == "close") %>%
  plot_belief_breakdown() +
  NULL
  # labs(
  #   title = "Disaggregated First Order Beliefs - Close"
  # )

ggsave(file.path(
  output_basepath,
  "disagg-fob-close.png"
  ),
  width = 7.5,
  height = 5.0,
  dpi = 500)



comp_belief_data %>%
  filter(belief_type == "1ord") %>%
  filter(assigned_dist_group == "far") %>%
  plot_belief_breakdown() +
  NULL
  # labs(
  #   title = "Disaggregated First Order Beliefs - Far"
  # )
  

ggsave(file.path(
  output_basepath,
  "disagg-fob-far.png"
  ),
  width = 7.5,
  height = 5.0,
  dpi = 500)


#### Rate of Change
roc_df = dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
  transmute(
    model_name,
    cluster_roc = map2(cluster_roc, stan_data, ~ { 
      mutate(.x,
        roc_distance = roc_distance / 1000,
        across(starts_with("per_"), divide_by, sd(.y$analysis_data$cluster.dist.to.pot)), 
        across(starts_with("per_"), multiply_by, 1000),
        across(starts_with("per_"), multiply_by, 100)
      ) 
    })
  ) %>%
  unnest(cluster_roc) %>% 
  filter(fct_match(assigned_treatment, c("control", "ink", "bracelet", "calendar"))) 



roc_plot = function(roc_df, treatment) {
    plot = roc_df %>%
        filter(assigned_treatment %in% c(treatment, "control")) %>%
        mutate(assigned_treatment = fct_relabel(assigned_treatment, str_to_title)) %>%
        ggplot(aes(roc_distance)) +
        geom_line(aes(y = per_0.5, color = assigned_treatment)) +
        geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = assigned_treatment), alpha = 0.4) +
        geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = assigned_treatment), alpha = 0.4) +
        geom_rug(
          aes(dist, color = assigned_treatment),
          alpha = 0.75,
          data = rf_analysis_data %>%
            filter(fct_match(assigned.treatment, c("control", treatment))) %>%
            mutate(assigned.treatment = fct_relabel(assigned.treatment, str_to_title)) %>%
            distinct(cluster_id, assigned_treatment = assigned.treatment, dist = cluster.dist.to.pot / 1000)) +
        scale_color_discrete(aesthetics = c("color", "fill")) +
        labs(
        # title = "Rates of Change",
        subtitle = str_glue("{str_to_title(treatment)} and Control"),
        x = "Distance to Treatment (d) [km]", y = latex2exp::TeX(r"{Rate of Change \[pp/km\]}") ,
        caption = "Line: Median. Outer ribbon: 80% credible interval. Inner ribbon: 50% credible interval." , 
        colour = "Treatment (z)", fill = "Treatment (z)" 
        ) +
        theme(legend.position = "bottom") +
        NULL

    return(plot)
}

treatments = c(
    "bracelet",
    "calendar",
    "ink"
)


single_roc_plots = map(treatments, ~ roc_df %>%
    roc_plot(treatment = .x)
)
iwalk(
    single_roc_plots,
    ~ggsave(plot = .x, 
    filename = file.path(output_basepath, str_glue("single-roc-{treatments[.y]}.png")),
    width = 7.5, height = 5.0, dpi = 500 )
)



#### Difference in Rate of Change ####
diff_roc_df = dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
  transmute(
    model, model_name,
    y_rate_of_change_diff = map2(cluster_roc_diff, stan_data, ~ {
      mutate(.x,
        roc_distance = roc_distance / 1000,
        across(starts_with("per_"), divide_by, sd(.y$analysis_data$cluster.dist.to.pot)), 
        across(starts_with("per_"), multiply_by, 1000),
        across(starts_with("per_"), multiply_by, 100)
      ) 
  })) %>%
  unnest(y_rate_of_change_diff) %>% 
  filter(fct_match(assigned_treatment, c("ink", "bracelet", "calendar"))) 

plot_roc_diff = function(data, treatment) {
    plot = 
        data %>%
            filter(assigned_treatment == treatment) %>%
            mutate(assigned_treatment = fct_relabel(assigned_treatment, str_to_title)) %>%
            ggplot(aes(roc_distance)) +
            geom_line(aes(y = per_0.5)) +
            geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75), alpha = 0.25) +
            geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9), alpha = 0.25) +
            geom_rug(
              aes(dist, color = assigned_treatment),
              alpha = 0.75,
              data = rf_analysis_data %>%
                filter(fct_match(assigned.treatment, c("control", treatment))) %>%
                mutate(assigned.treatment = fct_relabel(assigned.treatment, str_to_title)) %>%
                distinct(cluster_id, assigned_treatment = assigned.treatment, dist = cluster.dist.to.pot / 1000)) +
            scale_color_discrete("", aesthetics = c("color", "fill")) +
            labs(
            # title = "Difference in Takeup Derivative",
            # subtitle = str_glue("{str_to_title(treatment)}"),
            x = "Distance to Treatment [km]", y = "Difference in Takeup Derivative [pp/km]",
            caption = "Line: Median. Outer ribbon: 80% credible interval. Inner ribbon: 50% credible interval." 
            ) +
            #   facet_wrap(vars(assigned_treatment), labeller = labeller(.default = str_to_title)) +
            guides( colour = "none") +
            theme_bw() +
            theme(legend.position = "bottom") +
            # coord_cartesian(ylim = c(0, 0.15)) + 
            NULL 
    return(plot)
}  




single_roc_diff_plots = map(treatments, ~ diff_roc_df %>%
    plot_roc_diff(treatment = .x)
)


iwalk(
    single_roc_diff_plots,
    ~ggsave(plot = .x, 
    filename = file.path(output_basepath, str_glue("single-roc-diff-{treatments[.y]}.png")),
    width = 7.5, height = 5.0, dpi = 500 )
)



ggsave(file.path(output_basepath, str_glue("dist_fit{fit_version}-diff-roc-by-treat.png")), width = 7.5, height = 5.0, dpi = 500)

dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
  transmute(
    model, model_name,
    y_rate_of_change_diff = map2(cluster_roc_diff, stan_data, ~ {
      mutate(.x,
        roc_distance = roc_distance / 1000,
        across(starts_with("per_"), divide_by, sd(.y$analysis_data$cluster.dist.to.pot)), 
        across(starts_with("per_"), multiply_by, 1000),
        across(starts_with("per_"), multiply_by, 100)
      ) 
  })) %>%
  unnest(y_rate_of_change_diff) %>% 
  filter(fct_match(assigned_treatment, c("ink", "bracelet"))) %>% 
  ggplot(aes(roc_distance)) +
  geom_line(aes(y = per_0.5, color = model_name)) +
  geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = model_name), alpha = 0.25) +
  geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = model_name), alpha = 0.25) +
  scale_color_discrete("", aesthetics = c("color", "fill")) +
  labs(
    title = "Difference in Takeup Derivative",
    x = "Distance to Treatment [km]", y = "Difference in Takeup Derivative [pp/km]",
    caption = "Line: Median. Outer ribbon: 80% credible interval. Inner ribbon: 50% credible interval." 
  ) +
  facet_wrap(vars(assigned_treatment), labeller = labeller(.default = str_to_title)) +
  guides(fill = "none", colour = "none") +
  theme_bw() +
  theme(legend.position = "bottom") +
  # coord_cartesian(ylim = c(0, 0.15)) + 
  NULL 
ggsave(file.path(output_basepath, str_glue("dist_fit{fit_version}-diff-roc-facet-treat.png")), width = 7.5, height = 5.0, dpi = 500)


#### Difference Rep Returns by Dist ####


dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
  select(model_name, cluster_rep_return_dist) %>%
  unnest(cluster_rep_return_dist) %>% 
  mutate(
    roc_distance = roc_distance / 1000,
    across(starts_with("per_"), divide_by, 1000)
  ) %>% 
  ggplot(aes(roc_distance)) +
  geom_line(aes(y = per_0.5, color = model_name)) +
  geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = model_name), alpha = 0.4) +
  geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = model_name), alpha = 0.4) +
  scale_color_discrete("", aesthetics = c("color", "fill")) +
  labs(
    title = "Valuation of Reputational Returns in Terms of Distance",
    x = "Distance to Treatment [km]", y = "Distance Value [km]",
  ) +
  facet_wrap(vars(assigned_treatment), labeller = labeller(.default = str_to_title), nrow = 1) +
  theme(legend.position = "top") +
  NULL

dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
  select(model_name, cluster_rep_return_dist) %>%
  unnest(cluster_rep_return_dist) %>% 
  mutate(
    roc_distance = roc_distance / 1000,
    across(starts_with("per_"), divide_by, 1000)
  ) %>% 
  ggplot(aes(roc_distance)) +
  geom_line(aes(y = per_0.5, color = assigned_treatment)) +
  geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = assigned_treatment), alpha = 0.4) +
  geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = assigned_treatment), alpha = 0.4) +
  scale_color_discrete("", aesthetics = c("color", "fill")) +
  labs(
    title = "Valuation of Reputational Returns in Terms of Distance",
    x = "Distance to Treatment [km]", y = "Distance Value [km]",
  ) +
#   facet_wrap(vars(assigned_treatment), labeller = labeller(.default = str_to_title), nrow = 1) +
  theme(legend.position = "bottom") +
  NULL
ggsave(file.path(output_basepath, str_glue("dist_fit{fit_version}-rep-returns-dist-by-treat.png")), width = 7.5, height = 5.0, dpi = 500)

dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
  select(model_name, cluster_rep_return_dist) %>%
  unnest(cluster_rep_return_dist) %>% 
  mutate(
    roc_distance = roc_distance / 1000,
    across(starts_with("per_"), divide_by, 1000)
  ) %>% 
  ggplot(aes(roc_distance)) +
  geom_line(aes(y = per_0.5, color = assigned_treatment)) +
  geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = assigned_treatment), alpha = 0.4) +
  geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = assigned_treatment), alpha = 0.4) +
  scale_color_discrete("", aesthetics = c("color", "fill")) +
  labs(
    title = "Valuation of Reputational Returns in Terms of Distance",
    x = "Distance to Treatment [km]", y = "Distance Value [km]",
  ) +
  facet_wrap(vars(assigned_treatment), labeller = labeller(.default = str_to_title), nrow = 2) +
  theme_bw() +
  theme(legend.position = "top") + 
  guides(fill = "none", colour = "none") +
  NULL
ggsave(file.path(output_basepath, str_glue("dist_fit{fit_version}-rep-returns-dist-facet-treat.png")), width = 7.5, height = 5.0, dpi = 500)


## Rep returns one by one

plot_rep_returns_one_by_one = function(data, treatment) {
    subtitle_str = str_to_title(treatment)
    plot = data %>%
        filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
        select(model_name, cluster_rep_return_dist) %>%
        unnest(cluster_rep_return_dist) %>% 
        filter(assigned_treatment == treatment) %>%
        mutate(
        roc_distance = roc_distance / 1000,
        across(starts_with("per_"), divide_by, 1000)
        ) %>% 
        ggplot(aes(roc_distance)) +
        geom_line(aes(y = per_0.5)) +
        geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75), alpha = 0.4) +
        geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9), alpha = 0.4) +
        geom_rug(
          aes(dist, color = assigned_treatment),
          alpha = 0.75,
          data = rf_analysis_data %>%
            filter(fct_match(assigned.treatment, c(treatment))) %>%
            distinct(cluster_id, assigned_treatment = assigned.treatment, dist = cluster.dist.to.pot / 1000)) +
        scale_color_discrete("", aesthetics = c("color", "fill")) +
        labs(
        # title = "Valuation of Reputational Returns in Terms of Distance",
        x = "Distance to Treatment [km]", y = "Distance Value [km]", 
        subtitle = subtitle_str
        ) +
        theme(legend.position = "top") +
        guides(colour = "none") +
        NULL
    return(plot)
}


treatments = c(
    "bracelet",
    "calendar",
    "ink",
    "control"
)

single_rep_return_plots = map(treatments, ~ dist_fit_data %>%
    plot_rep_returns_one_by_one(treatment = .x)  
)
iwalk(
    single_rep_return_plots,
    ~ggsave(plot = .x, 
    filename = file.path(output_basepath, str_glue("single-rep-return-{treatments[.y]}.png")),
    width = 7.5, height = 5.0, dpi = 500 )
)

## Rep returns vs control
if (model_fit_by == "Ed") {

  rep_return_df = read_csv(str_interp("temp-data/processed_rep_return_dist_fit${fit_version}.csv")) 
  rep_return_plot = function(data, treatment) {
      plot = data %>%
          mutate(
            roc_distance = roc_distance / 1000,
            across(starts_with("per_"), divide_by, 1000)
          ) %>% 
          ggplot(aes(roc_distance)) +
          geom_line(aes(y = per_0.5)) +
          geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75), alpha = 0.4) +
          geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9), alpha = 0.4) +
          geom_rug(
            aes(dist, color = assigned_treatment),
            alpha = 0.75,
            data = rf_analysis_data %>%
              filter(fct_match(assigned.treatment, c("control", treatment))) %>%
              distinct(cluster_id, assigned_treatment = assigned.treatment, dist = cluster.dist.to.pot / 1000)) +
          scale_color_discrete("", aesthetics = c("color", "fill")) +
          labs(
              # title = str_glue("Valuation of Reputational Returns in Terms of Distance"),
              subtitle =  str_glue("{str_to_title(treatment)} Compared to Control"),
              x = "Distance to Treatment [km]", y = "Distance Value Compared to Control[km]",
              caption = "Line: Median. Outer ribbon: 80% credible interval. Inner ribbon: 50% credible interval."
          ) +
          theme(legend.position = "top") +
          geom_hline(yintercept = 0, linetype = "longdash") + 
          guides(colour = "none") +
          NULL
      return(plot)
  }


  treatments = c(
      "bracelet",
      "calendar",
      "ink"
  )

  comp_rep_return_plots = map(treatments, ~ rep_return_df %>%
      filter(assigned_treatment == .x) %>%
      rep_return_plot(treatment = .x)
  )

  iwalk(
      comp_rep_return_plots,
      ~ggsave(plot = .x, 
      filename = file.path(output_basepath, str_glue("comp-rep-return-{treatments[.y]}.png")),
      width = 7.5, height = 5.0, dpi = 500 )
  )
}


#### Different Distributions and Delta ####

find_truncated_normal_mean = function(mu, sigma, a, b){
    alpha = (a - mu)/sigma
    beta = (b - mu)/sigma

    Z = pnorm(beta) - pnorm(alpha)

    tmean = mu + sigma*(dnorm(alpha) - dnorm(beta))/Z
    return(tmean)
}


a = find_truncated_normal_mean(0, 1, -Inf, 0.5)
b = find_truncated_normal_mean(0, 1, 0.5, Inf)



calculate_normal_delta = function(cutoff, mu, sigma) {
    m_minus = find_truncated_normal_mean(mu = mu, sigma = sigma, a = -Inf, b = cutoff )
    m_plus = find_truncated_normal_mean(mu = mu, sigma = sigma, a = cutoff, b = Inf )
    return(m_plus - m_minus)
}

kappa <- function(df, a, b) {
  gamma((df + 1)/2)/((pt(b, df = df) - pt(a, df = df))*gamma(df/2)*(df*pi)^(1/2))
}

tau <- function(df, j) {
  (df - 2*j)/df
}

ex <- function(df, a, b) {
  ((kappa(df = df, a = a, b = b)*df)/(df - 1))*((1 + a^2/df)^(-(df - 1)/2) - (1 + b^2/df)^(-(df - 1)/2))
}

ex2 <- function(df, a, b) {
  ((df - 1)/tau(df = df, j = 1))*((pt(b*sqrt(tau(df = df, j = 1)), df = (df - 2)) - pt(a*sqrt(tau(df = df, j = 1)), df = (df - 2)))/(pt(b, df = df) - pt(a, df = df))) - df
}

calculate_t_delta = function(df, cutoff) {
    m_minus = ex(df, a = -Inf, b = cutoff )
    m_plus = ex(df, a = cutoff, b = Inf )
    return(m_plus - m_minus)
}

## Conditional Log Normal Expectation
E_ln_minus = function(k, mu, sigma){
    exp(mu + (sigma^2)/2) * pnorm((log(k) - mu - sigma^2)/sigma) / pnorm((log(k) - mu)/sigma)
}

E_ln_plus = function(k, mu, sigma){
    exp(mu + (sigma^2)/2) * pnorm((-log(k) +  mu + sigma^2)/sigma) / (1 - pnorm((log(k) - mu)/sigma))
}

calculate_ln_delta = function(cutoff, mu, sigma){
    m_minus = E_ln_minus(mu = mu, sigma = sigma, k = cutoff )
    m_plus = E_ln_plus(mu = mu, sigma = sigma, k = cutoff )
    return(m_plus - m_minus)
}



ml = 1
sl = 1
df = expand.grid(
    w = seq(from = 0, to = 6, length.out = 100)
) %>%
    as_tibble() %>%
    mutate(
        delta = calculate_ln_delta(mu = ml, sigma = sl, cutoff = w)
    )
df %>%
    ggplot(aes(
        x = w, 
        y = delta
    )) +
    geom_point() +
    stat_function(fun = dlnorm, n = 200, args = list(meanlog = ml, sdlog = sl)) + ylab("")  +
    geom_vline(xintercept = exp(ml - sl^2), linetype = "longdash") + 
    labs(
        title = "Log Normal Prosociality Distribution", 
        subtitle = "Density's mode doesn't correspond to Delta(w)'s minimum" )
ggsave("temp-data/skew-delta-w.png", width = 8, height = 6, dpi = 500)


#### Delta[w] ####

delta_w_plot = dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS"))) %>% 
  select(model_name, sim_delta) %>% 
  unnest(sim_delta) %>% 
  ggplot(aes(w)) +
  geom_line(aes(y = per_0.5)) +
  geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75), alpha = 0.3) +
  geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9), alpha = 0.3) +
  labs(
    x = "w", y = "", 
    caption = "Line: Median. Outer ribbon: 80% credible interval. Inner ribbon: 50% credible interval."
  ) +
  # coord_cartesian(ylim = c(0, 0.1))
  NULL

delta_w_plot

ggsave(
    plot= delta_w_plot, 
    filename = file.path(
        output_basepath,
        str_glue("delta-w-plot.png")
    ),
    width = 7.5, 
    height = 5.0,
    dpi = 500
)



## ATEs
grab_incentive_ate = function(data, model_type_to_plot, nested_data = est_takeup_te) {

  subset_data = data %>%
    filter(
      (fct_match(model_type, model_type_to_plot))
    ) %>% 
    mutate(
      {{ nested_data }} :=
        map_if({{ nested_data }}, fct_match(model_type, "structural"),
              filter, mu_assigned_treatment_left == assigned_treatment_left, mu_assigned_treatment_right == assigned_treatment_right) %>%
          map(filter,
              (is.na(assigned_dist_group_left) & is.na(assigned_dist_group_right)) | (assigned_dist_group_left == assigned_dist_group_right),
              assigned_treatment_left != assigned_treatment_right,
              fct_match(assigned_treatment_right, c("control")),
              fct_match(assigned_treatment_left, "bracelet") | !fct_match(assigned_treatment_right, "calendar"))
    )  %>%
    select(model, model_name, {{ nested_data }}, fit_type, model_color) %>% 
    mutate(
      {{ nested_data }} := map(
        {{ nested_data }},
        mutate,
        assigned_dist_group_left = fct_explicit_na(assigned_dist_group_left, "Combined") %>% 
          fct_relabel(str_to_title) %>% 
          fct_relevel("Combined"),
        assigned_treatment_left = fct_rev(factor(str_to_title(assigned_treatment_left)))
      )
    ) 
    return(subset_data)
}

reduced_incentive_level_plot = dist_fit_data %>%
    filter(
      (fct_match(model_type, "reduced form"))
    )  %>%
    mutate(est_takeup_level = map_if(est_takeup_level,
        fct_match(fit_type, "prior-predict"), 
        ~mutate(.x, 
          across(
            c(contains("per"), mean_est),
            ~if_else(
              !(fct_match(assigned_treatment,  "bracelet") & is.na(assigned_dist_group)),
              NA_real_,
              .x 
            )
          )
        )
      )
    ) %>%
    mutate(
      est_takeup_level = map(
        est_takeup_level,
        mutate,
        assigned_dist_group = fct_explicit_na(assigned_dist_group, "Combined") %>% 
          fct_relabel(str_to_title) %>% 
          fct_relevel("Combined"),
        assigned_treatment = fct_rev(factor(str_to_title(assigned_treatment), 
                                    levels = c("Bracelet",
                                                "Calendar",
                                                "Ink",
                                                "Control")))
        )
        ) %>%
  plot_estimands(.,est_takeup_level, assigned_treatment, results_group = fit_type) +
    scale_x_continuous("", breaks = seq(-1, 1, 0.1)) +
    scale_y_discrete("") +
    labs(
      # title = "Incentive Average Treatment Effect"
    ) +
    ggforce::facet_col(vars(assigned_dist_group), 
                space = "free",
                scales = "free_y") +
    NULL  +
    guides(colour = "none")
reduced_incentive_level_plot
ggsave(
  plot = reduced_incentive_level_plot,
  filename = file.path(
    output_basepath,
    "reduced-incentive-level-plot.png"
  ),
  width = 7.5,
  height = 5,
  dpi = 500
)


reduced_incentive_ate_plot = dist_fit_data %>%
  grab_incentive_ate(., model_type_to_plot = "reduced form") %>%
  plot_estimands(., 
                 est_takeup_te, 
                 assigned_treatment_left, 
                 results_group = fit_type, 
                 single_prior_predict = TRUE, 
                 top_levels = default_top_levels) +
    scale_x_continuous("", breaks = seq(-0.2, 0.2, 0.05)) +
    scale_y_discrete("") +
    labs(
      # title = "Incentive Average Treatment Effect"
    ) +
    ggforce::facet_col(vars(assigned_dist_group_left), 
                space = "free",
                scales = "free_y") +
    NULL  +
    guides(colour = "none")
reduced_incentive_ate_plot
ggsave(
  plot = reduced_incentive_ate_plot,
  filename = file.path(
    output_basepath,
    "reduced-incentive-ate-plot.png"
  ),
  width = 7.5,
  height = 5,
  dpi = 500
)

default_ate_limits = c(-0.1, 0.2)

structural_incentive_ate_plot = dist_fit_data %>%
  grab_incentive_ate(., model_type_to_plot = "structural") %>%
  plot_estimands(., 
                 est_takeup_te, 
                 assigned_treatment_left, 
                 results_group = fit_type, 
                 single_prior_predict = TRUE, 
                 top_levels = default_top_levels) +
    scale_x_continuous(
      "", 
      breaks = seq(-0.2, 0.2, 0.05), 
      limits = default_ate_limits) +
    scale_y_discrete("") +
    labs(
      # title = "Incentive Average Treatment Effect"
    ) +
    ggforce::facet_col(vars(assigned_dist_group_left), 
                space = "free",
                scales = "free_y") +
    NULL  +
    guides(colour = "none") 

structural_incentive_ate_plot

ggsave(
  plot = structural_incentive_ate_plot,
  filename = file.path(
    output_basepath,
    "structural-incentive-ate-plot.png"
  ),
  width = 7.5,
  height = 5.0,
  dpi = 500
)

## SIGNALLING

grab_signalling_ate = function(data, model_type_to_plot, models_we_want) {
  data = data %>% 
    filter(
      (fct_match(model_type, model_type_to_plot) & fct_match(model, models_we_want) )  
    ) %>% 
    select(model, model_name, est_takeup_te, fit_type, model_color) %>% 
    mutate(
      est_takeup_te = map(
        est_takeup_te,
        filter,
        (is.na(assigned_dist_group_left) & is.na(assigned_dist_group_right)) | (assigned_dist_group_left == assigned_dist_group_right),
        across(c(assigned_treatment_left, assigned_treatment_right), fct_match, "control"),
        !is.na(mu_assigned_treatment_left),
        fct_match(mu_assigned_treatment_left, "bracelet") | !fct_match(mu_assigned_treatment_right, "calendar"),
        fct_match(mu_assigned_treatment_right, "control"),
      ) %>% 
        map(
          mutate,
          assigned_dist_group_left = fct_explicit_na(assigned_dist_group_left, "Combined") %>% 
            fct_relabel(str_to_title) %>% 
            fct_relevel("Combined"),
          mu_assigned_treatment_left = fct_rev(factor(str_to_title(mu_assigned_treatment_left))),
        )
    ) 
  return(data)

}


structural_signalling_ate_plot = dist_fit_data %>%
  grab_signalling_ate(., "structural", models_we_want) %>%
      plot_estimands(., 
                     est_takeup_te, 
                     mu_assigned_treatment_left, 
                     results_group = fit_type,
                     single_prior_predict = TRUE, 
                     top_levels = default_top_levels) +
        scale_x_continuous("", breaks = seq(-0.2, 0.2, 0.05), 
          limits = default_ate_limits
        ) +
        scale_y_discrete("") +
        labs(
          # title = "Signaling Average Treatment Effect",
          # subtitle = str_glue("Holding private incentive at the control level.")
          ) +
        ggforce::facet_col(vars(assigned_dist_group_left), 
                   space = "free",
                   scales = "free_y") +
        guides(colour = "none") +
        NULL
structural_signalling_ate_plot
ggsave(
  plot = structural_signalling_ate_plot,
  filename = file.path(
    output_basepath,
    "structural-signaling-ate-plot.png"
  ),
  width = 7.5,
  height = 5.0,
  dpi = 500
)



grab_private_ate = function(data, model_type_to_plot, models_we_want) {
  data = data %>%  
    filter(
      (fct_match(model_type, model_type_to_plot) & fct_match(model, models_we_want) )  
    ) %>% 
    select(model, model_name, est_takeup_te, fit_type, model_color) %>% 
    mutate(
      est_takeup_te = map(
        est_takeup_te,
        filter,
        (is.na(assigned_dist_group_left) & is.na(assigned_dist_group_right)) | (assigned_dist_group_left == assigned_dist_group_right),
        !is.na(mu_assigned_treatment_left),
        !is.na(mu_assigned_treatment_right),
        across(c(mu_assigned_treatment_left, mu_assigned_treatment_right), fct_match, "control"),
        fct_match(assigned_treatment_left, "bracelet") | !fct_match(assigned_treatment_right, "calendar"),
        fct_match(assigned_treatment_right, "control"),
      ) %>% 
        map(
          mutate,
          assigned_dist_group_left = fct_explicit_na(assigned_dist_group_left, "Combined") %>% 
            fct_relabel(str_to_title) %>% 
            fct_relevel("Combined"),
          assigned_treatment_left = fct_rev(factor(str_to_title(assigned_treatment_left))),
        )
    )
  return(data)
}

structural_private_ate_plot = dist_fit_data %>%
  grab_private_ate(., model_type_to_plot = "structural", models_we_want) %>%
      plot_estimands(., 
                     est_takeup_te, 
                     assigned_treatment_left, 
                     results_group = fit_type, 
                     single_prior_predict = TRUE, 
                     top_levels = default_top_levels) +
        scale_x_continuous("", breaks = seq(-0.2, 0.2, 0.05), 
          limits =  default_ate_limits
          ) +
        scale_y_discrete("") +
        labs(
          # title = "Private Incentive Average Treatment Effect",
          # subtitle = str_glue("Holding signaling at the control level")
          ) +
        ggforce::facet_col(vars(assigned_dist_group_left), space = "free", scales = "free_y") +
      guides(colour = "none") +
        NULL
structural_private_ate_plot
ggsave(
  plot = structural_private_ate_plot,
  filename = file.path(
    output_basepath,
    "structural-private-ate-plot.png"
  ),
  width = 7.5,
  height = 5.0,
  dpi = 500
)



## Differences Between Reduced Form Treatment Arms

# TODO: remove est_takeup_te_diff - don't actually need can just filter here


#' Calculate Difference in ATEs across treatments, within distance groups
#'
#' Filtering process:
#' Ensure treat_left != treat_right. Ensure mu_assigned either matches treatment 
#' or is set to control in RF. Make sure distance groups match or are both NA
#'
#'
grab_incentive_ate_diff = function(data, base_comparison, model_type_to_plot, models_we_want){
  data =  data %>%
    filter(
      (fct_match(model_type, model_type_to_plot) & fct_match(model, models_we_want) )  
    ) %>% 
    select(model, model_name, est_takeup_te, fit_type, model_type, model_color) %>% 
      mutate(
        est_takeup_te = map_if(est_takeup_te,
        fct_match(model_type, "structural"), ~ {
          if (!is_null(.x)) {
            filter(.,
              (assigned_treatment_left != assigned_treatment_right) &
              mu_assigned_treatment_left == assigned_treatment_left &
              mu_assigned_treatment_right == assigned_treatment_right &
              ((is.na(assigned_dist_group_left) & is.na(assigned_dist_group_right)) |
              (assigned_dist_group_left == assigned_dist_group_right) &
             fct_match(assigned_treatment_right, base_comparison))
            )
          }
        },
        .else = ~ filter(.,
          (assigned_treatment_left != assigned_treatment_right) &
          mu_assigned_treatment_left == "control" &
          mu_assigned_treatment_right == "control" &
          ((is.na(assigned_dist_group_left) & is.na(assigned_dist_group_right)) |
          (assigned_dist_group_left == assigned_dist_group_right)) &
          fct_match(assigned_treatment_right,  base_comparison))
      ) %>%
        map(
          mutate,
          assigned_dist_group_left = fct_explicit_na(assigned_dist_group_left, "Combined") %>% 
            fct_relabel(str_to_title) %>% 
            fct_relevel("Combined"),
          assigned_treatment_left = fct_rev(factor(str_to_title(assigned_treatment_left)))
        )
      )
  return(data)
}


plot_ate_diff = function(data, base_comparison, model_type_to_plot, models_we_want) {
  if (base_comparison == "bracelet") {
    top_levels = c("Calendar", "Combined")
  } else {
    top_levels = c("Bracelet", "Combined")
  }
  p = data %>%
    grab_incentive_ate_diff(data = .,
                            base_comparison = base_comparison, 
                            model_type_to_plot = model_type_to_plot, 
                            models_we_want =  models_we_want
                            ) %>%
    plot_estimands(., 
                   est_takeup_te, 
                   assigned_treatment_left, 
                   results_group = fit_type,
                   single_prior_predict = TRUE,
                   top_levels = top_levels) +
      scale_x_continuous("", breaks = seq(-1, 1, 0.05)) +
      scale_y_discrete("") +
      labs(
        title = str_glue("Incentive Average Treatment Effect Difference"), 
        subtitle =  str_glue("Treatment Effect Difference Compared to {str_to_title(base_comparison)}")) +
      ggforce::facet_col(vars(assigned_dist_group_left), space = "free", scales = "free_y") +
    guides(colour = "none") +
      NULL
  return(p)
}





p_ate_diffs = map(
  treatments, 
  ~plot_ate_diff(
    data = dist_fit_data,
    base_comparison = .x, 
    model_type_to_plot = "reduced form",
    models_we_want = "REDUCED_FORM_NO_RESTRICT"
  )
)

iwalk(
    p_ate_diffs,
    ~ggsave(plot = .x, 
    filename = file.path(output_basepath, str_glue("reduced-form-incentive-ate-diffs-{treatments[.y]}.png")),
    width = 7.5, height = 5.0, dpi = 500 )
)



## Social Multiplier ##


# dist_fit_data %>%
#   filter(
#     fct_match(model_type, "structural")
#     ) %>%
#     select(fit_type, cluster_social_multiplier) %>%
#     unnest() %>%
#     filter(assigned_treatment == "control") %>%
#     select( 
#       fit_type, roc_distance, mean_est
#     ) %>%
#     spread(fit_type, mean_est) %>%
#     mutate(diff = fit - `prior-predict`) %>%
#     ggplot(aes(
#       x = roc_distance, 
#       y = diff
#     )) +
#     geom_point()


# dist_fit_data %>%
#   filter(
#     fct_match(model_type, "structural")
#     ) %>%
#     select(fit_type, cluster_social_multiplier) %>%
#     filter(fit_type == "fit") %>%
#     unnest() %>%
#     filter(assigned_treatment == "control") %>%
#     mutate(roc_distance = roc_distance / 1000) %>%
#     ggplot(aes(roc_distance)) +
#     geom_line(aes( 
#       y = mean_est,
#       group = fit_type, 
#       linetype = fit_type
#     ),  colour = "black") 

# dist_fit_data %>%
#   filter(
#     fct_match(model_type, "structural")
#     ) %>%
#     select(fit_type, cluster_social_multiplier) %>%
#     unnest() %>%
#     filter(assigned_treatment == "control") %>%
#     mutate(roc_distance = roc_distance / 1000) %>%
#     ggplot(aes(roc_distance)) +
#     geom_line(aes( 
#       y = mean_est,
#       group = fit_type, 
#       linetype = fit_type
#     ),  colour = "black")  +
#     geom_ribbon(
#       data = . %>% filter(fit_type == "prior-predict"),
#       aes(
#         ymin = per_0.25,
#         y = mean_est,
#         ymax = per_0.75, 
#         fill = fit_type
#       ), 
#       alpha = 0.2
#     ) +
#     geom_ribbon(
#       data = . %>% filter(fit_type != "prior-predict"),
#       aes(
#         ymin = per_0.25,
#         y = mean_est,
#         ymax = per_0.75, 
#         fill = fit_type
#       ), 
#       alpha = 0.4
#     ) +
#     geom_hline(yintercept = -1, linetype = "longdash") +
#     geom_rug(
#       aes(dist),
#       alpha = 0.75,
#       data = rf_analysis_data %>%
#         mutate(assigned.treatment = fct_relabel(assigned.treatment, str_to_title)) %>%
#         distinct(cluster_id, assigned_treatment = assigned.treatment, dist = cluster.dist.to.pot / 1000)) +
#     labs(
#       x = "Distance to Treatment (d) [km]", 
#       y = "Social Multiplier" , 
#       fill = "", 
#       colour = "", 
#       linetype = "",
#       caption = "Line: Median. Ribbon: 50% credible interval."
#     ) +
#     scale_fill_discrete(
#       breaks = c("fit", "prior-predict"),
#       labels = c("Structural Model", "Prior")
#     ) +
#     scale_linetype_discrete(
#       breaks = c("fit", "prior-predict"),
#       labels = c("Structural Model", "Prior")
#     ) +
#     xlim(0, 2.5) +
#     NULL


# ggsave(
#   file.path(
#     output_basepath, 
#     str_glue("social-multiplier-estim-plot.png")), 
#     width = 7.5, height = 5.0, dpi = 500)


# dist_fit_data %>%
#   filter(
#     fct_match(model_type, "structural")
#     ) %>%
#     select(fit_type, cluster_social_multiplier) %>%
#     filter(fit_type == "fit") %>%
#     unnest() %>%
#     # filter(assigned_treatment == "control") %>%
#     mutate(roc_distance = roc_distance / 1000) %>%
#     ggplot(aes(roc_distance)) +
#     geom_line(aes( 
#       y = mean_est,
#       colour = assigned_treatment
#     ))  +
#     # geom_hline(yintercept = -1, linetype = "longdash")
#     # geom_ribbon(
#     #   data = . %>% filter(fit_type != "prior-predict"),
#     #   aes(
#     #     ymin = per_0.25,
#     #     y = mean_est,
#     #     ymax = per_0.75, 
#     #     fill = fit_type
#     #   ), 
#     #   alpha = 0.4
#     # ) +
#     geom_rug(
#       aes(dist),
#       alpha = 0.75,
#       data = rf_analysis_data %>%
#         mutate(assigned.treatment = fct_relabel(assigned.treatment, str_to_title)) %>%
#         distinct(cluster_id, assigned_treatment = assigned.treatment, dist = cluster.dist.to.pot / 1000)) +
#     labs(
#       x = "Distance to Treatment (d) [km]", 
#       y = "Social Multiplier" , 
#       fill = "", 
#       colour = "", 
#       linetype = "",
#       caption = "Line: Median"
#     )  +
#     xlim(0, 2.5) +
#     facet_wrap(~assigned_treatment)

# ggsave(
#   file.path(
#     output_basepath, 
#     str_glue("social-multiplier-point-estim-plot.png")), 
#     width = 7.5, height = 5.0, dpi = 500)



# dist_fit_data %>%
#   filter(
#     fct_match(model_type, "structural")
#     ) %>%
#     select(fit_type, cluster_social_multiplier) %>%
#     unnest() %>%
#     filter(assigned_treatment == "control") %>%
#     ggplot(aes(roc_distance)) +
#     geom_line(aes(y = per_0.5, color = fit_type)) +
#     geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = fit_type), alpha = 0.4) +
#     geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = fit_type), alpha = 0.4) +
#     # geom_rug(
#     #   aes(dist, color = assigned_treatment),
#     #   alpha = 0.75,
#     #   data = rf_analysis_data %>%
#     #     mutate(assigned.treatment = fct_relabel(assigned.treatment, str_to_title)) %>%
#     #     distinct(cluster_id, assigned_treatment = assigned.treatment, dist = cluster.dist.to.pot / 1000)) +
#     scale_color_discrete(aesthetics = c("color", "fill"))  +




