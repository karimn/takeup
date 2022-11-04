
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
    ) 

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
    ) 

combined_belief_prop %>%
  plot_single_beliefs_est(
    width = 0.7, 
    order = 1,
    crossbar_width = 0.4, 
    vline = FALSE) +
    labs(subtitle = "Proportion") +
    NULL



combined_belief_ate %>%
    plot_single_beliefs_est(order = 1)  +
    labs(subtitle = "Treatment Effects")
ggsave(
    file.path(output_basepath, str_glue("combined-ate-fob.png")), 
    width = 7.5, 
    height = 5.0, 
    dpi = 500
    )

combined_belief_ate %>%
    plot_single_beliefs_est(order = 2)  +
    labs(subtitle = "Treatment Effects")
ggsave(
    file.path(output_basepath, str_glue("combined-ate-sob.png")), 
    width = 7.5, 
    height = 5.0, 
    dpi = 500
    )


#### Rate of Change
roc_df = dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST", "STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
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
        ggplot(aes(roc_distance)) +
        geom_line(aes(y = per_0.5, color = assigned_treatment)) +
        geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = assigned_treatment), alpha = 0.4) +
        geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = assigned_treatment), alpha = 0.4) +
        scale_color_discrete("", aesthetics = c("color", "fill")) +
        labs(
        title = "Rates of Change",
        subtitle = str_glue("{str_to_title(treatment)} and Control"),
        x = "Distance to Treatment (d) [km]", y = latex2exp::TeX(r"{Rate of Change \[pp/km\]}") ,
        caption = "Line: Median. Outer ribbon: 80% credible interval. Inner ribbon: 50% credible interval." 
        ) +
        theme(legend.position = "top") +
        NULL

    return(plot)
}

treatments = c(
    "bracelet",
    "calendar",
    "ink"
)

single_roc_plots[[2]]

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
dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST", "STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
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
  geom_line(aes(y = per_0.5, color = assigned_treatment)) +
  geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75, fill = assigned_treatment), alpha = 0.25) +
  geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9, fill = assigned_treatment), alpha = 0.25) +
  scale_color_discrete("", aesthetics = c("color", "fill")) +
  labs(
    title = "Difference in Rate of Change",
    x = "Distance to Treatment [km]", y = "Difference in Rate of Change [pp/km]",
    caption = "Line: Median. Outer ribbon: 80% credible interval. Inner ribbon: 50% credible interval." 
  ) +
#   facet_wrap(vars(assigned_treatment), labeller = labeller(.default = str_to_title)) +
  guides( colour = "none") +
  theme_bw() +
  theme(legend.position = "bottom") +
  # coord_cartesian(ylim = c(0, 0.15)) + 
  NULL 
ggsave(file.path(output_basepath, str_glue("dist_fit{fit_version}-diff-roc-by-treat.png")), width = 7.5, height = 5.0, dpi = 500)

dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST", "STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
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
    title = "Difference in Rate of Change",
    x = "Distance to Treatment [km]", y = "Difference in Rate of Change [pp/km]",
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
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS", "STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST"))) %>%
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
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS", "STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST"))) %>%
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
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS", "STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST"))) %>%
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
        filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS", "STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST"))) %>%
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
        scale_color_discrete("", aesthetics = c("color", "fill")) +
        labs(
        title = "Valuation of Reputational Returns in Terms of Distance",
        x = "Distance to Treatment [km]", y = "Distance Value [km]", 
        subtitle = subtitle_str
        ) +
        theme(legend.position = "top") +
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

rep_return_df = read_csv(str_interp("temp-data/processed_rep_return_dist_fit${fit_version}.csv")) %>%
  mutate(
    roc_distance = roc_distance / 1000
  ) 

rep_return_plot = function(data, treatment) {
    plot = data %>%
        ggplot(aes(roc_distance)) +
        geom_line(aes(y = per_0.5)) +
        geom_ribbon(aes(ymin = per_0.25, ymax = per_0.75), alpha = 0.4) +
        geom_ribbon(aes(ymin = per_0.1, ymax = per_0.9), alpha = 0.4) +
        scale_color_discrete("", aesthetics = c("color", "fill")) +
        labs(
            title = str_glue("Valuation of Reputational Returns in Terms of Distance"),
            subtitle =  str_glue("{str_to_title(treatment)} Compared to Control"),
            x = "Distance to Treatment [km]", y = "Distance Value Compared to Control[km]",
        ) +
        theme(legend.position = "top") +
        geom_hline(yintercept = 0, linetype = "longdash") + 
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


## Delta[w]

delta_w_plot = dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS", "STRUCTURAL_LINEAR_U_SHOCKS_NO_BELIEFS_DIST"))) %>% 
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
