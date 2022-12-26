
library(tidyverse)
library(posterior)
library(cmdstanr)
library(tidybayes)

files = fs::dir_ls("data/stan_analysis_data", regexp = "72_REDUCED")
files = files[!str_detect(files, "rds")]
fit = as_cmdstan_fit(files)



sms_beta_draws = fit %>%
    gather_draws(
        reduced_beta_phone_group[phone_group_index, treatment_index]
    )

summ_sms_beta_draws = sms_beta_draws %>%
    median_qi() %>%
    to_broom_names()

summ_sms_beta_draws %>%
    ggplot(aes(
        x = estimate,
        xmin = conf.low, 
        xmax = conf.high, 
        colour = factor(treatment_index), 
        y = factor(treatment_index)
    )) +
    facet_wrap(~phone_group_index) +
    geom_pointrange() +
    geom_vline(xintercept = 0, linetype = "longdash")

dist_fit_data %>%
  select(fit) %>%
  unnest(fit) %>%
  filter(str_detect(variable, "reduced_beta_phone_group")) %>%
  filter(!str_detect(variable, "raw|alpha")) %>%
  unnest(iter_data)  %>%
  mutate(
    phone_group_index = str_extract(variable, "(?<=\\[)\\d(?=,)") %>% as.factor(), 
    treatment_index = str_extract(variable, "(?<=,)\\d(?=\\])") %>% as.factor()
  ) %>%
  ggplot(aes(
    x = iter_est, 
    fill = treatment_index
  )) +
  geom_histogram() +
  facet_wrap(~phone_group_index)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

# Data --------------------------------------------------------------------

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

monitored_sms_data <- analysis.data %>% 
  filter(mon_status == "monitored") %>% 
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



analysis_data <- monitored_sms_data %>% 
mutate(
assigned_treatment = assigned.treatment, 
assigned_dist_group = dist.pot.group, 
sms_treatment = sms.treatment.2, 
phone_owner = if_else(phone_owner == TRUE, "phone", "nophone"), 
sms_treatment = str_replace_all(sms_treatment, "\\.", "")) %>%
# reminder.only only present in control condition
filter(
    sms_treatment != "reminderonly"
)  %>%
filter(phone_owner == "phone") %>%
mutate(sms_treatment = factor(sms_treatment))

library(fixest)


base_fe_fit = analysis_data %>%
    feglm(
        dewormed ~ assigned_treatment*assigned_dist_group, 
        data = ., 
        family = binomial(link = "probit")
    )



sms_fe_fit = analysis_data %>%
    feglm(
        dewormed ~ assigned_treatment*assigned_dist_group*sms_treatment, 
        data = ., 
        family = binomial(link = "probit")
    )

install.packages("marginaleffects")
library(marginaleffects)

plot_cap(
    base_fe_fit,
    c("assigned_treatment", "assigned_dist_group")
)

plot_cap(
    sms_fe_fit,
    c(
    "assigned_treatment", 
    "assigned_dist_group", 
    "sms_treatment"
    )
)

pred_df = predictions(
    sms_fe_fit,
    newdata = datagrid(
        assigned_treatment = unique(analysis_data$assigned_treatment), 
        assigned_dist_group = unique(analysis_data$assigned_dist_group), 
        sms_treatment = unique(analysis_data$sms_treatment)
    )
) %>%
    as_tibble()

summ_sms_comp = comparisons(sms_fe_fit) %>%
    summary() 

summ_sms_comp



breakdown_sms_comp = comparisons(
    sms_fe_fit, 
    variables = "sms_treatment",
    newdata = datagrid(
        assigned_treatment = unique(analysis_data$assigned_treatment), 
        assigned_dist_group = unique(analysis_data$assigned_dist_group))
)



pred_df %>%
    ggplot(aes(
        x = predicted, 
        xmin = conf.low, 
        xmax = conf.high, 
        y = assigned_treatment,
        colour = sms_treatment
    )) +
    geom_pointrange(position = position_dodge2(0.4))  +
    facet_wrap(~assigned_dist_group, ncol = 1) +
    theme_bw() +
    theme(legend.position = "bottom")  +
    labs(
        title = "Estimated Takeup Level by SMS Treatment", 
        x = "Takeup", 
        y = ""
    )

ggsave(
    "temp-plots/takeup-level-by-sms.png",
    width = 8,
    height = 6,
    dpi = 500
)

breakdown_sms_comp %>%
    ggplot(aes(
        x = comparison, 
        xmin = conf.low, 
        xmax = conf.high, 
        y = assigned_treatment,
        colour =  assigned_dist_group
    )) + 
    geom_pointrange(position = position_dodge2(0.4)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(
        title = "Estimated SMS Takeup TE, by Close vs Far", 
        x = "Estimated TE"
    )

ggsave(
    "temp-plots/takeup-TE-by-sms.png",
    width = 8,
    height = 6,
    dpi = 500
)
