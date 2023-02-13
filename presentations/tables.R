#!/usr/bin/Rscript
script_options <- docopt::docopt(
  stringr::str_glue(
"Usage:
  tables.R <fit-version> [options] [<tables>...]
  
Options:
  --cores=<num-cores>  Number of cores to use [default: 12]
  --input-path=<input-path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
  --output-path=<output-path>  Path to find results [default: temp-data]
  --models=<models>
  "), 
  args = if (interactive()) "
    71  \
    --cores=1 \
    --models=STRUCTURAL_LINEAR_U_SHOCKS
    balance_indiv
    " else commandArgs(trailingOnly = TRUE)
)


library(tidyverse)
library(broom)
library(knitr)
library(kableExtra)
library(ggthemes)
library(fixest)
library(margins)

source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))
source(file.path("multilvlr", "multilvlr_util.R"))

fit_version = script_options$fit_version

treat_levels_c = c("control", "ink", "calendar", "bracelet")
treat_levels = c("ink", "calendar", "bracelet")
dist_levels = c("close", "far")

quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

output_basepath = file.path(
  script_options$output_path,
  str_glue("output_dist_fit{script_options$fit_version}")
)


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



if ("balance_indiv" %in% script_options$tables) {


  balance_variables = c(
    "age",
    "ethnicity",
    "phone_owner", 
    "Pupil.Teacher.Ratio",
    "Pupil.Classroom.Ratio",
    "Pupil.Toilet.Ratio",
    "Total.Number.of.Classrooms",
    "Total.Toilets",
    "Total.Boys",
    "Total.Girls",
    "Total.Enrolment", 
    "GOK.TSC.Male",
    "GOK.TSC.Female",
    "Non.Teaching.Staff.Male",
    "Non.Teaching.Staff.Female",
    "Local.Authority.Male",
    "Local.Authority.Female",
    "phone",
    "gender"
  )

  rct_school_df = rct.schools.data %>% 
    as_tibble()

  rct_school_df %>%
    colnames()
  analysis_data %>%
    colnames()

  analysis_school_data = left_join(
    analysis_data,
    rct_school_df %>% mutate(cluster.id = as.numeric(cluster.id)) ,
    by = "cluster.id"
  )







  analysis_school_data = analysis_school_data %>%
    mutate(
      treat_dist = paste0(
        "treat: ", 
        assigned.treatment,
        ", dist: ", dist.pot.group
        ) %>% factor()
  )  %>%
  mutate(
    female = fct_match(gender, "female")
  )

  indiv_balance_vars = c(
    "female",
    "age",
    "phone_owner"
  )

  analysis_data %>%
    select(phone_owner) %>%
    unique()

  indiv_balance_fit = analysis_school_data %>%
    feols(
      data = ., 
      .[indiv_balance_vars] ~ 0 + treat_dist
        ) 


  col_order = c(
    "lhs", 
    paste0(treat_levels_c, "_close"),
    paste0(treat_levels_c, "_far")
  )

  indiv_balance_tidy_df = indiv_balance_fit %>%
    map_dfr(tidy, .id = "lhs") %>%
    mutate(
      lhs = str_remove(lhs, "lhs: ")
    ) %>%
    select(
      lhs, term, estimate, std.error
    )  

  indiv_balance_tidy_df

  create_balance_input = function(tidy_df, col_order, digits = 3) {
    bal_input_df = tidy_df %>%
      mutate(
        estimate = round(estimate, digits),
        std.error = round(std.error, digits)
        ) %>%
      mutate(
        estim_std = linebreak(paste0(estimate, "\n", str_glue("({std.error})"))) 
      )  %>%
      mutate(lhs_treat = str_extract(term, "(?<=disttreat: ).*(?=\\,)")) %>%
      mutate(lhs_dist = str_extract(term, "(?<=dist: ).*$")) %>%
      select(lhs, lhs_treat, lhs_dist, estim_std) %>%
      mutate(
        lhs_treat = factor(lhs_treat, treat_levels_c) %>% fct_rev(), 
        lhs_dist = factor(lhs_dist, dist_levels) %>% fct_rev()
        ) %>%
      pivot_wider(
        names_from = c(lhs_treat, lhs_dist),
        values_from = estim_std, 
      )    %>%
      select(
        all_of(col_order)
        )  %>%
        mutate(
          lhs = str_to_title(lhs),
          lhs = str_replace_all(lhs, "_", " ")
        )
    return(bal_input_df)
  }


  wide_indiv_bal_df = create_balance_input(
    indiv_balance_tidy_df, 
    col_order = col_order, 
    digits = 3
  ) 

  wide_indiv_bal_df %>% 
    knitr::kable(
      format = "latex",
      escape = FALSE, 
      col.names = c("", rep(str_to_title(treat_levels_c), 2)), 
      booktabs = TRUE
    ) %>%
    kable_styling() %>%
    add_header_above(c("", "Close" = 4, "Far" = 4 ))



  library(gt)
  wide_df %>%
    gt() %>%
    as_kable_extra()

}



## Fit Loading

load(file.path("temp-data", str_interp("processed_dist_fit${fit_version}_lite.RData")))

dist_fit_data %<>%
  left_join(tribble(
    ~ fit_type,        ~ model_color,
      "fit",           "black", 
      "prior-predict", "darkgrey",
  ), by = "fit_type") %>%
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS", "REDUCED_FORM_NO_RESTRICT")))

delta <- function(v, ...) dnorm(v, ...) / ((pnorm(v, ...) * pnorm(v, ..., lower.tail = FALSE)))

belief_data = dist_fit_data %>%
  filter(fct_match(fit_type, "fit"), fct_match(model, c("STRUCTURAL_LINEAR_U_SHOCKS"))) %>%
  pull(beliefs_results) %>%
  first() 


cols_we_want = c(
    "est_takeup_level",
    "beliefs_results",
    "wtp_results",
    "est_takeup_te",
    "est_takeup_dist_te"
)


dist_fit_data %>%
    filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) %>%
    select(est_takeup_te) %>%
    unnest() %>%
    filter(
        mu_assigned_treatment_left == assigned_treatment_left,
        mu_assigned_treatment_right == assigned_treatment_right, 
        assigned_treatment_right == "control", 
        is.na(assigned_dist_group_right))

dist_fit_data %>%
    filter(fit_type == "fit") %>%
    select(any_of(cols_we_want)) %>%
    unnest(cols = c(est_takeup_level))





incentive_te = dist_fit_data %>%
  mutate(
    est_takeup_te =
      map_if(est_takeup_te, fct_match(model_type, "structural"),
             filter, mu_assigned_treatment_left == assigned_treatment_left, mu_assigned_treatment_right == assigned_treatment_right) %>%
        map(filter,
            (is.na(assigned_dist_group_left) & is.na(assigned_dist_group_right)) | (assigned_dist_group_left == assigned_dist_group_right),
            assigned_treatment_left != assigned_treatment_right,
            fct_match(assigned_treatment_right, c("control")),
            fct_match(assigned_treatment_left, "bracelet") | !fct_match(assigned_treatment_right, "calendar")),
    model_color = canva_pal(canva_palette_vibrant)(n())
  ) %>% 
  select(model, model_type, model_name, est_takeup_te, fit_type, model_color) %>% 
  mutate(
    est_takeup_te = map(
      est_takeup_te,
      mutate,
      assigned_dist_group_left = fct_explicit_na(assigned_dist_group_left, "Combined") %>% 
        fct_relabel(str_to_title) %>% 
        fct_relevel("Combined"),
      assigned_treatment_left = str_to_title(assigned_treatment_left)
    )) %>%
    unnest(est_takeup_te)  %>%
    select( 
      model_type,
      assigned_treatment = assigned_treatment_left, 
      assigned_dist = assigned_dist_group_left,
      mean_est, 
      per_0.05, 
      per_0.95
    )



signalling_te = dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) %>%
  select(model, model_type, model_name, est_takeup_te, fit_type) %>% 
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
        mu_assigned_treatment_left = str_to_title(mu_assigned_treatment_left),
      ),
    model_color = canva_pal(canva_palette_vibrant)(n())) %>%
    unnest(est_takeup_te)  %>%
    select(
      assigned_treatment = mu_assigned_treatment_left, 
      assigned_dist = assigned_dist_group_left, 
      model_type,
      mean_est,
      per_0.05,
      per_0.95
    )




private_te = dist_fit_data %>% 
  filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) %>%
  select(model, model_type, model_name, est_takeup_te, fit_type) %>% 
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
        assigned_treatment_left = str_to_title(assigned_treatment_left),
      ),
    model_color = canva_pal(canva_palette_vibrant)(n())) %>%
    unnest(est_takeup_te) %>%
    select(
      assigned_treatment = assigned_treatment_left,
      assigned_dist = assigned_dist_group_left, 
      model_type,
      mean_est, 
      per_0.05, 
      per_0.95
    )



te_df = bind_rows(
  incentive_te %>% mutate(estimand = "incentive"),
  private_te %>% mutate(estimand = "private"),
  signalling_te %>% mutate(estimand = "signal")
)
library(kableExtra)
te_df %>%
  kbl()

rf_te_input_df = te_df %>%
  filter(model_type == "reduced form") %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate(
    estim_value = linebreak(
      paste0(mean_est, "\n", "(", per_0.05, ", ", per_0.95, ")"), 
      align = "c"
      )
  ) %>%
  select(
    model_type, 
    assigned_treatment, 
    assigned_dist, 
    estim_value, 
    estimand
  ) %>%
  pivot_wider(
    names_from = assigned_dist, 
    values_from = estim_value
  )  %>%
  rename(
    rf_Close = Close, 
    rf_Far = Far, 
    rf_Combined = Combined
  )

structural_te_input_df = te_df %>%
  filter(model_type == "structural") %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  mutate(
    estim_value = linebreak(
      paste0(mean_est, "\n", "(", per_0.05, ", ", per_0.95, ")"), 
      align = "c"
      )
  ) %>%
  select(
    model_type, 
    assigned_treatment, 
    assigned_dist, 
    estim_value, 
    estimand
  ) %>%
  pivot_wider(
    id_cols = c(assigned_treatment, estimand),
    names_from = assigned_dist, 
    values_from = estim_value
  )  %>%
  mutate(across(
    c(Close, Far), 
    ~if_else(
      estimand == "private", 
      NA_character_, 
      .x
    )
  ))


te_input_df = left_join(
  structural_te_input_df,
  rf_te_input_df,
  by = c("estimand", "assigned_treatment")
)  %>%
  select(-model_type) %>%
  select(
    Estimand = estimand, 
    Treatment = assigned_treatment, 
    everything()
  ) %>%
  mutate(across(c(contains("rf_"), Close, Far), replace_na, "-")) %>%
  mutate(Estimand = factor(Estimand, levels = c("incentive", "signal", "private", ""))) %>%
  arrange(Estimand)



latex_group_gap_space = "2em"

te_kbl_df = te_input_df %>%
  mutate(
    Estimand = if_else(
      Treatment == "Ink", 
      Estimand, 
      factor("")
    )
  ) %>%
  select(-Estimand) %>%
  kbl(
    col.names = c(
      # "Estimand", 
      # "Treatment", 
      "",
      "Close", 
      "Far", 
      "Combined", 
      "Close", 
      "Far",
      "Combined"
    ), 
    format = "latex", 
    booktabs = TRUE, 
    escape = FALSE, 
    align = "lcccccc", 
    caption = "Caption goes here"
  ) %>%
  kable_styling(
    latex_options = c("scale_down")
  ) %>%
  add_header_above(
    c(
      " " = 1,
      "Structural" = 3, 
      "Reduced Form" = 3
      )
  ) %>%
  pack_rows(
    index = c(
      "Incentive" = 3, 
      "Signalling" = 3, 
      "Private" = 3
    ), 
    # latex_gap_space = latex_group_gap_space, 
    hline_after = TRUE, 
    # hline_before = TRUE,
    bold = TRUE
  ) %>%
  footnote(
    general = "The footnote goes here."
  )

te_kbl_df %>%
  save_kable(
    file.path(
      script_options$output_path,
      "te-table.tex"
    )
  )

