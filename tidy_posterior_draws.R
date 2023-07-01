#!/usr/bin/Rscript
script_options <- docopt::docopt(
  stringr::str_glue(
"Usage:
  tidy_posterior_draws.R <fit-version> [options] [<params>...]
  
Options:
  --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
  --output-path=<path>  Path to find results [default: temp-data]
  --model=<model>  Which model to postprocess
  --full-inputname=<full-inputname>  Input name to use.
  --prior  Postprocess the prior predictive
  --chain=<chain>  Which chains to tidy [default: 1-4]
  --exclude-params=<exclude-params>  Parameters to exclude [default: cluster_error]
  "), 
  args = if (interactive()) "
  95
  --input-path=temp-data
  --output-path=temp-data
  --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP_FOB
  --chain=1-4
  --exclude-params=cluster_error

  " else commandArgs(trailingOnly = TRUE)
)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(tidybayes)

fit_type_str = if_else(script_options$prior, "prior", "fit")

all_params = c(
    "ates",
    "levels",
    "cluster_error",
    "rep_return_dist_draws",
    "prop_draws",
    "roc_draws",
    "belief_ates",
    "belief_probs",
    "wtp_params",
    "sm_draws"
)
params_we_want = if (length(script_options$params) == 0) all_params else script_options$params

exclude_params = if (!is.null(script_options$exclude_params)) str_split_1(script_options$exclude_params, ",") else ""

params_we_want = params_we_want[!(params_we_want %in% exclude_params)]

if (is.null(script_options$full_inputname)) {
    input_filename = str_glue(
        "rvar_processed_dist_{fit_type_str}{script_options$fit_version}_{params_we_want}_{script_options$model}_{script_options$chain}.rds"
    )
} else {
    input_filename = script_options$full_inputname
}



param_df = tibble(
    input_path = file.path(
        script_options$input_path,
        input_filename
    ),
    version = script_options$fit_version,
    model = script_options$model,
    param = params_we_want,
    fit_type = if_else(script_options$prior, "prior-predict", "fit")
)

param_df = param_df %>%
    mutate(
        draws = map(input_path, read_rds)
    ) %>%
    select(-input_path) %>%
    mutate(draws = map(draws, ~ {
        .x %>%
            ungroup() %>%
            select(-fit_version, -model, -fit_type)
    }))


param_df = param_df %>%
    mutate(
        tidy_draws = map(
            draws, 
            ~median_qi(.x, value, .width = c(0.95, 0.9, 0.8, 0.5), na.rm = TRUE) %>%
                to_broom_names()
        ),
    )

param_df %>%
    saveRDS(
        file.path(
            script_options$output_path,
            str_glue(
                "tidy_processed_dist_{fit_type_str}{script_options$fit_version}_{script_options$model}_{script_options$chain}.rds"
            )

        )
    )


# wide_struct_tes = all_tes %>% 
#     filter(estimand == "overall") %>%
#     select(
#       dist_treatment,
#       dist_group,
#       pr_takeup
#     ) %>%
#     pivot_wider(
#       names_from = dist_group,
#       values_from = pr_takeup
#     ) %>%
#     select(dist_treatment, combined, close, far) %>%
#     arrange(dist_treatment) %>%
#     bind_rows(
#       # bracelet minus calendar row
#       all_tes %>%
#         filter(dist_treatment %in% c("bracelet", "calendar")) %>%
#         filter(estimand == "overall") %>%
#         pivot_wider(names_from = dist_treatment, values_from = pr_takeup) %>%
#         mutate(
#           bracelet_minus_calendar = bracelet - calendar
#         ) %>%
#         select(dist_group, bracelet_minus_calendar) %>%
#         pivot_wider(
#           names_from = dist_group,
#           values_from = bracelet_minus_calendar
#         ) %>%
#         mutate(dist_treatment = "bracelet_minus_calendar")
#     )