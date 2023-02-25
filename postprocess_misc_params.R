
#!/usr/bin/Rscript
#
# This script is used to postprocess Stan fit for the various models, reduced form and structural. In addition to putting our analysis in a format that 
# allows for easy extraction of all levels and treatment effects, it allows handles imputing take-up levels for counterfactuals using Stan-generated cost-benefits 
# and reputational returns parameters: using these we can calculate the probability of take-up after calculating the v^* fixed point solution.
#

script_options <- docopt::docopt(
  stringr::str_glue(
"Usage:
  postprocess_misc_params.R <fit-version> [--full-outputname --cores=<num-cores> --output-path=<path> --input-path=<path> --load-from-csv --no-prior --no-rate-of-change --keep-fit --models=<models>]
  
Options:
  --cores=<num-cores>  Number of cores to use [default: 12]
  --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
  --output-path=<path>  Path to find results [default: temp-data]
  --keep-fit "), 
  # args = if (interactive()) "29" else commandArgs(trailingOnly = TRUE)
  # args = if (interactive()) "30" else commandArgs(trailingOnly = TRUE)
  # args = if (interactive()) "test --full-outputname --load-from-csv --cores=1" else commandArgs(trailingOnly = TRUE)
  # args = if (interactive()) "31 --cores=6" else commandArgs(trailingOnly = TRUE) 
  # args = if (interactive()) "test --full-outputname --cores=1 --input-path=/tigress/kn6838/takeup --output-path=/tigress/kn6838/takeup" else commandargs(trailingonly = true)
  args = if (interactive()) "71  --cores=1 --models=STRUCTURAL_LINEAR_U_SHOCKS  --load-from-csv " else commandArgs(trailingOnly = TRUE)
)

library(magrittr)
library(tidyverse)
library(rlang)
library(cmdstanr)
library(tidybayes)
library(furrr)
library(posterior)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

fit_version <- script_options$fit_version
postprocess_cores <- as.integer(script_options$cores)

if (postprocess_cores > 1) {
  if (interactive()) {
    plan(multisession, workers = postprocess_cores)
  } else {
    plan(multicore, workers = postprocess_cores)
  }
} else {
  plan(sequential)
}

# # Analysis Data --------------------------------------------------------------------

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

# analysis_data <- monitored_nosms_data


string_to_match = str_glue(
    "dist_fit{script_options$fit_version}_{script_options$models}-\\d+\\.csv$"
)

cmdstan_files = fs::dir_ls(
    script_options$input_path,
    regexp = string_to_match 
)

if (interactive()) {
    cmdstan_files = cmdstan_files[1]
}

cmdstan_fit = as_cmdstan_fit(cmdstan_files)


var_names = get_variables(cmdstan_fit) %>%
    str_extract("\\D+") %>%
    str_remove_all("\\[") %>%
    unique()



wtp_rvars = spread_rvars(cmdstan_fit, hyper_wtp_mu, wtp_sigma)

rvar_pnorm = rfun(pnorm)

wtp_rvars = wtp_rvars %>%
    mutate(
        mu_std = hyper_wtp_mu / wtp_sigma, 
        phi_mu_std = rvar_pnorm(mu_std)
        )

wtp_rvars %>%
    select(phi_mu_std) %>%
    median_qi(.width = 0.9) %>%
    to_broom_names() %>%
    mutate(variable = "wtp_phi_mu_std") %>%
    write_csv(
        file.path(
            script_options$output_path, 
            str_interp(
                "wtp_phi_summary_draws_dist_fit${fit_version}.csv"
            )
        )
    )



draw_rvars = gather_rvars(cmdstan_fit, wtp_sigma)


rm(cmdstan_fit)
gc()


draw_pis = draw_rvars %>%
    point_interval(
        .value,
        .width = c(0.95, 0.9, 0.8, 0.5)
    ) 

draw_point_summ = draw_rvars %>%
    group_by(variable = .variable) %>%
    summarise(
        mean_est = mean(.value), 
        per_0.5 = median(.value)
    )
    
draw_per_df = draw_pis %>%
    pivot_wider(
        id_cols = .variable, 
        names_from = .width, 
        values_from = .lower:.upper
    )

per_colnames = c(
    "variable",
    "per_0.025",
    "per_0.05",
    "per_0.1",
    "per_0.25",
# upper
    "per_0.975",
    "per_0.95",
    "per_0.90",
    "per_0.75"
)

colnames(draw_per_df) = per_colnames



clean_draw_output = left_join(
    draw_point_summ,
    draw_per_df,
    by = "variable"
)


clean_draw_output %>%
    write_csv(
        file.path(
            script_options$output_path,
            str_interp("misc_processed_dist_fit${fit_version}.csv")
        )
    )

# var_names

# cluster_rep_return_rvar_df = gather_rvars(cmdstan_fit, cluster_rep_return_dist[i,j,k])

#   cluster_rep_return_dist_control = dist_fit_data %>%
#     select(cluster_rep_return_dist) %>%
#     unnest(cols = c(cluster_rep_return_dist)) %>%
#     filter(assigned_treatment == "control")

#   cluster_rep_return_te_df = dist_fit_data %>%
#       select(cluster_rep_return_dist) %>%
#       unnest(cols = c(cluster_rep_return_dist))  %>%
#       # filter(assigned_treatment != "control")  %>%
#       unnest(iter_data) %>%
#       left_join(
#         cluster_rep_return_dist_control %>% 
#           select(roc_distance_index, iter_data) %>%
#           unnest(iter_data) %>%
#           rename(iter_est_right = iter_est),
#         by = c("roc_distance_index", "iter_id")
#       ) %>%
#       mutate( 
#         iter_est_te = iter_est - iter_est_right
#       ) %>%
#       nest(iter_data = c(iter_id, iter_est_te, iter_est, iter_est_right)) %>%
#       select(-contains("per"), -mean_est) %>%
#       mutate(
#         mean_est = map_dbl(iter_data, ~ mean(.$iter_est_te)),
#         takeup_quantiles = map(iter_data, quantilize_est, iter_est_te, wide = TRUE, quant_probs = c(quant_probs))
#       ) %>% 
#       unnest(takeup_quantiles)

#   # cluster_rep_return_te_df %>%
#   #   select(-iter_data) %>%
#   #   write_csv(file.path(script_options$output_path, str_interp("processed_rep_return_dist_fit${fit_version}.csv")))
