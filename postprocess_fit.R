#!/usr/bin/Rscript

if (!interactive()) withr::with_dir("..", source(file.path("packrat", "init.R")))

script_options <- if (interactive()) {
  NULL
} else {
  docopt::docopt(
    
"Usage:
   postprocess_fit <output-name> [-s <num-samples> -c <chains> --outcomes=<outcomes> --use-stanfit --no-rhats --summarize --hist-only]"

)
}

print(script_options)

library(magrittr)
library(tidyverse)

source(file.path("multilvlr", "multilvlr_util.R"))

output_name <- script_options$`output-name`

cat("Loading pre-processed analysis data...")
load(file.path("stan_analysis_data", str_c(output_name, ".RData")))
cat("done\n")

if (script_options$`use-stanfit`) {
  library(rstan)
  
  cat(str_interp("Loading stanfit CSV (${chains} chains)..."))
  
  chains <- script_options$chains %||% 8
  
  model_fit <- dir("stanfit", pattern = str_interp("${output_name}_[1-${chains}]\\.csv"), full.names = TRUE) %>%
    read_stan_csv()
  
  cat("done\n")
  
  check_hmc_diagnostics(model_fit)
  
  if (!script_options$`no-rhats`) {
    cat("Extracting R hats...")
    
    model_fit_rhats <- try(bayesplot::rhat(model_fit))
    summary(model_fit_rhats)
    
    cat("done\n")
  } else {
    model_fit_rhats <- NULL
  }
   
  if (!is_null(script_options$`num-samples`)) {
    num_samples <- as.integer(script_options$`num-samples`)
  
    cat("restricting to ", num_samples, "samples...")
  } else {
    num_samples <- NULL
  }
  
  db_src <- NULL
  
  cat("done\n")
} else {
  library(dbplyr)
  cat("Connecting to iterations database...")
  
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = file.path("stan_analysis_data", str_c(output_name, ".sqlite")), flags = RSQLite::SQLITE_RO)
  db_src <- src_dbi(con, auto_disconnect = TRUE)
  
  model_fit <- NULL
  
  cat("done\n")
}

if (!is_null(script_options$outcomes)) {
  outcomes <- script_options$outcomes %>% str_split(",") %>% unlist()

  cat("[ Restricting analysis to outcomes", str_c(outcomes, collapse = ", "), "]\n")
} else {
  outcomes <- NULL
}

processed_data <- postprocess_model_fit(model_fit, stan_data, db_src, num_samples, verbose = TRUE, summarize_est = script_options$summarize %||% FALSE, 
                                        hist_only = script_options$`hist-only` %||% FALSE,
                                        treatment_variable = c("cluster_treatment", "dist_pot_treatment"),
                                        outcomes = outcomes) 

cat("Saving post-processed data...")

if (!script_options$`no-rhats` && script_options$`use-stanfit`) {
  processed_data %<>% update_list(model_fit_rhats) 
}

processed_data %>% 
  as.environment() %>% 
  save(list = names(.), envir = .,
       file = file.path("stan_analysis_data", str_interp("${output_name}_processed.RData")))
cat("done\n")