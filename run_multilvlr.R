#!/usr/bin/Rscript

if (!interactive()) withr::with_dir("..", source(file.path("packrat", "init.R")))

"Usage:
   run_multilvlr fit [options]
   run_multilvlr cv [options] [--folds=<folds> --sequential]
  
Options:
  --num-chains=<num-chains>, -c <num-chains>  Number of sampling chaings [default: 4]
  --iterations=<iterations>, -i <iterations>  Number of sampling iterations [default: 300]
  --warmup=<iterations>, -w <iterations>  Number of warmup iterations
  --adapt-delta=<adapt-delta> [default: 0.8]
  --max-treedepth=<max-treedepth> [default: 10]
  --analysis-data-only  Generate analysis data only; do not sample
  --debug-output  Run sampling for debug output
  --exclude-outliers  Exclude outlier observations from analysis data
  --impute-na  Impute unmeasured outcomes (not for counterfactuals but for households that have invalid or unmeasured outcomes)
  --always-predict
  --no-post-process  Do not post-process posterior samples
  --num-samples=<num-samples>, -s <num-samples>  Number of samples to post-process
  --hist-only  Calculate outcome histograms only
  --output-name=<output-name>, -o <output-name>  Name to use in stanfit .csv files and analysis data .RData file
  --no-db  Don't extract stanfit iterations to SQLite database
  --prior-predict  Run a prior posterior predictive distribution check
  --folds=<folds>  Number of folds in k-fold cross valiation [default: 10]" -> script_desc

script_options <- if (interactive()) {
  # docopt::docopt(script_desc, "-c 4 -i 4000 -w 2000")
  docopt::docopt(script_desc, "fit -c 2 -i 20 -w 10")
} else {
  docopt::docopt(script_desc)
}

# Setup -------------------------------------------------------------------

library(magrittr)
library(tidyverse)
library(parallel)
library(pbmcapply)
library(modelr)
library(rstan)
library(loo)

library(econometr)

source(file.path("multilvlr", "multilvlr_util.R"))

options(mc.cores = 12) # max(1, parallel::detectCores()))
rstan_options(auto_write = TRUE)

load(file.path("data", "analysis.RData"))

output_name <- script_options$`output-name`

if (is_null(output_name)) {
  output_name <- str_interp("takeup_model")
}

print(script_options)

# Prepare Data ------------------------------------------------------------

scale_normal_outcomes <- TRUE

all_ate <- #list(
  # tribble(
  #   ~ cluster_treatment_left,    ~ cluster_treatment_right,
  # 
  #   "calendar",                  "control",
  #   "ink",                       "control",
  #   "bracelet",                  "control",
  # ),
  
  tribble(
    ~ cluster_treatment_left,    ~ cluster_treatment_right, ~ dist_pot_treatment_left, ~ dist_pot_treatment_right,
    
    "ink",                       "control",                "close",                    "close", 
    "ink",                       "control",                "far",                      "far", 
    "calendar",                  "control",                "close",                    "close", 
    "calendar",                  "control",                "far",                      "far", 
    "bracelet",                  "control",                "close",                    "close", 
    "bracelet",                  "control",                "far",                      "far", 
    "bracelet",                  "ink",                    "close",                    "close", 
    "bracelet",                  "ink",                    "far",                      "far", 
    "bracelet",                  "calendar",               "close",                    "close", 
    "bracelet",                  "calendar",               "far",                      "far", 
  )
# ) 

exclude_outliers <- script_options$`exclude-outliers` %||% TRUE

get_stan_data <- function(.data, ...) { 
  .data %>% 
    mutate(cluster_treatment = assigned.treatment,
           dist_pot_treatment = dist.pot.group,
           cluster = cluster.id) %>% 
    prepare_bayesian_analysis_data(
     model_levels_metadata = tribble(
       ~ level_name,        ~ level_id_name,     ~ with_treatment_corr, ~ default_coef_scale, ~ default_coef_corr_lkj_df, ~ contained_in,                                ~ subgroup_containers,
       
       # "wave",              "wave_id",           TRUE,                 2,                     3,                          NA,                                            NA,
       # "county",            "county_id",         TRUE,                 0.8,                   3,                          "wave",                                        NA,
       # "cluster",           "cluster_id",        FALSE,                0.3,                   100,                        c("county", "wave"),                           NA,
       # "household",         "obs_index",         FALSE,                0,                     0,                          c("cluster", "county", "wave"),                c("cluster", "county"),
       # "cluster",           "cluster_id",        FALSE,                0.3,                   100,                        NA,                                            NA,
       # "household",         "obs_index",         FALSE,                0,                     0,                          c("cluster"),                                  c("cluster"), 
       "county",            "county_id",         TRUE,                 0.8,                   3,                          NA,                                            NA,
       "cluster",           "cluster_id",        FALSE,                0.3,                   100,                        c("county"),                                   NA,
       "household",         "obs_index",         FALSE,                0,                     0,                          c("cluster", "county"),                        "cluster",
      
       # TODO  Add subcounty 
       
     ) %>%  
       bind_rows(
         tribble(
           ~ level_name,             ~ with_treatment_corr, ~ default_coef_scale, ~ default_coef_corr_lkj_df, ~ covar_for_level,
           
           # Only for covar subgroups
           # "hh_covar_subgroup",      TRUE,                  1,                    3,                          "household",
           
         ) %>% 
           mutate(contained_in = NA)
       ),
     
     outcome_model_metadata = tribble( ####### ENDOGENOUS OUTCOMES #########
       ~ outcome_type,              ~ outcome_model_type,      ~ outcome_model_scaled, ~ variable_name,               ~ hyper_intercept_mean,  ~ exclude_outliers,
       
       "dewormed",                  "logit",                   FALSE,                  NA,                            0,                     FALSE,
       
     ) %>% 
       mutate(treatment_formula = lst(~ cluster_treatment * dist_pot_treatment),
              obs_level = "household",
              ate_pairs = rep(lst(all_ate), n()),
              variable_type = "modeled endogenous") %>%
      bind_rows(
        tribble(                     #####  MODELED COVARIATES ######
          ~ outcome_type,          ~ outcome_model_type,     ~ outcome_model_scaled, ~ variable_name,           ~ obs_level,
          
          # "phone_owner",          "logit",                   FALSE,                  NA,                        "household",

        ) %>% 
          mutate(variable_type = "modeled exogenous"),

        tribble(                    ##### UNMODELED COVARIATES ######
          ~ outcome_type,                ~ obs_level, ~variable_name, ~ calculator,

        ) %>% 
          mutate(variable_type = "unmodeled exogenous")
      ) %>%
       mutate(
        treatment_outcome_sigma_scale = if_else(!outcome_model_type %in% c("logit", "ordered_logit"), 1, 0), # Assuming all normal outcomes are scaled
        hyper_param_df = 7,
        hyper_intercept_mean = 0,
        hyper_intercept_scale = 0.75,
        hyper_treatment_coef_scale = 1.25),
    
      # Cross Validation 
      cv_level = if (script_options$cv) "cluster" else NA, 
      excluded_ids = array(dim = 0),
      num_excluded_ids = 0,
      
      debug_output = as.integer(script_options$`debug-output` %||% FALSE),
      
      estimate_ate = !(script_options$`no-ate` %||% FALSE),
      
      impute_na = script_options$`impute-na` %||% FALSE,
     
      run_type = if (script_options$`prior-predict` %||% FALSE) "prior_predict" else "fit", 
     
      obs_effects_scale = 0.5,
      obs_effects_corr_lkj_df = 2, 
     
      use_obs_effects = !(script_options$`no-obs-effects` %||% FALSE),
     
      always_predict = script_options$`always-predict`,
      
      ...
    )
}

stan_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  get_stan_data()

num_chains <- if (script_options$`debug-output` %||% FALSE) 2 else as.integer(script_options$`num-chains` %||% 2)

if (!interactive()) {
  if (!script_options$`debug-output`) {
    save(stan_data, file = file.path("stan_analysis_data", str_c(output_name, ".RData")))
  }
  
  if (script_options$`analysis-data-only`) {
    cat(str_interp("Output: ${output_name}\n"))
    quit()
  }
}

 # Sampling ----------------------------------------------------------------

multilvlr_model <- stan_model(file.path("multilvlr", "multilvlr.stan"))

if (script_options$fit) {
  param_to_save <- c(
    # "hyper_predictor",
     "iter_model_level_predictor_with_containers",
     "iter_level_mean", "iter_level_quantiles", "iter_te_mean", "iter_te_quantiles", "iter_model_level_treatment_residuals",
     "iter_model_level_treatment_residual_variance", "iter_model_level_te_residuals", "iter_model_level_te_residual_variance",
     "log_lik"
     # "obs_effects_tau"
     # "iter_level_ecdf", "iter_te_ecdf_diff"
  )
  
  model_fit <- sampling(
    multilvlr_model,
    data = stan_data, 
    chains = num_chains,
    iter = if (script_options$`debug-output` %||% FALSE) 2 else as.integer(script_options$iterations %||% 20),
    warmup = as.integer(script_options$warmup) %||% (as.integer(script_options$iterations %||% 300) %/% 2),
    control = lst(max_treedepth = as.integer(script_options$`max-treedepth` %||% 10),
                  adapt_delta = as.numeric(script_options$`adapt-delta` %||% 0.8)), 
    include = TRUE,
    pars = param_to_save,
    sample_file = if (script_options$`debug-output` %||% TRUE) NULL else file.path("stanfit", str_c(output_name, ".csv")), 
    save_warmup = FALSE)
  
  check_hmc_diagnostics(model_fit)
} else {
  param_to_save <- c("log_lik", "log_lik_heldout")
  
  kfold_groups <- kfold_split_stratified(K = as.integer(script_options$folds), 
                                         x = stan_data$prepared_analysis_data %>% 
                                           distinct(cluster_id, assigned.treatment, dist.pot.group) %>% 
                                           group_indices(assigned.treatment, dist.pot.group))
  
  log_lik_list <- map(seq(as.integer(script_options$folds)), ~ which(kfold_groups == .)) %>% 
    pbmclapply(function(excluded_clusters, stan_data, multilvlr_model) {
      cat("Excluding: ", excluded_clusters)
      
      curr_stan_data <- stan_data %>%
        list_modify(
          excluded_ids = excluded_clusters,
          num_excluded_ids = length(excluded_clusters))
      
    sampling(
      multilvlr_model,
      data = curr_stan_data, 
      chains = num_chains,
      iter = if (script_options$`debug-output` %||% FALSE) 2 else as.integer(script_options$iterations %||% 20),
      warmup = as.integer(script_options$warmup) %||% (as.integer(script_options$iterations %||% 300) %/% 2),
      control = lst(max_treedepth = as.integer(script_options$`max-treedepth` %||% 10),
                    adapt_delta = as.numeric(script_options$`adapt-delta` %||% 0.8)), 
      include = TRUE,
      pars = param_to_save,
      sample_file = NULL, 
      save_warmup = FALSE) %>% 
      extract_log_lik("log_lik_heldout") 
  },
  stan_data = stan_data,
  multilvlr_model = multilvlr_model,
  mc.silent = !script_options$sequential,
  mc.cores = if (script_options$sequential) 1 else 3)
  
  takeup_kfold <- tryCatch(
    kfold(log_lik_list), 
    error = function(err) { 
      print(err)
      return(log_lik_list)
    })
}

# Post processing ---------------------------------------------------------

if (script_options$fit) {
  if (!(script_options$`no-db` %||% FALSE) && require(dbplyr)) {
    db_filename <- file.path("stan_analysis_data", str_c(output_name, ".sqlite"))
    
    suppressWarnings(file.remove(db_filename))
    
    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_filename)
    
    extract_model_fit_to_db(db_src = con, stan_data = stan_data, model_fit = model_fit, verbose = TRUE)
  }
  
  if (!script_options$`no-post-process`) {
    if (!is_null(script_options$`num-samples`) && script_options$`no-db`) {
      num_samples <- as.integer(script_options$`num-samples`)
  
      cat("restricting to", num_samples, "samples...done.\n")
    } else {
      num_samples <- NULL
    }
    
    if (!script_options$`no-db`) {
      model_fit <- NULL
      db_src <- src_dbi(con, auto_disconnect = TRUE) 
    } else {
      db_src <- NULL
    }
    
    processed_data <- postprocess_model_fit(model_fit, stan_data, db_src, num_samples, verbose = TRUE, summarize_est = TRUE, 
                                            hist_only = script_options$`hist-only` %||% FALSE, no_hist = TRUE,
                                            treatment_variable = c("cluster_treatment", "dist_pot_treatment"))
    
    cat("Saving post-processed data...")
    
    processed_data %>% 
      as.environment() %>% 
      save(list = names(.), envir = .,
           file = file.path("stan_analysis_data", str_interp("${output_name}_processed.RData")))
    
    cat("done.\n")
  }
  
  if (!(script_options$`debug-output` %||% FALSE)) {
    cat(str_interp("Output: ${output_name}\n"))
  }
} else {
  save(takeup_kfold, 
       file = file.path("stan_analysis_data", str_interp("${output_name}_cv.RData")))
}
