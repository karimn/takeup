#!/usr/bin/Rscript

if (!interactive()) withr::with_dir("..", source(file.path("packrat", "init.R")))

script_options <- docopt::docopt(
"Usage:
  run_stan_dist fit [--structural=<structural type> --no-save --sequential --iter=<iter> --thin=<thin> --force-iter --models=<models> --outputname=<output file name> --update-output --no-dist-cluster-effects --model-file=<stan model>]
  run_stan_dist cv [--structural=<structural type> --folds=<number of folds> --no-save --sequential --iter=<iter> --thin=<thin> --force-iter --models=<models> --outputname=<output file name> --no-dist-cluster-effects --model-file=<stan model>]
  
Options:
  --model-file<stan model>  Stan model file to use [default: dist_spline3.stan]
  --folds=<number of folds>  Cross validation folds [default: 10]
  --iter=<iter>  Number of (warmup + sampling) iterations [default: 8000]
  --thin=<thin>  Thin samples [default: 1]
  --structural=<structural type>  Specify structural model type [default: 0]", 

  args = if (interactive()) "fit --structural=1 --sequential --iter=4000 --thin=1 --force-iter --models=6 --model-file=dist_splines4.stan" else commandArgs(trailingOnly = TRUE) 
) 

library(magrittr)
library(tidyverse)
library(parallel)
library(pbmcapply)
library(HRW)
library(splines2)
library(rstan)
library(loo)

options(mc.cores = 12)
rstan_options(auto_write = TRUE)

folds <- as.integer(script_options$folds %||% 10)
iter <- as.integer(script_options$iter)
output_name <- if (!is_null(script_options$outputname)) { script_options$outputname } else if (script_options$fit) { "dist_fit" } else { "dist_kfold" }
output_file_name <- file.path("stan_analysis_data", str_c(output_name, ".RData"))
use_dist_cluster_effects <- !(script_options$`no-dist-cluster-effects` %||% FALSE)
models_to_run <- if (!is_null(script_options$models)) c(str_split(script_options$models, ",", simplify = TRUE)) %>% as.integer()
structural_type <- as.integer(script_options$structural) 
thin_by <- as.integer(script_options$thin)
model_file <- file.path("stan_models", script_options$`model-file`)

source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))

# Data --------------------------------------------------------------------

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (. - mean(.)) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original) + mean(original)

monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot),
         cluster_id = group_indices(., cluster.id))

cluster_analysis_data <- monitored_nosms_data %>% 
  group_by(cluster.id, assigned.treatment, cluster.dist.to.pot, dist.pot.group) %>% 
  summarize(prop_takeup = mean(dewormed)) %>% 
  ungroup() %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot),
         standard_prop_takeup = standardize(prop_takeup))


# Functions ---------------------------------------------------------------

get_spline_range <- function(x) {
  lst(lower_range = 1.01 * min(x) - 0.01 * max(x),
      upper_range = 1.01 * max(x) - 0.01 * min(x))
}

calculate_splines <- function(x, num_interior_knots, splines_for = x, spline_type = c("osullivan", "i-spline", "b-spline"), ...) {
  spline_type <- match.arg(spline_type)
  
  spline_range <- get_spline_range(x)
  
  interior_knots <- quantile(unique(x), seq(0, 1, length = num_interior_knots + 2)) %>% 
    magrittr::extract(2:(num_interior_knots + 1))

  switch(
    spline_type,
    
    "osullivan" = ZOSull(splines_for, range.x = unlist(spline_range), intKnots = interior_knots, ...),
    "i-spline" = iSpline(splines_for, knots = interior_knots,  Boundary.knots = unlist(spline_range), ...),
    "b-spline" = bSpline(splines_for, knots = interior_knots,  Boundary.knots = unlist(spline_range), ...)
  )
}

extract_sim_level <- function(fit, par, grid_dist) {
  analysis_data <- monitored_nosms_data
  
  fit %>% 
    as.data.frame(par = par) %>% 
    mutate(iter_id = seq_len(n())) %>% 
    gather(indices, iter_est, -iter_id) %>% 
    tidyr::extract(indices, c("grid_index", "assigned_treatment"), "(\\d+),(\\d+)", convert = TRUE) %>% 
    group_nest(grid_index, assigned_treatment, .key = "iter_data") %>% 
    mutate(mean_est = map_dbl(iter_data, ~ mean(.$iter_est)),
           quantiles_est = map(iter_data, 
                               function(iter, probs) 
                                 quantile(iter$iter_est, probs = probs, names = FALSE) %>% 
                                 enframe(name = NULL, value = "est") %>% 
                                 mutate(per = probs),
                               probs = c(0.05, 0.1, 0.5, 0.9, 0.95))) %>% 
    left_join(tibble(grid_dist) %>% mutate(grid_index = seq_len(n())), by = "grid_index") %>% 
    mutate(assigned_treatment = factor(assigned_treatment, labels = levels(analysis_data$assigned.treatment)),
           grid_dist = unstandardize(grid_dist, analysis_data$cluster.dist.to.pot)) 
}

extract_obs_fit_level <- function(fit, par, stan_data) {
  analysis_data <- monitored_nosms_data
  
  fit %>% 
    as.data.frame(par = par) %>% 
    mutate(iter_id = seq_len(n())) %>% 
    gather(indices, iter_est, -iter_id) %>% 
    tidyr::extract(indices, "obs_index", "(\\d+)", convert = TRUE) %>% 
    group_nest(obs_index, .key = "iter_data") %>% 
    mutate(mean_est = map_dbl(iter_data, ~ mean(.$iter_est)),
           quantiles_est = map(iter_data, 
                               function(iter, probs) 
                                 quantile(iter$iter_est, probs = probs, names = FALSE) %>% 
                                 enframe(name = NULL, value = "est") %>% 
                                 mutate(per = probs),
                               probs = c(0.05, 0.1, 0.5, 0.9, 0.95)),
           assigned_treatment = stan_data$assigned_treatment,
           assigned_dist = unstandardize(stan_data$standard_dist, analysis_data$cluster.dist.to.pot)) 
}

extract_sim_diff <- function(level) {
  level %>% 
    select(-mean_est, -quantiles_est) %>% 
    rename(assigned_treatment_left = assigned_treatment) %>% 
    mutate(assigned_treatment_right = "control") %>% 
    left_join(select(., -grid_dist, -assigned_treatment_right), 
              by = c("assigned_treatment_right" = "assigned_treatment_left", "grid_index"), 
              suffix = c("_left", "_right")) %>% 
    filter(assigned_treatment_left != assigned_treatment_right) %>% 
    mutate(iter_diff = map2(iter_data_left, iter_data_right, inner_join, by = "iter_id", suffix = c("_left", "_right")) %>% 
             map(mutate, iter_est = iter_est_left - iter_est_right)) %>% 
    mutate(mean_est = map_dbl(iter_diff, ~ mean(.$iter_est)),
           quantiles_est = map(iter_diff, 
                                       function(iter, probs) 
                                         quantile(iter$iter_est, probs = probs, names = FALSE) %>% 
                                         enframe(name = NULL, value = "est") %>% 
                                         mutate(per = probs),
                                                # interval_size = round(2 * abs(probs - 0.5), 4)),
                                       probs = c(0.05, 0.1, 0.5, 0.9, 0.95))) 
}

stan_list <- function(models_info, stan_data) {
  get_sample_file_name <- . %>% 
    sprintf("dist_model_%s.csv") %>% 
    file.path("stanfit", .) 
  
  dist_model <- stan_model(model_file)
  
  if (script_options$fit) {
    inner_sampler <- function(curr_model, dist_model, stan_data) {
      curr_stan_data <- stan_data %>%
        list_modify(!!!curr_model) %>%
        map_if(is.factor, as.integer)
      
      browser() 
      
      fit <- sampling(
        dist_model,
        iter = if (script_options$`force-iter`) iter else (curr_stan_data$iter %||% iter),
        thin = thin_by,
        chains = 4,
        control = lst(adapt_delta = 0.9),
        save_warmup = FALSE,
        refresh = 100,
        init = curr_model$init %||% "random",
        data = curr_stan_data)
      
      return(fit)
    }
    
    if (script_options$sequential) {
      models_info %>% 
        map(inner_sampler,
            dist_model = dist_model,
            stan_data = stan_data)
      
    } else {
      models_info %>% 
        # future_imap(function(curr_model, model_name, dist_model) {
        pbmclapply(inner_sampler,
                   dist_model = dist_model,
                   stan_data = stan_data,
                   ignore.interactive = TRUE,
                   mc.silent = TRUE,
                   mc.cores = 3)
    }
    
    models_info %>% 
      # imap(function(curr_model, model_name, dist_model) {
      # future_imap(function(curr_model, model_name, dist_model) {
      pbmclapply(,
      dist_model = dist_model,
      stan_data = stan_data,
      ignore.interactive = TRUE,
      mc.silent = !script_options$sequential,
      mc.cores = if (script_options$sequential) 1 else 3)
      # .progress = TRUE)
  } else if (script_options$cv) {
    models_info %>% 
      imap(function(curr_model, model_name, dist_model, stan_data) {
        kfold_groups <- kfold_split_stratified(K = folds, x = stan_data$cluster_assigned_dist_group_treatment)
        
        log_lik_list <- map(seq(folds), ~ which(kfold_groups == .)) %>% 
          pbmclapply(function(excluded_clusters, model_name, dist_model, stan_data) {
            curr_stan_data <- stan_data %>%
              list_modify(!!!curr_model,
                          excluded_clusters = excluded_clusters,
                          num_excluded_clusters = length(excluded_clusters)) %>%
              map_if(is.factor, as.integer)
            
            curr_iter <- if (script_options$`force-iter`) iter else (curr_stan_data$iter %||% iter)
            curr_chains <- 4
            
            sampling(dist_model,
                     iter = curr_iter,
                     thin = thin_by, # (curr_chains * curr_iter) %/% 1000,
                     chains = curr_chains,
                     control = lst(adapt_delta = 0.9),
                     save_warmup = FALSE,
                     refresh = 1000,
                     init = curr_model$init %||% "random",
                     pars = "log_lik_heldout", 
                     data = curr_stan_data) %>% 
              extract_log_lik("log_lik_heldout") 
          },
          dist_model = dist_model,
          stan_data = stan_data,
          model_name = model_name,
          ignore.interactive = TRUE,
          mc.silent = !script_options$sequential,
          mc.cores = if (script_options$sequential) 1 else 3)
          
        return(tryCatch(kfold(log_lik_list), 
                        error = function(err) { 
                          print(err)
                          return(log_lik_list)
                        }))
      },
      dist_model = dist_model,
      stan_data = stan_data) 
  }
}

# Splines -----------------------------------------------------------------

num_interior_knots <- 100

cluster_standard_dist <- distinct(monitored_nosms_data, cluster_id, standard_cluster.dist.to.pot) %>% 
  arrange(cluster_id) %>% 
  pull(standard_cluster.dist.to.pot)

Z_osullivan <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, spline_type = "osullivan")
Z_i_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, spline_type = "i-spline")
Z_b_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, spline_type = "b-spline")

grid_dist <- get_spline_range(cluster_standard_dist) %>% unname() %>% list_modify(length = 1001) %>% do.call(seq, .)

Z_grid_osullivan <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, splines_for = grid_dist, spline_type = "osullivan")
Z_grid_i_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, splines_for = grid_dist, spline_type = "i-spline")
Z_grid_b_spline <- calculate_splines(cluster_standard_dist, num_interior_knots = num_interior_knots, splines_for = grid_dist, spline_type = "b-spline")



# Models ------------------------------------------------------------------

num_treatments <- n_distinct(monitored_nosms_data$assigned.treatment)
num_clusters <- n_distinct(monitored_nosms_data$cluster_id)

generate_initializer <- function(base_init = function() lst(beta_cluster_sd = abs(rnorm(num_treatments))), 
                                 structural_type = 0) {
  base_list <- base_init()
  
  if (structural_type > 0 || !is_empty(base_list)) {
    function() {
      if (structural_type > 0) {
        base_list %>% 
          list_modify(
            mu_rep = abs(rnorm(num_treatments)),
            mu_rep_raw = abs(rnorm(num_treatments)),
            structural_beta = rnorm(num_treatments),
            # dist_cost_k = abs(rnorm(num_treatments)),
            dist_cost_k = abs(rnorm(1)),
            structural_beta_cluster_raw = matrix(rnorm(num_clusters * num_treatments), nrow = num_clusters, ncol = num_treatments),
            structural_beta_cluster = matrix(rnorm(num_clusters * num_treatments), nrow = num_clusters, ncol = num_treatments),
            structural_beta_cluster_sd = abs(rnorm(num_treatments)),
            
            cluster_takeup_prob = runif(num_clusters),
            structural_cluster_takeup_prob = runif(num_clusters)
          )
                      
      } else {
        return(ret_list) 
      }
    }
  } else return(NULL)
}

models <- lst(
  NO_DIST = lst(model_type = 1),
  
  DISCRETE_DIST = lst(model_type = 2),
  
  LINEAR_DIST = lst(model_type = 3,
                    use_dist_cluster_effects = FALSE),
  
  QUADRATIC_NONLINEAR_DIST = lst(model_type = 4,
                                 use_dist_cluster_effects = FALSE),
  
  CUBIC_NONLINEAR_DIST = lst(model_type = 5,
                             use_dist_cluster_effects = FALSE),
  
  SEMIPARAM_NONLINEAR_DIST_OSULLIVAN = lst(model_type = 6,
                                           use_dist_cluster_effects = FALSE,
                                           init = generate_initializer(structural_type = structural_type)),
                                           # init = function() lst(u_splines_cluster_v_sigma = abs(rnorm(num_clusters, 0, 0.25)))),
  
  SEMIPARAM_NONLINEAR_DIST_BSPLINE = lst(model_type = 7, 
                                         Z_splines_v = Z_b_spline,
                                         num_knots_v = ncol(Z_splines_v),
                                         Z_grid_v = Z_grid_b_spline,
                                         iter = 16000,
                                         init = function() lst(u_splines_cluster_v_sigma = abs(rnorm(num_clusters, 0, 0.25)))),
  
  SEMIPARAM_NONLINEAR_DIST_ISPLINE = lst(model_type = 8,
                                         Z_splines_v = Z_i_spline,
                                         num_knots_v = ncol(Z_splines_v),
                                         Z_grid_v = Z_grid_i_spline,
                                         init = function() lst(u_splines_cluster_v_sigma = abs(rnorm(num_clusters, 0, 0.25)))),
  
  LINEAR_DIST_CE = lst(model_type = 3,
                    use_dist_cluster_effects = TRUE,
                    init = function() lst(dist_beta_cluster_linear_sd = abs(rnorm(num_treatments, 0, 0.25)))),
  
  QUADRATIC_NONLINEAR_DIST_CE = lst(model_type = 4,
                                 use_dist_cluster_effects = TRUE,
                                 init = function() lst(dist_beta_cluster_linear_sd = abs(rnorm(num_treatments, 0, 0.25)),
                                                       dist_beta_cluster_quadratic_sd = abs(rnorm(num_treatments, 0, 0.25)))),
  
  CUBIC_NONLINEAR_DIST_CE = lst(model_type = 5,
                             use_dist_cluster_effects = TRUE,
                             init = function() lst(dist_beta_cluster_linear_sd = abs(rnorm(num_treatments, 0, 0.25)),
                                                   dist_beta_cluster_quadratic_sd = abs(rnorm(num_treatments, 0, 0.25)),
                                                   dist_beta_cluster_cubic_sd = abs(rnorm(num_treatments, 0, 0.25)))),
)


# Stan Run ----------------------------------------------------------------

stan_data <- lst(
  num_obs = nrow(monitored_nosms_data),
  num_grid_obs = length(grid_dist),
  num_treatments = 4,
  num_clusters = n_distinct(monitored_nosms_data$cluster_id),
  obs_cluster_id = monitored_nosms_data$cluster_id,
  cluster_assigned_treatment = distinct(monitored_nosms_data, cluster_id, assigned.treatment) %>% 
    arrange(cluster_id) %>% 
    pull(assigned.treatment),
  takeup = monitored_nosms_data$dewormed,
  cluster_standard_dist = distinct(monitored_nosms_data, cluster_id, standard_cluster.dist.to.pot) %>% 
    arrange(cluster_id) %>% 
    pull(standard_cluster.dist.to.pot),
  cluster_assigned_dist_group_treatment = distinct(monitored_nosms_data, cluster_id, assigned.treatment, dist.pot.group) %>% 
    arrange(cluster_id) %>% 
    group_indices(assigned.treatment, dist.pot.group),
  num_dist_group_treatments = n_distinct(cluster_assigned_dist_group_treatment),
  grid_dist,
  
  Z_splines_v = Z_osullivan,
  num_knots_v = ncol(Z_splines_v),
  
  Z_grid_v = Z_grid_osullivan,
  
  u_splines_v_sigma_sd = 1,
 
  num_excluded_clusters = 0,
  excluded_clusters = array(dim = 0),
  
  use_binomial = 0,
  use_cluster_effects = 1,
  use_dist_cluster_effects,
  is_structural = structural_type,
) %>% 
  list_modify(!!!map(models, pluck, "model_type") %>% set_names(~ str_c("MODEL_TYPE_", .)))

models <- if (!is_null(models_to_run)) {
  models[models_to_run]
} else {
  models[-c(7,8)]
} 

if (script_options$fit) {
  # dist_fit <- models[-c(7,8)] %>% stan_list(stan_data)
  dist_fit <- models %>% stan_list(stan_data)
  
  if (script_options$`update-output`) {
    new_dist_fit <- dist_fit
  
    dist_fit <- tryCatch({
      load(output_file_name)
      
      list_modify(dist_fit, !!!new_dist_fit)
    })  
  }
} else if (script_options$cv) {
  dist_kfold <- models %>% stan_list(stan_data)
  
  if (script_options$`update-output`) {
    new_dist_kfold <- dist_kfold
  
    dist_kfold <- tryCatch({
      load(output_file_name)
      
      list_modify(dist_kfold, !!!new_dist_kfold)
    })  
  }
}

if (!script_options$`no-save`) {
  if (script_options$fit) {
    save(dist_fit, models, monitored_nosms_data, cluster_analysis_data, grid_dist, extract_sim_level, extract_obs_fit_level, extract_sim_diff, standardize, unstandardize,
         file = output_file_name)
  } else if (script_options$cv) {
    save(dist_kfold, models, monitored_nosms_data, 
         file = output_file_name)
  }
}

