
# Data --------------------------------------------------------------------

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)
# standardize <- as_mapper(~ (. - mean(.)) / sd(.))
# unstandardize <- function(standardized, original) standardized * sd(original) + mean(original)

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

quantilize_est <- function(iter, var, quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95)) {
  iter %>% 
    pull({{ var }}) %>% 
    quantile(probs = quant_probs, names = FALSE) %>% 
    enframe(name = NULL, value = "est") %>% 
    mutate(per = quant_probs)
}

extract_sim_level <- function(fit, par, grid_dist, quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95)) {
  analysis_data <- monitored_nosms_data
  dist_as_data <- is.data.frame(grid_dist)
  
  fit %>% 
    as.data.frame(par = par) %>% 
    mutate(iter_id = seq_len(n())) %>% 
    gather(indices, iter_est, -iter_id) %>% 
    tidyr::extract(indices, c("grid_index", "assigned_treatment"), "(\\d+),(\\d+)", convert = TRUE) %>% 
    group_nest(grid_index, assigned_treatment, .key = "iter_data") %>% 
    mutate(mean_est = map_dbl(iter_data, ~ mean(.$iter_est)),
           quantiles_est = map(iter_data, quantilize_est, 
                               iter_est,
                               # function(iter, quant_probs) 
                               #   quantile(iter$iter_est, probs = quant_probs, names = FALSE) %>% 
                               #   enframe(name = NULL, value = "est") %>% 
                               #   mutate(per = quant_probs),
                               quant_probs = quant_probs),
           assigned_treatment = factor(assigned_treatment, labels = levels(analysis_data$assigned.treatment))) %>% 
    {
      if (dist_as_data) {
        left_join(., grid_dist, by = c("grid_index", "assigned_treatment"))
      } else {
        left_join(., tibble(grid_dist) %>% mutate(grid_index = seq_len(n())), by = "grid_index")  
      }
    } %>% 
    mutate(grid_dist = unstandardize(grid_dist, analysis_data$cluster.dist.to.pot)) 
}

# extract_obs_fit_level <- function(fit, par, stan_data, cluster_level = FALSE, mix = FALSE, quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95)) {
extract_obs_fit_level <- function(fit, par, stan_data, iter_level = c("obs", "cluster", "treatment"), mix = FALSE, quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95)) {
  analysis_data <- monitored_nosms_data
  
  iter_level <- rlang::arg_match(iter_level)
  
  obs_index_col <- if (iter_level == "treatment") "treatment_index" else "obs_index"
  
  fit %>% 
    as.data.frame(par = par) %>% 
    mutate(iter_id = seq_len(n())) %>% 
    gather(indices, iter_est, -iter_id) %>% 
    tidyr::extract(indices, 
                   if (mix) c("mix_index", obs_index_col) else obs_index_col, 
                   if (mix) "(\\d+),(\\d+)" else "(\\d+)", 
                   convert = TRUE) %>% 
    group_by_at(vars(obs_index_col, one_of("mix_index"))) %>% 
    group_nest(.key = "iter_data") %>% 
    ungroup() %>% 
    mutate(mean_est = map_dbl(iter_data, ~ mean(.$iter_est)),
           quantiles_est = map(iter_data, quantilize_est, 
                               iter_est,
                               # function(iter, quant_probs) 
                               #   quantile(iter$iter_est, probs = quant_probs, names = FALSE) %>% 
                               #   enframe(name = NULL, value = "est") %>% 
                               #   mutate(per = quant_probs),
                               quant_probs = quant_probs),
           # assigned_treatment = if (cluster_level) stan_data$cluster_assigned_treatment else stan_data$assigned_treatment,
           assigned_treatment = switch(iter_level,
                                       cluster = stan_data$cluster_assigned_treatment,
                                       obs = stan_data$assigned_treatment,
                                       treatment = factor(treatment_index, levels = 1:4, labels = levels(stan_data$cluster_assigned_treatment))), 
           # assigned_dist = unstandardize(if (cluster_level) stan_data$cluster_standard_dist else stan_data$standard_dist, analysis_data$cluster.dist.to.pot)) 
           assigned_dist = switch(iter_level,
                                  cluster = unstandardize(stan_data$cluster_standard_dist, analysis_data$cluster.dist.to.pot),
                                  obs = unstandardize(stan_data$standard_dist, analysis_data$cluster.dist.to.pot),
                                  treatment = NA))
}

extract_sim_diff <- function(level, quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95)) {
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
           quantiles_est = map(iter_diff, quantilize_est, iter_est, 
                                       # function(iter, quant_probs) 
                                       #   quantile(iter$iter_est, probs = quant_probs, names = FALSE) %>% 
                                       #   enframe(name = NULL, value = "est") %>% 
                                       #   mutate(per = quant_probs),
                                                # interval_size = round(2 * abs(probs - 0.5), 4)),
                                       probs = quant_probs)) 
}

stan_list <- function(models_info, stan_data) {
  get_sample_file_name <- . %>% 
    sprintf("dist_model_%s.csv") %>% 
    file.path("stanfit", .) 
  
  if (script_options$fit) {
    inner_sampler <- function(curr_model, stan_data) {
      control <- lst(adapt_delta = 0.8) 
      
      if (!is_null(curr_model$control)) {
        control %<>% list_modify(!!!curr_model$control)
        curr_model$control <- NULL
      } 
      
      curr_stan_data <- stan_data %>%
        list_modify(!!!curr_model) %>%
        map_if(is.factor, as.integer)
        
      dist_model <- stan_model(file.path("stan_models", curr_model$model_file))
      
      fit <- sampling(
        dist_model,
        iter = if (script_options$`force-iter`) iter else (curr_stan_data$iter %||% iter),
        thin = thin_by,
        chains = chains,
        control = control,
        save_warmup = FALSE,
        refresh = 100,
        init = curr_model$init %||% "random",
        data = curr_stan_data)
      
      return(fit)
    }
    
    if (script_options$sequential) {
      models_info %>% 
        map(inner_sampler,
            stan_data = stan_data)
      
    } else {
      models_info %>% 
        pbmclapply(inner_sampler,
                   stan_data = stan_data,
                   ignore.interactive = TRUE,
                   mc.silent = TRUE,
                   mc.cores = 3)
    }
  } else if (script_options$cv) {
    models_info %>% 
      imap(function(curr_model, model_name, stan_data) {
        kfold_groups <- kfold_split_stratified(K = folds, x = stan_data$cluster_assigned_dist_group_treatment)
        
        dist_model <- stan_model(file.path("stan_models", curr_model$model_file))
        
        log_lik_list <- map(seq(folds), ~ which(kfold_groups == .)) %>% 
          pbmclapply(function(excluded_clusters, model_name, dist_model, stan_data) {
            curr_stan_data <- stan_data %>%
              list_modify(!!!curr_model,
                          excluded_clusters = excluded_clusters,
                          num_excluded_clusters = length(excluded_clusters)) %>%
              map_if(is.factor, as.integer)
            
            curr_iter <- if (script_options$`force-iter`) iter else (curr_stan_data$iter %||% iter)
            curr_chains <- chains
            
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
      stan_data = stan_data) 
  }
}
