
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

extract_sim_level <- function(fit, par, stan_data, quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95)) {
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
                               quant_probs = quant_probs),
           assigned_treatment = factor(assigned_treatment, labels = levels(analysis_data$assigned.treatment))) %>% 
    arrange(assigned_treatment, grid_index) %>%
    mutate(assigned_dist_standard = rep(stan_data$grid_dist, nlevels(assigned_treatment)),
           assigned_dist = unstandardize(assigned_dist_standard, analysis_data$cluster.dist.to.pot))
}

extract_obs_fit_level <- function(fit, par, stan_data, iter_level = c("obs", "cluster", "treatment", "none"), mix = FALSE, quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95)) {
  analysis_data <- monitored_nosms_data
  
  iter_level <- rlang::arg_match(iter_level)
  
  obs_index_col <- if (iter_level == "treatment") "treatment_index" else "obs_index"
  
  fit_data <- fit %>% 
    as.data.frame(par = par) %>% 
    mutate(iter_id = seq_len(n())) %>% 
    gather(indices, iter_est, -iter_id)
  
  if (iter_level == "none") {
    fit_data %>% 
      select(-indices)
  } else {
    fit_data %>% 
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
                                 quant_probs = quant_probs),
             assigned_treatment = switch(iter_level,
                                         cluster = stan_data$cluster_assigned_treatment,
                                         obs = stan_data$assigned_treatment,
                                         treatment = factor(treatment_index, levels = 1:4, labels = levels(stan_data$cluster_assigned_treatment))), 
             assigned_dist_standard = switch(iter_level,
                                             cluster = stan_data$cluster_standard_dist,
                                             obs = stan_data$standard_dist,
                                             treatment = NULL),
             assigned_dist = switch(iter_level,
                                    cluster = unstandardize(assigned_dist_standard, analysis_data$cluster.dist.to.pot),
                                    obs = unstandardize(assigned_dist_standard, analysis_data$cluster.dist.to.pot),
                                    treatment = NULL))
  }
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
                                       probs = quant_probs)) 
}

generate_initializer <- function(num_treatments,
                                 num_clusters,
                                 base_init = function() lst(beta_cluster_sd = abs(rnorm(num_treatments))), 
                                 structural_type = 0,
                                 num_mix = 1,
                                 use_cluster_effects = use_cluster_effects,
                                 use_mu_cluster_effects = use_cluster_effects,
                                 restricted_private_incentive = FALSE,
                                 cost_model_type = NA,
                                 num_knots = NA,
                                 # semiparam = FALSE,
                                 # param_kappa = FALSE,
                                 suppress_reputation = FALSE) {
  base_list <- base_init()
  
  if (structural_type > 0 || !is_empty(base_list)) {
    function() {
      if (structural_type > 0) {
        num_beta_param <- if (restricted_private_incentive) num_treatments - 1 else num_treatments
        # num_beta_param <- num_treatments
        
        salience <- cost_model_type %in% cost_model_types[c("param_linear_salience", "param_quadratic_salience", "semiparam_salience")] 
        param_kappa <- cost_model_type == cost_model_types["param_kappa"]
        linear <- !param_kappa
        quadratic <- cost_model_type %in% cost_model_types[c("param_quadratic", "param_quadratic_salience")] 
        semiparam <- cost_model_type %in% cost_model_types[c("semiparam", "semiparam_salience")] 
        
        init <- lst(
          mu_rep_raw = if (suppress_reputation) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
          mu_rep = if (suppress_reputation) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
          beta_control = rnorm(1, 0, 0.1), 
          beta_ink_effect = rnorm(1, 0, 0.1), 
          beta_calendar_effect = if (restricted_private_incentive) abs(rnorm(1, 0, 0.1)) else rnorm(1, 0, 0.1), 
          beta_bracelet_effect = if (restricted_private_incentive) abs(rnorm(1, 0, 0.1)) else rnorm(1, 0, 0.1), 
          beta_salience = abs(rnorm(1, 0, 0.1)),
          dist_beta_salience = abs(rnorm(1, 0, 0.1)),
          dist_quadratic_beta_salience = abs(rnorm(1, 0, 0.1)),
          # onesided_structural_beta = as.array(- abs(rnorm(num_treatments - num_beta_param, 0, 0.1))),
          dist_cost_k_raw = if (!param_kappa) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
          dist_cost_k = if (!param_kappa) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
          dist_beta_v = if (!linear) array(dim = 0) else if (salience) as.array(abs(rnorm(1, 0, 0.1))) else rnorm(num_treatments, 0, 0.1), 
          dist_quadratic_beta_v = if (!quadratic) array(dim = 0) else if (salience) as.array(abs(rnorm(1, 0, 0.1))) else rnorm(num_treatments, 0, 0.1), 
          u_splines_v_raw = if (!semiparam) {
            array(dim = c(0, num_knots)) 
          } else {
            num_semiparam_treatments <- if (salience) 1 else num_treatments
            array(rnorm(num_semiparam_treatments * num_knots, 0, 0.1), dim = c(num_semiparam_treatments, num_knots))
          },
          u_splines_v = u_splines_v_raw,
          # structural_beta_cluster = if (use_cluster_effects) matrix(rnorm(num_clusters * num_beta_param), nrow = num_clusters, ncol = num_beta_param) else array(dim = 0),
          # structural_beta_cluster_raw = if (use_cluster_effects) matrix(rnorm(num_clusters * num_beta_param), nrow = num_clusters, ncol = num_beta_param) else array(dim = c(0, num_beta_param)),
          # structural_beta_cluster_sd = if (use_cluster_effects) abs(rnorm(num_beta_param, 0, 0.25)) else array(dim = 0),
          structural_beta_cluster = if (use_cluster_effects) matrix(rnorm(num_clusters * num_treatments, 0, 0.1), nrow = num_clusters, ncol = num_treatments) else array(dim = 0),
          structural_beta_cluster_raw = if (use_cluster_effects) matrix(rnorm(num_clusters * num_treatments, 0, 0.1), nrow = num_clusters, ncol = num_treatments) else array(dim = c(0, num_treatments)),
          structural_beta_cluster_sd = if (use_cluster_effects) abs(rnorm(num_treatments, 0, 0.1)) else array(dim = 0),
          structural_cluster_takeup_prob = matrix(rbeta(num_clusters * num_mix, 10, 10), nrow = num_mix),
          # mu_cluster_effects_raw = if (use_cluster_effects && !suppress_reputation) matrix(0, num_clusters, num_treatments) else array(dim = c(0, num_treatments)),
          mu_cluster_effects_raw = if (use_mu_cluster_effects && !suppress_reputation) matrix(rnorm(num_clusters * num_treatments, 0, 0.1), num_clusters, num_treatments) else array(dim = c(0, num_treatments)),
          mu_cluster_effects_sd = if (use_mu_cluster_effects && !suppress_reputation) abs(rnorm(num_treatments, 0, 0.1)) else array(dim = 0),
          cluster_mu_rep = if (use_mu_cluster_effects && !suppress_reputation) matrix(mu_rep, num_clusters, num_treatments, byrow = TRUE) else array(0, num_treatments),
          lambda_v_mix = rep(1 / num_mix, num_mix),
          v_mix_mean = as.array(0.1),
          v_sd = rbeta(1, 8, 1),
          
          # mu_rep_raw = rep(0.05, num_treatments),
          # mu_rep = rep(0.1, num_treatments),
          # structural_beta = rep(0.2, num_beta_param), 
          # onesided_structural_beta = as.array(-0.01), 
          # dist_cost_k_raw = rep(0.05, num_treatments),
          # dist_cost_k = rep(0.05, num_treatments),
          # structural_beta_cluster = if (use_cluster_effects) matrix(rep(0.1, num_clusters * num_treatments), nrow = num_clusters, ncol = num_treatments) else array(dim = 0),
          # structural_beta_cluster_raw = if (use_cluster_effects) matrix(rep(0.2, num_clusters * num_treatments), nrow = num_clusters, ncol = num_treatments) else array(dim = c(0, 4)),
          # structural_beta_cluster_sd = if (use_cluster_effects) rep(0.1, num_treatments) else array(dim = 0),
          # structural_cluster_takeup_prob = rep(0.5, num_clusters) 
          
          # cluster_mu_rep = if (suppress_reputation) array(dim = c(0, num_treatments)) else matrix(mu_rep, num_clusters, num_treatments, byrow = TRUE),
        )
        
        base_list %>% 
          list_modify(!!!init)
                      
      } else {
        return(base_list) 
      }
    }
  } else return(NULL)
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

create_stan_enum <- function(enum_names, prefix) {
  enum_names %>% 
    set_names(seq(length(.)), .) %>% 
    structure(prefix = prefix, class = "stan_enum") 
}

enum2stan <- function(x) UseMethod("enum2stan")
enum2stan.stan_enum <- function(enum) { 
  upper_prefix <- str_to_upper(attr(enum, 'prefix'))
  min_const_name <- str_glue("MIN_{upper_prefix}_VALUE")
  max_const_name <- str_glue("MAX_{upper_prefix}_VALUE")
  values <- str_glue("int<lower = {min_const_name}, upper = {max_const_name}> {upper_prefix}_{str_to_upper(names(enum))};")
  
  print(
    str_glue("int {min_const_name};"),
    str_glue("int {max_const_name};"),
    values)
}

enum2stan_data <- function(x) UseMethod("enum2stan_data")
enum2stan_data.stan_enum <- function(enum) { 
  upper_prefix <- str_to_upper(attr(enum, 'prefix'))
  min_const_name <- str_glue("MIN_{upper_prefix}_VALUE")
  max_const_name <- str_glue("MAX_{upper_prefix}_VALUE")
  
  as.list(enum) %>% 
    set_names(str_c(upper_prefix, "_", str_to_upper(names(enum)))) %>% 
    list_modify(!!min_const_name := min(enum),
                !!max_const_name := max(enum))
}

rep_normal <- function(v, ...) dnorm(v, ...) / ((pnorm(v, ...) * pnorm(v, ..., lower.tail = FALSE)))

generate_v_cutoff_fixedpoint <- function(b, mu) {
  function(v_cutoff) {
    v_cutoff + b + mu * rep_normal(v_cutoff)
  }
}


# Constants ---------------------------------------------------------------

cost_model_types <- create_stan_enum(c("param_kappa", "param_linear", "param_quadratic", "semiparam", 
                                       "param_linear_salience", "param_quadratic_salience", "semiparam_salience"), "cost_model_type") 
