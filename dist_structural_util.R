
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

# Model Classes -----------------------------------------------------------

# setClass("TakeUpModel", 
#          contains = "stanmodel",
#          slots = c("pars" = "character"))
# 
# setClass("ReducedFormTakeUpModel",  
#          contains = "TakeUpModel")
# 
# setClass("StructuralTakeUpModel",  
#          contains = "TakeUpModel")
# 
# create_initializer <- function(model_file, pars) {
#   compiled_model <- stan_model(model_file)
#   
#   function(.Object, ...) {
#     .Object <- callNextMethod()
#     
#     .Object@pars <- pars 
#     return(.Object)
#   }
# }
# 
# setMethod("initialize", "ReducedFormTakeUpModel", 
#           create_initializer(file.path("stan_models", "takeup_reduced.stan"), 
#                              pars = c("structural_cluster_benefit_cost", "structural_cluster_obs_v", "structural_cluster_takeup_prob", "beta", "cluster_cf_benefit_cost", "cluster_rep_benefit_cost")))
# 
# setMethod("initialize", "StructuralTakeUpModel", 
#           create_initializer(file.path("stan_models", "dist_struct_fixedpoint.stan"), 
#                              pars = c("total_error_sd", "cluster_dist_cost", "structural_cluster_benefit_cost", "structural_cluster_obs_v", "structural_cluster_takeup_prob",
#                                       "beta", "dist_beta_v", "mu_rep", "cluster_cf_benefit_cost", "mu_cluster_effects_raw", "mu_cluster_effects_sd", "cluster_mu_rep", 
#                                       "cluster_rep_benefit_cost", "sim_benefit_cost",
#                                       "group_dist_mean", "group_dist_sd", "group_dist_mix",
#                                       "dist_beta_county_raw", "dist_beta_county_sd")))

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

extract_sim_level <- function(fit, par, stan_data, quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95)) {
  analysis_data <- stan_data$analysis_data
  dist_as_data <- is.data.frame(stan_data$grid_dist)
  
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

extract_rep_level <- function(fit, par, stan_data, quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95), thin = 1, dewormed_var = dewormed) {
  analysis_data <- stan_data$analysis_data
  
  fit_data <- fit %>% 
    as.array(par = par) %>% {
      if (thin > 1) magrittr::extract(., (seq_len(nrow(.)) %% thin) == 0,,) else .
    } %>% 
    plyr::adply(3, function(cell) tibble(iter_data = list(cell)) %>% mutate(ess_bulk = if (thin > 1) ess_bulk(cell), 
                                                                            ess_tail = if (thin > 1) ess_tail(cell),
                                                                            rhat = if (thin > 1) Rhat(cell),)) %>% 
    tidyr::extract(parameters, "cluster_id", "(\\d+)", convert = TRUE) %>%  
    mutate(iter_data = map(iter_data, ~ tibble(iter_est = c(.), iter_id = seq(nrow(.) * ncol(.)))),
           assigned_dist_standard = stan_data$cluster_standard_dist[cluster_id],
           assigned_dist = unstandardize(assigned_dist_standard, analysis_data$cluster.dist.to.pot)) %>% 
    # left_join(stan_data$analysis_data %>% count(cluster_id, assigned_treatment = assigned.treatment, name = "cluster_size"), by = "cluster_id") 
    left_join(stan_data$analysis_data %>%
                select(cluster_id, assigned_treatment = assigned.treatment, dewormed) %>%
                add_count(cluster_id, name = "cluster_size") %>% 
                group_by(cluster_id, assigned_treatment, cluster_size) %>%
                summarize(obs_num_takeup = sum({{ dewormed_var }})) %>%
                ungroup(),
              by = c("cluster_id")) 
  
  fit_data %<>% 
    as_tibble()
}

extract_obs_cf <- function(fit, par, stan_data, iter_level = c("obs", "cluster"), quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95), thin = 1, dewormed_var = dewormed, always_diagnose = TRUE) {
  analysis_data <- stan_data$analysis_data
  
  iter_level <- rlang::arg_match(iter_level)
  
  obs_index_col <- switch(iter_level, obs = "obs_index", cluster = "cluster_id")
  
  cluster_treatment_map <- stan_data$cluster_treatment_map %>% 
    mutate(treatment_index = seq(n()))
  
  is_stanfit <- is(fit, "stanfit")
  
  fit %>% {
    if (is_stanfit) {
      as.array(par = par) 
    } else {
      .$draws(par)
    }
  } %>% {
      if (thin > 1) magrittr::extract(., (seq_len(nrow(.)) %% thin) == 0,,) else .
    } %>%  
    plyr::adply(3, function(cell) tibble(iter_data = list(cell)) %>% mutate(ess_bulk = if (is_stanfit && (thin > 1 || always_diagnose)) ess_bulk(cell), 
                                                                            ess_tail = if (is_stanfit && (thin > 1 || always_diagnose)) ess_tail(cell),
                                                                            rhat = if (is_stanfit && (thin > 1 || always_diagnose)) Rhat(cell),)) %>% {
      if (is_stanfit) rename(., variable = parameters) else .
    } %>% 
    tidyr::extract(variable, c("treatment_index", obs_index_col), "(\\d+),(\\d+)", convert = TRUE) %>%  
    mutate(iter_data = map(iter_data, ~ tibble(iter_est = c(.), iter_id = seq(nrow(.) * ncol(.)))),
           cluster_id = switch(iter_level,
                               cluster = cluster_id, 
                               obs = stan_data$obs_cluster_id[!!sym(obs_index_col)]),
           treatment_index_obs = stan_data$cluster_assigned_dist_group_treatment[cluster_id],
           assigned_dist_standard_obs = stan_data$cluster_standard_dist[cluster_id],
           assigned_dist_obs = unstandardize(assigned_dist_standard_obs, analysis_data$cluster.dist.to.pot)) %>% 
    left_join(cluster_treatment_map, by = "treatment_index") %>% 
    left_join(cluster_treatment_map, by = c("treatment_index_obs" = "treatment_index"), suffix = c("", "_obs")) %>% 
    left_join(stan_data$analysis_data %>% count(cluster_id, name = "cluster_size"), by = "cluster_id") %>% 
    left_join(stan_data$analysis_data %>%
                select(cluster_id, assigned_treatment = assigned.treatment, assigned_dist_group = dist.pot.group, {{ dewormed_var }}) %>%
                group_by(cluster_id, assigned_treatment, assigned_dist_group) %>%
                summarize(obs_num_takeup = sum({{ dewormed_var }})) %>%
                ungroup(),
              by = c("cluster_id", "assigned_treatment", "assigned_dist_group")) %>% 
    as_tibble()
}

extract_obs_fit_level <- function(fit, par, stan_data, iter_level = c("obs", "cluster", "none"), by_treatment = FALSE, mix = FALSE, summarize_est = TRUE, quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95)) {
  analysis_data <- stan_data$analysis_data
  
  iter_level <- rlang::arg_match(iter_level)
  
  # obs_index_col <- if (iter_level == "treatment") "treatment_index" else "obs_index"
  obs_index_col <- if (iter_level != "none") "obs_index"
  
  if (by_treatment) obs_index_col %<>% c("treatment_index") 
  if (mix) obs_index_col %<>% c("mix_index", .) 
  
  is_stanfit <- is(fit, "stanfit")
  
  fit_data <- tryCatch({
    fit_data <- if (is_stanfit) {
      fit %>% 
        as.array(par = par) %>% 
        plyr::adply(3, function(cell) tibble(iter_data = list(cell)) %>% mutate(ess_bulk = ess_bulk(cell), 
                                                                                ess_tail = ess_tail(cell),
                                                                                rhat = Rhat(cell),)) 
    } else {
      fit_data <- fit$draws(par) 
      
      if (length(dim(fit_data)) == 3) {
        plyr::adply(fit_data, 3, function(cell) tibble(iter_data = list(cell))) 
      } else {
        tibble(iter_data = list(fit_data)) %>% 
          mutate(variable = par)
      }
    }
    
    fit_data %>% 
      mutate(iter_data = map(iter_data, ~ tibble(iter_est = c(.), iter_id = seq(nrow(.) * ncol(.))))) %>% 
      as_tibble()
  }, error = function(err) NULL)
  
  if (is_null(fit_data)) return(NULL) 
  
  if (iter_level == "none" && !by_treatment) {
    if (is_stanfit) {
      fit_data %>% 
        select(-parameters) 
    } else {
      fit_data %>% 
        select(-variable)
    }
  } else {
    extracted <- fit_data %>% {
      if (is_stanfit) rename(., variable = parameters) else .
    } %>% 
      tidyr::extract(variable, 
                     obs_index_col, 
                     str_c(rep_along(obs_index_col, r"{(\d+)}"), collapse = ","),
                     convert = TRUE)
    
    if (summarize_est) {
      extracted %<>% 
        mutate(mean_est = map_dbl(iter_data, ~ mean(.$iter_est)),
               quantiles_est = map(iter_data, quantilize_est, 
                                   iter_est,
                                   quant_probs = quant_probs))
    }
    
    if (by_treatment) {
      extracted %<>% 
        mutate(
          assigned_treatment = factor(treatment_index, levels = 1:4, labels = levels(stan_data$cluster_assigned_treatment)), 
          assigned_dist_standard = NULL,
          assigned_dist = NULL 
        )
        
    } else {
      extracted %<>% 
        mutate(
          assigned_treatment = switch(iter_level,
                                      cluster = stan_data$cluster_assigned_treatment,
                                      obs = stan_data$assigned_treatment),
          assigned_dist_standard = switch(iter_level,
                                          cluster = stan_data$cluster_standard_dist,
                                          obs = stan_data$standard_dist),
          assigned_dist = switch(iter_level,
                                 cluster = unstandardize(assigned_dist_standard, analysis_data$cluster.dist.to.pot),
                                 obs = unstandardize(assigned_dist_standard, analysis_data$cluster.dist.to.pot)),
        )
    }
    
    return(extracted)
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
                                 num_counties,
                                 base_init = function() lst(beta_cluster_sd = abs(rnorm(num_treatments))), 
                                 structural_type = 0,
                                 num_mix = 1,
                                 use_cluster_effects = use_cluster_effects,
                                 use_county_effects = use_cluster_effects,
                                 use_param_dist_cluster_effects = use_cluster_effects,
                                 use_mu_cluster_effects = use_cluster_effects,
                                 use_mu_county_effects = FALSE,
                                 use_single_cost_model = FALSE,
                                 restricted_private_incentive = FALSE,
                                 cost_model_type = NA,
                                 num_knots = NA,
                                 name_matched = FALSE,
                                 suppress_reputation = FALSE) {
  base_list <- base_init()
  
  if (structural_type > 0 || !is_empty(base_list)) {
    function() {
      if (structural_type > 0) {
        num_beta_param <- if (restricted_private_incentive) num_treatments - 1 else num_treatments
        # num_beta_param <- num_treatments
        
        salience <- cost_model_type %in% cost_model_types[c("param_linear_salience", "param_quadratic_salience", "semiparam_salience")] 
        param_kappa <- cost_model_type == cost_model_types["param_kappa"]
        discrete_cost <- cost_model_type == cost_model_types["discrete"]
        linear <- !param_kappa && !discrete_cost
        quadratic <- cost_model_type %in% cost_model_types[c("param_quadratic", "param_quadratic_salience")] 
        semiparam <- cost_model_type %in% cost_model_types[c("semiparam", "semiparam_salience")] 
        
        num_discrete_dist <- if (cost_model_type %in% cost_model_types["discrete"]) 2 else 1
        num_incentive_treatments <- num_treatments
        num_treatments <- if (cost_model_type %in% cost_model_types["discrete"]) num_treatments * num_discrete_dist else num_treatments
        
        init <- lst(
          mu_rep_raw = if (suppress_reputation && !salience) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
          mu_rep = if (suppress_reputation && !salience) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
          v_mu = rnorm(1, 0, 0.1),
          beta_control = rnorm(num_discrete_dist, 0, 0.1) %>% { if (num_discrete_dist > 1) as.array(.) else . }, 
          beta_ink_effect = rnorm(num_discrete_dist, 0, 0.1) %>% { if (num_discrete_dist > 1) as.array(.) else . },
          beta_calendar_effect = { 
            if (restricted_private_incentive) abs(rnorm(num_discrete_dist, 0, 0.1)) else rnorm(num_discrete_dist, 0, 0.1) 
          } %>% { if (num_discrete_dist > 1) as.array(.) else . },
          beta_bracelet_effect = { 
            if (restricted_private_incentive) abs(rnorm(num_discrete_dist, 0, 0.1)) else rnorm(num_discrete_dist, 0, 0.1)
          } %>% { if (num_discrete_dist > 1) as.array(.) else . },
          beta_salience = abs(rnorm(1, 0, 0.1)),
          dist_beta_salience = abs(rnorm(1, 0, 0.1)),
          dist_quadratic_beta_salience = abs(rnorm(1, 0, 0.1)),
          beta_nm_effect = if (name_matched) rnorm(num_treatments, 0, 0.1) else array(dim = 0),
          dist_cost_k_raw = if (!param_kappa) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
          dist_cost_k = if (!param_kappa) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
          dist_beta_v = if (!linear) array(dim = 0) else if (salience || use_single_cost_model) as.array(abs(rnorm(1, 0, 0.1))) else rnorm(num_treatments, 0, 0.1), 
          dist_quadratic_beta_v = if (!quadratic) array(dim = 0) else if (salience || use_single_cost_model) as.array(abs(rnorm(1, 0, 0.1))) else abs(rnorm(num_treatments, 0, 0.1)), 
          u_splines_v_raw = if (!semiparam) {
            array(dim = c(0, num_knots)) 
          } else {
            num_semiparam_treatments <- if (salience || use_single_cost_model) 1 else num_treatments
            array(rnorm(num_semiparam_treatments * num_knots, 0, 0.1), dim = c(num_semiparam_treatments, num_knots))
          },
          u_splines_v = u_splines_v_raw,
         
          u_sd = abs(rnorm(1, 0, 0.5)),
          # u_sd = 0.01,
          # test_w = 0.461662,
          
          structural_beta_cluster = if (use_cluster_effects) matrix(rnorm(num_clusters * num_treatments, 0, 0.1), nrow = num_clusters, ncol = num_treatments) else array(dim = 0),
          structural_beta_cluster_raw = if (use_cluster_effects) matrix(rnorm(num_clusters * num_treatments, 0, 0.1), nrow = num_clusters, ncol = num_treatments) else array(dim = c(0, num_treatments)),
          structural_beta_cluster_sd = if (use_cluster_effects) abs(rnorm(num_treatments, 0, 0.1)) else array(dim = 0),
          
          structural_beta_county = if (use_county_effects) matrix(rnorm(num_counties * num_treatments, 0, 0.1), nrow = num_counties, ncol = num_treatments) else array(dim = 0),
          structural_beta_county_raw = if (use_county_effects) matrix(rnorm(num_counties * num_treatments, 0, 0.1), nrow = num_counties, ncol = num_treatments) else array(dim = c(0, num_treatments)),
          structural_beta_county_sd = if (use_county_effects) abs(rnorm(num_treatments, 0, 0.1)) else array(dim = 0),
          
          beta_nm_effect_cluster = if (use_cluster_effects && name_matched) matrix(rnorm(num_clusters * num_treatments, 0, 0.1), nrow = num_clusters, ncol = num_treatments) else array(dim = c(0, num_treatments)),
          beta_nm_effect_cluster_raw = if (use_cluster_effects && name_matched) matrix(rnorm(num_clusters * num_treatments, 0, 0.1), nrow = num_clusters, ncol = num_treatments) else array(dim = c(0, num_treatments)),
          beta_nm_effect_cluster_sd = if (use_cluster_effects && name_matched) abs(rnorm(num_treatments, 0, 0.1)) else array(dim = 0),
          
          structural_cluster_takeup_prob = matrix(rbeta(num_clusters * num_mix, 10, 10), nrow = num_mix),
          
          dist_beta_cluster_raw = if (discrete_cost) {
            array(0, dim = c(0, 0))
          } else if (!use_param_dist_cluster_effects) { 
            if (salience || use_single_cost_model) array(dim = c(0, 1)) else array(dim = c(0, num_treatments)) 
          } else {
            if (salience || use_single_cost_model) matrix(rnorm(num_clusters, 0, 0.1), nrow = num_clusters, ncol = 1)
            else matrix(rnorm(num_clusters * num_treatments, 0, 0.1), nrow = num_clusters, ncol = num_treatments)
          },
          dist_beta_cluster_sd = if (!use_param_dist_cluster_effects) array(dim = 0) 
            else if (salience || use_single_cost_model) as.array(abs(rnorm(1, 0, 0.1)))
            else abs(rnorm(num_treatments, 0, 0.1)),
          # mu_cluster_effects_raw = if (use_mu_cluster_effects && !suppress_reputation) matrix(rnorm(num_clusters * num_treatments, 0, 0.1), num_clusters, num_treatments) else array(dim = c(0, num_incentive_treatments)),
          # mu_cluster_effects_sd = if (use_mu_cluster_effects && !suppress_reputation) abs(rnorm(num_treatments, 0, 0.1)) else array(dim = 0),
          # mu_county_effects_raw = if (use_mu_county_effects && !suppress_reputation) matrix(rnorm(num_counties * num_treatments, 0, 0.1), num_counties, num_treatments) else array(dim = c(0, num_incentive_treatments)),
          # mu_county_effects_sd = if (use_mu_county_effects && !suppress_reputation) abs(rnorm(num_treatments, 0, 0.1)) else array(dim = 0),
          mu_cluster_effects_raw = if (!use_mu_cluster_effects || (suppress_reputation && !salience)) array(dim = c(0, num_incentive_treatments)) else matrix(rnorm(num_clusters * num_treatments, 0, 0.1), num_clusters, num_treatments),
          mu_cluster_effects_sd = if (!use_mu_cluster_effects || (suppress_reputation && !salience)) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)) ,
          mu_county_effects_raw = if (!use_mu_cluster_effects || (suppress_reputation && !salience)) array(dim = c(0, num_incentive_treatments)) else matrix(rnorm(num_counties * num_treatments, 0, 0.1), num_counties, num_treatments) ,
          mu_county_effects_sd = if (!use_mu_cluster_effects || (suppress_reputation && !salience)) array(dim = 0) else abs(rnorm(num_treatments, 0, 0.1)),
          cluster_mu_rep = if (!use_mu_cluster_effects || (suppress_reputation && !salience)) array(0, num_treatments) else matrix(mu_rep, num_clusters, num_treatments, byrow = TRUE),
          lambda_v_mix = rep(1 / num_mix, num_mix),
          v_mix_mean = as.array(0.1),
          v_sd = rbeta(1, 8, 1),
          
          group_dist_mix = MCMCpack::rdirichlet(2, rep(10, 2))
        )
        
        base_list %>% 
          list_modify(!!!init)
                      
      } else {
        return(base_list) 
      }
    }
  } else return(NULL)
}

stan_list <- function(models_info, stan_data, use_cmdstanr = FALSE, include_paths = ".") {
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
      
      
      iter <- if (script_options$force_iter) iter else (curr_stan_data$iter %||% iter)
        
      fit <- if (use_cmdstanr) {
        dist_model <- cmdstan_model(file.path("stan_models", curr_model$model_file), include_paths = include_paths)
        
        dist_model$sample(
          data = curr_stan_data %>% 
            discard(~ is_function(.x) | is.character(.)) %>% 
            list_modify(analysis_data = NULL),
          chains = chains,
          parallel_chains = chains, 
          iter_warmup = iter %/% 2, 
          iter_sampling = iter %/% 2, 
          save_warmup = FALSE,
          thin = (curr_stan_data$thin %||% 1),
          init = curr_model$init,
          adapt_delta = control$adapt_delta,
          max_treedepth = control$max_treedepth
        )
      } else { 
        stan_model(file.path("stan_models", curr_model$model_file)) %>% 
          sampling(
            iter = iter, 
            thin = (curr_stan_data$thin %||% 1),
            chains = chains,
            control = control,
            save_warmup = FALSE,
            pars = curr_model$pars,
            init = curr_model$init %||% "random",
            data = curr_stan_data
          )
      }
      
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
      imap(function(curr_model, model_name, stan_data, use_cmdstanr) {
        kfold_groups <- kfold_split_stratified(K = folds, x = stan_data$cluster_assigned_dist_group_treatment)
        
        dist_model <- if (use_cmdstanr) {
          cmdstan_model(file.path("stan_models", curr_model$model_file), include_paths = include_paths)
        } else {
          stan_model(file.path("stan_models", curr_model$model_file))
        }
        
        log_lik_list <- map(seq(folds), ~ which(kfold_groups == .)) %>% 
          pbmclapply(function(excluded_clusters, model_name, dist_model, stan_data, use_cmdstanr) {
            curr_stan_data <- stan_data %>%
              list_modify(!!!curr_model,
                          excluded_clusters = excluded_clusters,
                          num_excluded_clusters = length(excluded_clusters)) %>%
              map_if(is.factor, as.integer)
            
            curr_iter <- if (script_options$force_iter) iter else (curr_stan_data$iter %||% iter)
            curr_chains <- chains
            
            if (use_cmdstanr) {
              fit <- dist_model$sample(
                data = curr_stan_data,
                iter_warmup = curr_iter %/% 2,
                iter_sampling = curr_iter %/% 2,
                save_warmup = FALSE,
                thin = thin_by,
                chains = curr_chains,
                parallel_chains = chains,
                adapt_delta = 0.9,,
                init = curr_model$init
              )
              
              fit$draws("log_lik_heldout")
            } else {
              sampling(dist_model,
                       iter = curr_iter,
                       thin = thin_by,
                       chains = curr_chains,
                       control = lst(adapt_delta = 0.9),
                       save_warmup = FALSE,
                       init = curr_model$init %||% "random",
                       pars = "log_lik_heldout", 
                       data = curr_stan_data) %>% 
                extract_log_lik("log_lik_heldout") 
            }
          },
          dist_model = dist_model,
          stan_data = stan_data,
          use_cmdstanr = use_cmdstanr,
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
      stan_data = stan_data, use_cmdstanr = use_cmdstanr) 
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

rep_normal_1std <- function(v, ...) {
  phi_v <- dnorm(v, ...)
  Phi_1mPhi_v <- ((pnorm(v, ...) * pnorm(v, ..., lower.tail = FALSE)))
  
  (phi_v / Phi_1mPhi_v^2) * (-(v * Phi_1mPhi_v) + (phi_v * (2 * pnorm(v, ...) - 1))) 
}

social_multiplier <- function(v, mu) {
  - 1 / (1 + mu * rep_normal_1std(v))
}

expect_y_partial_bbar <- function(v, mu, sigma) {
  - dnorm(v, sd = sigma) * social_multiplier(v, mu) 
}

expect_y_partial_c <- function(v, mu, sigma) {
  dnorm(v, sd = sigma) * social_multiplier(v, mu) 
}

expect_y_partial_d <- function(v, mu, sigma, linear_dist_cost) {
  expect_y_partial_c(v, mu, sigma) %>% 
    multiply_by(map_dbl(linear_dist_cost, ~ weighted.mean(.x$iter_est, .x$cluster_size))) 
}

# Constants ---------------------------------------------------------------

cost_model_types <- create_stan_enum(c("param_kappa", "param_linear", "param_quadratic", "semiparam", 
                                       "param_linear_salience", "param_quadratic_salience", "semiparam_salience",
                                       "discrete"), "cost_model_type") 
