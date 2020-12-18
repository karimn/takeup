data {
  int MIN_COST_MODEL_TYPE_VALUE;
  int MAX_COST_MODEL_TYPE_VALUE;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_KAPPA;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_LINEAR;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_QUADRATIC;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_SEMIPARAM;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_SEMIPARAM_SALIENCE;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_DISCRETE;

  int<lower = 1> num_v_mix;
  
  // Multilevel Configuration 
  int<lower = 0, upper = 1> use_cluster_effects;
  int<lower = 0, upper = 1> use_county_effects;
  int<lower = 0, upper = 1> use_mu_cluster_effects;
  int<lower = 0, upper = 1> use_mu_county_effects;
  int<lower = 0, upper = 1> use_param_dist_cluster_effects; // These are used for parameteric (linear, quadratic) distance cost models only
  int<lower = 0, upper = 1> use_param_dist_county_effects;
  
  int<lower = 0, upper = 1> use_binomial;
  int<lower = 0, upper = 1> use_cost_k_restrictions;
  int<lower = 0, upper = 1> use_private_incentive_restrictions;
  int<lower = 0, upper = 1> use_salience_effect;
  int<lower = 0, upper = 1> use_single_cost_model;
  int<lower = 0, upper = 1> use_name_matched_obs;
  int<lower = 0, upper = 1> use_shifting_v_dist;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> use_cost_model;
  int<lower = 0, upper = 1> suppress_reputation;
  int<lower = 0, upper = 1> suppress_shocks; 
  int<lower = 0, upper = 1> use_u_in_delta;
  int<lower = 0, upper = 1> generate_rep;
  int<lower = 0, upper = 1> generate_sim;
  int<lower = 0, upper = 1> predict_prior;
  
  real<lower = 0> alg_sol_f_tol;
  real<lower = 0> alg_sol_rel_tol;
  int<lower = 0> alg_sol_max_steps;

  int<lower = 1> num_obs; // Actual observations
  int<lower = 1> num_clusters;
  int<lower = 1> num_counties;
  int<lower = 1> num_treatments;
  int<lower = 1> num_discrete_dist;
  
  int<lower = 1, upper = num_clusters> obs_cluster_id[num_obs];
  int<lower = 1, upper = num_counties> obs_county_id[num_obs];
  int<lower = 1, upper = num_counties> cluster_county_id[num_clusters];
 
  int<lower = 0, upper = 1> takeup[num_obs]; // Observed outcome variable
  int<lower = 0, upper = 1> is_name_matched[num_obs];
  vector[num_clusters] cluster_standard_dist; // Standardized distance to treatment
  
  int<lower = 1> cluster_treatment_map[num_treatments * num_discrete_dist, 2];
  int<lower = 1, upper = num_treatments * num_discrete_dist> cluster_assigned_dist_group_treatment[num_clusters]; 
  
  // Reputation
  
  // Semiparametric Cost Model (Splines)
  
  int<lower = 3> num_knots_v;
  matrix[num_clusters, num_knots_v] Z_splines_v; 
  
  real<lower = 0> u_splines_v_sigma_sd;
  
  // K-fold CV 
 
  int<lower = 0, upper =1> cluster_log_lik; 
  int<lower = 0> num_excluded_clusters;
  int<lower = 1, upper = num_clusters> excluded_clusters[num_excluded_clusters];
  
  // Simulation
  
  int<lower = 1> num_grid_obs; // Simulation observations
  int<lower = 1> num_small_grid_obs; // Simulation observations
  vector[num_grid_obs] grid_dist; // Simulation distances
  vector[num_small_grid_obs] small_grid_dist; // Simulation distances
  matrix[num_grid_obs, num_knots_v] Z_grid_v;
  
  // int<lower = 0> num_sim_sm_v;
  // vector[num_sim_sm_v] sim_sm_v;
  
  // Prior hyperparameters
  
  real<lower = 0> mu_rep_sd;
  real<lower = 0> beta_control_sd;
  real<lower = 0> beta_far_effect_sd;
  real<lower = 0> beta_ink_effect_sd;
  real<lower = 0> beta_calendar_effect_sd;
  real<lower = 0> beta_bracelet_effect_sd;
  
  real<lower = 0> structural_beta_county_sd_sd;
  real<lower = 0> structural_beta_cluster_sd_sd;
}

transformed data {
  int num_dist_group_treatments = num_treatments * num_discrete_dist; 
    
  matrix<lower = 0, upper = 1>[num_dist_group_treatments, num_dist_group_treatments] treatment_map_design_matrix = rep_matrix(0, num_dist_group_treatments, num_dist_group_treatments);
  matrix<lower = 0, upper = 1>[num_dist_group_treatments, num_dist_group_treatments] restricted_treatment_map_design_matrix;
  matrix<lower = 0, upper = 1>[num_clusters, num_dist_group_treatments] cluster_treatment_design_matrix = rep_matrix(0, num_clusters, num_dist_group_treatments);
  
  int<lower = 1> cluster_size[num_clusters] = count(num_clusters, obs_cluster_id);
  int<lower = 0> cluster_takeup_count[num_clusters] = count_by_group_test(takeup, obs_cluster_id, { 1 }, 1);
  
  int<lower = 1, upper = num_dist_group_treatments> assigned_treatment[num_obs] = cluster_assigned_dist_group_treatment[obs_cluster_id]; 
    
  int<lower = 0, upper = num_obs> num_name_matched = sum(is_name_matched);
  int<lower = 0, upper = num_obs> num_monitored = num_obs - sum(is_name_matched);
  
  int<lower = 1, upper = num_obs> monitored_obs[num_monitored] = which(is_name_matched, { 0 }, 1);
  int<lower = 1, upper = num_obs> name_matched_obs[num_name_matched] = which(is_name_matched, { 1 }, 1);
 
  int<lower = 0, upper = num_clusters> num_included_clusters = num_clusters - num_excluded_clusters;
  int<lower = 1, upper = num_clusters> included_clusters[num_included_clusters] = num_excluded_clusters > 0 ? which(seq(1, num_clusters, 1), excluded_clusters, 0) : seq(1, num_clusters, 1);
  int<lower = 0, upper = num_obs> num_included_obs = sum(cluster_size[included_clusters]);
  int<lower = 0, upper = num_obs> num_excluded_obs = num_obs - num_included_obs;
  int<lower = 1, upper = num_obs> included_obs[num_included_obs] = num_included_obs != num_obs ? which(obs_cluster_id, included_clusters, 1) : seq(1, num_obs, 1);
  int<lower = 1, upper = num_obs> excluded_obs[num_excluded_obs]; 
  
  int<lower = 0, upper = num_obs> num_included_monitored_obs = num_equals(monitored_obs, included_obs);
  int<lower = 0, upper = num_obs> num_excluded_monitored_obs = num_monitored - num_included_monitored_obs;
  int<lower = 0, upper = num_obs> num_included_name_matched_obs = num_equals(name_matched_obs, included_obs);
  int<lower = 0, upper = num_obs> num_excluded_name_matched_obs = num_name_matched - num_included_name_matched_obs;
  int<lower = 1, upper = num_obs> included_monitored_obs[num_included_monitored_obs] = which(monitored_obs, included_obs, 1);
  int<lower = 1, upper = num_obs> excluded_monitored_obs[num_excluded_monitored_obs]; 
  int<lower = 1, upper = num_obs> included_name_matched_obs[num_included_name_matched_obs] = which(name_matched_obs, included_obs, 1);
  int<lower = 1, upper = num_obs> excluded_name_matched_obs[num_excluded_name_matched_obs];
  
  int<lower = 0, upper = 1> use_dist_salience = in_array(use_cost_model, 
                                                         { COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE, COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE, COST_MODEL_TYPE_SEMIPARAM_SALIENCE });
  int<lower = 0, upper = 1> use_semiparam = in_array(use_cost_model, { COST_MODEL_TYPE_SEMIPARAM, COST_MODEL_TYPE_SEMIPARAM_SALIENCE }); 
  
  int num_treatments_param_kappa = use_cost_model == COST_MODEL_TYPE_PARAM_KAPPA ? num_treatments : 0;
  int num_treatments_param = 0;
  int num_treatments_param_quadratic = 0; 
  int num_treatments_semiparam = 0;
  
  vector[2] sim_grid_mu[num_grid_obs];
  
  int<lower = 1, upper = num_clusters * num_dist_group_treatments> long_cluster_by_treatment_index[num_clusters];
  
  int<lower = 1> num_dist_group_mix = 2;
  
  // long_cluster_by_treatment_index = array_add(array_product(array_subtract(cluster_treatment_map[cluster_assigned_dist_group_treatment, 1], { 1 }),
  long_cluster_by_treatment_index = array_add(array_product(array_subtract(cluster_assigned_dist_group_treatment, { 1 }),
            { num_clusters }),
            seq(1, num_clusters, 1));
  
  for (grid_index in 1:num_grid_obs) {
    sim_grid_mu[grid_index] = rep_vector(0, 2);
  }
  
  if (use_single_cost_model || in_array(use_cost_model, { COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE, 
                                                          COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE, 
                                                          COST_MODEL_TYPE_SEMIPARAM_SALIENCE })) {
    num_treatments_param = 1;
  } else if (in_array(use_cost_model, { COST_MODEL_TYPE_SEMIPARAM, COST_MODEL_TYPE_PARAM_LINEAR, COST_MODEL_TYPE_PARAM_QUADRATIC })) {
    num_treatments_param = num_treatments;
  } 
  
  if ((use_single_cost_model && use_cost_model == COST_MODEL_TYPE_PARAM_QUADRATIC) || use_cost_model == COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE) {
    num_treatments_param_quadratic = 1;
  } else if (in_array(use_cost_model, { COST_MODEL_TYPE_PARAM_QUADRATIC })) {
    num_treatments_param_quadratic = num_treatments;
  } 
  
  if (use_semiparam) {
    if (use_single_cost_model || use_cost_model == COST_MODEL_TYPE_SEMIPARAM_SALIENCE) {
      num_treatments_semiparam = 1;
    } else {
      num_treatments_semiparam = num_treatments;
    }
  }
  
  if (num_excluded_obs > 0) {
    excluded_obs = which(obs_cluster_id, included_clusters, 0);
    excluded_monitored_obs = which(monitored_obs, excluded_obs, 1); 
    excluded_name_matched_obs = which(name_matched_obs, excluded_obs, 1);
  }
  
  treatment_map_design_matrix[, 1] = rep_vector(1, num_dist_group_treatments);
  
  for(treatment_index in 1:num_dist_group_treatments) {
    treatment_map_design_matrix[treatment_index, treatment_index] = 1;
  }
  
  restricted_treatment_map_design_matrix = treatment_map_design_matrix;
  
  if (use_private_incentive_restrictions) {
    // restricted_treatment_map_design_matrix[3,4] = 1; // Calendar effect is > then bracelet
    // restricted_treatment_map_design_matrix[7,8] = 1; // Calendar effect is > then bracelet
  }
 
  cluster_treatment_design_matrix = 
    treatment_map_design_matrix[cluster_assigned_dist_group_treatment]; 
}