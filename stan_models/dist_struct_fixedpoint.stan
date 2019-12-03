functions {
#include /../multilvlr/util.stan
  
  real reputational_returns_normal(real v, vector lambda, vector mix_mean, vector mix_sd) {
    int num_mix = num_elements(lambda);
    real rep = 0; 
    
    for (mix_index in 1:num_mix) {
      real mix_Phi_v = Phi((v - mix_mean[mix_index]) / mix_sd[mix_index]);
      
      rep += lambda[mix_index] * exp(normal_lpdf(v | mix_mean[mix_index], mix_sd[mix_index])) / (mix_Phi_v * (1 - mix_Phi_v));
    }
    
    return rep;
  }
  
  vector param_kappa_dist_cost(vector dist, vector k) {
    if (num_elements(k) == 1) {
      return (k[1] * square(dist)) / 2;
    } else {
      return (k .* square(dist)) / 2;
    }
  }
  
  vector param_dist_cost(vector dist, vector linear_dist_cost, vector quadratic_dist_cost, matrix u_splines, matrix Z_splines) {
    int num_cost = num_elements(dist);
    vector[num_cost] cost = rows_dot_product(u_splines, Z_splines);
    
    if (num_elements(linear_dist_cost) == 1) {
      cost += dist * linear_dist_cost[1];
    } else {
      cost += dist .* linear_dist_cost;
    }
    
    if (num_elements(quadratic_dist_cost) == 1) {
      cost += square(dist) * quadratic_dist_cost[1];
    } else {
      cost += square(dist) .* quadratic_dist_cost;
    }
    
    return cost;
  }
  
  vector v_fixedpoint_solution_normal(vector model_param, vector theta, real[] x_r, int[] x_i) {
    real v_cutoff = model_param[1];
    
    real benefit_cost = theta[1];
    real mu = theta[2];
    
    int num_v_mix = x_i[1];
    
    vector[num_v_mix] lambda = theta[3:(3 + num_v_mix - 1)];
    vector[num_v_mix] mix_mean = theta[(3 + num_v_mix):(3 + 2 * num_v_mix - 1)];
    vector[num_v_mix] mix_sd = theta[(3 + 2 * num_v_mix):(3 + 3 * num_v_mix - 1)];
    
    return [ v_cutoff + benefit_cost + mu * reputational_returns_normal(v_cutoff, lambda, mix_mean, mix_sd) ]';
  }
  
  real mixed_binomial_lpmf(int[] outcomes, vector lambda, int[] N, matrix prob) {
    int num_obs = num_elements(outcomes);
    int num_mix = num_elements(lambda);
    real logp = 0;
    
    for (obs_index in 1:num_obs) {
      vector[num_mix] curr_logp;
      
      for (mix_index in 1:num_mix) {
        curr_logp[mix_index] = log(lambda[mix_index]) + binomial_lpmf(outcomes[obs_index] | N[obs_index], prob[mix_index, obs_index]); 
      }
    
      logp += log_sum_exp(curr_logp);  
    }
    
    return logp;
  }
  
  real mixed_bernoulli_lpmf(int[] outcomes, vector lambda, matrix prob) {
    int num_obs = num_elements(outcomes);
    int num_mix = num_elements(lambda);
    real logp = 0;
    
    for (obs_index in 1:num_obs) {
      vector[num_mix] curr_logp;
      
      for (mix_index in 1:num_mix) {
        curr_logp[mix_index] = log(lambda[mix_index]) + bernoulli_lpmf(outcomes[obs_index] | prob[mix_index, obs_index]); 
      }
   
      logp += log_sum_exp(curr_logp);  
    }
    
    return logp;
  }
  
  vector prepare_solver_theta(real benefit_cost, real mu_rep, real v_mu, vector lambda, vector mix_mean, real v_sd, vector mix_sd) {
    int num_v_mix = num_elements(lambda);
    vector[2 + 3 * num_v_mix] solver_theta;
    
    solver_theta[1:2] = [ benefit_cost, mu_rep ]';
    solver_theta[3:(3 + num_v_mix - 1)] = lambda;
    solver_theta[(3 + num_v_mix):(3 + num_v_mix + num_v_mix - 1)] = append_row(v_mu, mix_mean);
    solver_theta[(3 + num_v_mix + num_v_mix):(3 + 3 * num_v_mix - 1)] = append_row(v_sd, mix_sd);
 
    return solver_theta; 
  }
  
  int[] prepare_cluster_assigned_dist_group_treatment(int[] cluster_assigned_treatment, int[] cluster_assigned_dist_group) {
    int num_clusters = num_elements(cluster_assigned_treatment);
    int num_treatments = max(cluster_assigned_treatment);
    int num_dist = max(cluster_assigned_dist_group);
    
    int treatment_id[num_clusters];
   
    int treatment_id_map[num_treatments, num_dist];
    
    for (dist_index in 1:num_dist) {
      for (treatment_index in 1:num_treatments) {
        treatment_id_map[treatment_index, dist_index] = treatment_index + dist_index - 1;
      }
    }
    
    for (cluster_index in 1:num_clusters) {
      treatment_id[cluster_index] = treatment_id_map[cluster_assigned_treatment[cluster_index], cluster_assigned_dist_group[cluster_index]];
    } 
    
    return treatment_id;
  }
  
  real normal_lb_rng(real mu, real sigma, real lb) {
    real p = normal_cdf(lb, mu, sigma);  // cdf for bounds
    real u = uniform_rng(p, 1);
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf for value
  }
}

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
  int<lower = 0, upper = 1> generate_rep;

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
  
  int<lower = 1, upper = num_treatments> cluster_assigned_treatment[num_clusters]; // Actual assigned treatments 
  int<lower = 1, upper = num_discrete_dist> cluster_assigned_dist_group[num_clusters];
  int<lower = 1, upper = num_treatments * num_discrete_dist> cluster_assigned_dist_group_treatment[num_clusters]; 
  
  // Reputation
  
  real<lower = 0> mu_rep_sd;
  
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
}

transformed data {
  int num_actual_treatments = use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_treatments * num_discrete_dist : num_treatments; 
    
  matrix<lower = 0, upper = 1>[num_actual_treatments, num_actual_treatments] treatment_map_design_matrix = rep_matrix(0, num_actual_treatments, num_actual_treatments);
  matrix<lower = 0, upper = 1>[num_actual_treatments, num_actual_treatments] restricted_treatment_map_design_matrix;
  matrix<lower = 0, upper = 1>[num_clusters, num_actual_treatments] cluster_treatment_design_matrix = rep_matrix(0, num_clusters, num_actual_treatments);
  
  int<lower = 1> cluster_size[num_clusters] = count(num_clusters, obs_cluster_id);
  int<lower = 0> cluster_takeup_count[num_clusters] = count_by_group_test(takeup, obs_cluster_id, { 1 }, 1);
  
  int<lower = 1, upper = num_actual_treatments> assigned_treatment[num_obs] = 
    use_cost_model == COST_MODEL_TYPE_DISCRETE ? cluster_assigned_dist_group_treatment[obs_cluster_id] : cluster_assigned_treatment[obs_cluster_id]; 
    
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
  int num_treatments_name_matched = use_name_matched_obs ? num_actual_treatments : 0;
  
  int<lower = 0> num_shocks = suppress_shocks ? 0 : (use_cost_model == COST_MODEL_TYPE_DISCRETE ? 2 : 3); // U_v, U_b, U_r: shocks to altruism, private benefits and reputational benefits
  
  vector[num_shocks] sim_grid_mu[num_grid_obs];
  
  int<lower = 0, upper = num_obs> num_treatment_obs[num_actual_treatments] = count(num_actual_treatments, assigned_treatment);
  
  int<lower = 0, upper = num_obs> num_treatment_cf[num_actual_treatments] = array_subtract(rep_array(num_obs, num_actual_treatments), num_treatment_obs);
  
  int num_cf_obs = num_obs * (num_actual_treatments - 1);
  
  int<lower = 1, upper = num_clusters * num_actual_treatments> long_cluster_by_treatment_index[num_clusters];
  
  if (use_cost_model == COST_MODEL_TYPE_DISCRETE) {
    long_cluster_by_treatment_index = array_add(array_product(array_subtract(cluster_assigned_dist_group_treatment, { 1 }), 
              { num_clusters }), 
              seq(1, num_clusters, 1));
  } else {
    long_cluster_by_treatment_index = array_add(array_product(array_subtract(cluster_assigned_treatment, { 1 }), 
              { num_clusters }), 
              seq(1, num_clusters, 1));
  }
  
  for (grid_index in 1:num_grid_obs) {
    sim_grid_mu[grid_index] = rep_vector(0, num_shocks);
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
  
  treatment_map_design_matrix[, 1] = rep_vector(1, num_actual_treatments);
  
  
  {
    int treatment_cf_size_pos = 1;
    int cluster_treatment_cf_size_pos = 1;
    
    for(treatment_index in 1:num_actual_treatments) {
      int treatment_cf_size_end = treatment_cf_size_pos + num_treatment_cf[treatment_index] - 1;
      
      treatment_map_design_matrix[treatment_index, treatment_index] = 1;
      
      treatment_cf_size_pos = treatment_cf_size_end + 1;
    }
  }
  
  
  restricted_treatment_map_design_matrix = treatment_map_design_matrix;
  
  if (use_private_incentive_restrictions) {
    restricted_treatment_map_design_matrix[3,4] = 1; // Calendar effect is > then bracelet
    
    if (use_cost_model == COST_MODEL_TYPE_DISCRETE) {
      restricted_treatment_map_design_matrix[7,8] = 1; // Calendar effect is > then bracelet
    }
  }
 
  if (use_cost_model == COST_MODEL_TYPE_DISCRETE) {
    cluster_treatment_design_matrix = 
      treatment_map_design_matrix[cluster_assigned_dist_group_treatment]; 
  } else {
    cluster_treatment_design_matrix = 
      treatment_map_design_matrix[cluster_assigned_treatment]; 
  } 
}

parameters {
  // Levels: control ink calendar bracelet
  vector[use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_discrete_dist : 1] beta_control;
  vector[use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_discrete_dist : 1] beta_ink_effect;
  vector<lower = (use_private_incentive_restrictions ? 0 : negative_infinity())>[use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_discrete_dist : 1] beta_calendar_effect;
  vector<lower = (use_private_incentive_restrictions ? 0 : negative_infinity())>[use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_discrete_dist : 1] beta_bracelet_effect;
  
  matrix[use_cluster_effects ? num_clusters : 0, num_actual_treatments] structural_beta_cluster_raw;
  row_vector<lower = 0>[use_cluster_effects ? num_actual_treatments : 0] structural_beta_cluster_sd;
  
  matrix[use_county_effects ? num_counties : 0, num_actual_treatments] structural_beta_county_raw;
  row_vector<lower = 0>[use_county_effects ? num_actual_treatments : 0] structural_beta_county_sd;
  
  // Salience
  
  real<lower = 0> beta_salience;
  real<lower = 0> dist_beta_salience;
  real<lower = 0> dist_quadratic_beta_salience;
  
  // Name Matched
  
  vector[num_treatments_name_matched] beta_nm_effect;
  
  matrix[use_cluster_effects && use_name_matched_obs ? num_clusters : 0, num_actual_treatments] beta_nm_effect_cluster_raw;
  row_vector<lower = 0>[use_cluster_effects && use_name_matched_obs ? num_actual_treatments : 0] beta_nm_effect_cluster_sd;
  
  matrix[use_county_effects && use_name_matched_obs ? num_clusters : 0, num_actual_treatments] beta_nm_effect_county_raw;
  row_vector<lower = 0>[use_county_effects && use_name_matched_obs ? num_actual_treatments : 0] beta_nm_effect_county_sd;
  
  // V Mixture
  
  vector<lower = 0.5, upper = 1>[num_v_mix - 1] recursive_lambda_v_mix;
  
  vector<lower = 0>[num_v_mix - 1] v_mix_mean_diff;
  vector<lower = 0>[num_v_mix - 1] v_mix_sd;
  
  // Reputational Returns
  
  real v_mu;
  
  row_vector<lower = 0>[suppress_reputation ? 0 : num_actual_treatments] mu_rep_raw;
  // real<lower = 0, upper = 1> v_sd;
  vector<lower = 0>[suppress_shocks ? 0 : num_shocks - 1] ub_ur_sd;
  cholesky_factor_corr[suppress_shocks ? 0 : num_shocks] L_all_u_corr;
  
  matrix[use_mu_cluster_effects && !suppress_reputation ? num_clusters : 0, num_actual_treatments] mu_cluster_effects_raw;
  row_vector<lower = 0>[use_mu_cluster_effects && !suppress_reputation ? num_actual_treatments : 0] mu_cluster_effects_sd;
  
  matrix[use_mu_county_effects && !suppress_reputation ? num_counties : 0, num_actual_treatments] mu_county_effects_raw;
  row_vector<lower = 0>[use_mu_county_effects && !suppress_reputation ? num_actual_treatments : 0] mu_county_effects_sd;
  
  // Parameteric Cost Model
  
  vector<lower = 0>[num_treatments_param_kappa] dist_cost_k_raw;
  
  // Linear Parametric Cost
  
  vector<lower = (use_dist_salience ? 0 : negative_infinity())>[num_treatments_param] dist_beta_v; // Linear distance*treatment effects
  
  matrix[use_param_dist_cluster_effects ? num_clusters : 0, num_treatments_param] dist_beta_cluster_raw;
  row_vector<lower = 0>[use_param_dist_cluster_effects ? num_treatments_param : 0] dist_beta_cluster_sd;
  
  matrix[use_param_dist_county_effects ? num_counties : 0, num_treatments_param] dist_beta_county_raw;
  row_vector<lower = 0>[use_param_dist_county_effects ? num_treatments_param : 0] dist_beta_county_sd;
  
  // Quadratic Cost Model
  
  vector<lower = 0>[num_treatments_param_quadratic] dist_quadratic_beta_v; // Quadratic distance*treatment effects
  
  // Semiparameteric Cost Model
  
  matrix[num_treatments_semiparam, num_knots_v] u_splines_v_raw;
  real<lower = 0> u_splines_v_sigma;
  
  // Assigned Distance Model
  
  vector<lower = 0>[num_discrete_dist] group_dist_mean;
  vector<lower = 0>[num_discrete_dist] group_dist_sd;
  
  matrix<lower = 0>[num_clusters, num_discrete_dist - 1] missing_cluster_standard_dist; 
}

transformed parameters {
  vector[num_actual_treatments] beta; 
  vector[num_actual_treatments] structural_treatment_effect; 
  vector[num_treatments_name_matched] structural_treatment_effect_nm;
  vector[num_clusters] structural_cluster_benefit_cost;
  matrix[num_clusters, num_actual_treatments] structural_beta_cluster = rep_matrix(0, num_clusters, num_actual_treatments);
  matrix[num_clusters, num_actual_treatments] structural_beta_cluster_nm = rep_matrix(0, num_clusters, num_actual_treatments);
  matrix[num_counties, num_actual_treatments] structural_beta_county = rep_matrix(0, num_counties, num_actual_treatments);
  matrix[num_counties, num_actual_treatments] structural_beta_county_nm = rep_matrix(0, num_counties, num_actual_treatments);
  
  row_vector<lower = 0>[suppress_reputation ? 0 : num_treatments] mu_rep = mu_rep_raw;
  matrix<lower = 0>[!suppress_reputation ? num_clusters : 0, num_treatments] cluster_mu_rep;
  
  vector[num_clusters] structural_cluster_obs_v = rep_vector(0, num_clusters);
  vector[use_name_matched_obs ? num_clusters : 0] structural_cluster_obs_v_nm;
  matrix<lower = 0, upper = 1>[num_v_mix, num_clusters] structural_cluster_takeup_prob;
  matrix<lower = 0, upper = 1>[num_v_mix, use_name_matched_obs ? num_clusters : 0] structural_cluster_takeup_prob_nm;
  
  simplex[num_v_mix] lambda_v_mix;
  ordered[num_v_mix - 1] v_mix_mean = cumulative_sum(v_mix_mean_diff);
  
  vector<lower = 0>[num_treatments_param_kappa] dist_cost_k;
  vector[num_treatments] linear_dist_cost = rep_vector(0, num_treatments);
  vector[num_treatments] quadratic_dist_cost = rep_vector(0, num_treatments);
  matrix[num_clusters, num_treatments] cluster_linear_dist_cost = rep_matrix(0, num_clusters, num_treatments);
  matrix[num_treatments, num_knots_v] u_splines_v = rep_matrix(0, num_treatments, num_knots_v);
  vector[num_clusters] cluster_dist_cost = rep_vector(0, num_clusters);
  
  matrix<lower = 0>[num_clusters, num_discrete_dist] cf_cluster_standard_dist; 
  
  // cholesky_factor_cov[num_shocks] L_all_u_vcov;
  matrix[num_shocks, num_shocks] L_all_u_vcov;
  real total_error_sd;
  
  for (cluster_index in 1:num_clusters) {
    int dist_group_pos = 1;
    
    for (dist_group_index in 1:num_discrete_dist) {
      if (cluster_assigned_dist_group[cluster_index] == dist_group_index) {
        cf_cluster_standard_dist[cluster_index, dist_group_index] = cluster_standard_dist[cluster_index];
      } else {
        cf_cluster_standard_dist[cluster_index, dist_group_index] = missing_cluster_standard_dist[cluster_index, dist_group_pos];
        
        dist_group_pos += 1;
      }
    }
  }
 
  if (use_cost_model == COST_MODEL_TYPE_DISCRETE) {
    for (dist_index in 1:num_discrete_dist) {
      beta[(num_treatments * (dist_index - 1) + 1):(num_treatments * dist_index)] = 
        [ beta_control[dist_index], beta_ink_effect[dist_index], beta_calendar_effect[dist_index], beta_bracelet_effect[dist_index] ]';
    }
  } else {
    beta = [ beta_control[1], beta_ink_effect[1], beta_calendar_effect[1], beta_bracelet_effect[1] ]';
  }
  
  structural_treatment_effect = restricted_treatment_map_design_matrix * beta;
  
  if (num_treatments_name_matched > 0) {
    structural_treatment_effect_nm = treatment_map_design_matrix * beta_nm_effect;
    structural_cluster_obs_v_nm = rep_vector(0, num_clusters);
  }
  
  if (num_v_mix > 1) {
    lambda_v_mix[1] = recursive_lambda_v_mix[1];
    
    for (mix_index in 2:(num_v_mix - 1)) {
      lambda_v_mix[mix_index] = recursive_lambda_v_mix[mix_index] * (1 - sum(lambda_v_mix[1:(mix_index - 1)]));
    }
    
    lambda_v_mix[num_v_mix] = 1 - sum(lambda_v_mix[1:(num_v_mix - 1)]);
  } else {
    lambda_v_mix[1] = 1;
  }

  // Levels: control ink calendar bracelet
 
  if (!suppress_reputation) { 
    mu_rep[2] += mu_rep[1];
    mu_rep[3] += mu_rep[1];
    mu_rep[4] += mu_rep[1];
    
    cluster_mu_rep = rep_matrix(mu_rep, num_clusters);
    
    if (use_mu_cluster_effects) {
      matrix[num_clusters, num_treatments] mu_cluster_effects =  mu_cluster_effects_raw .* rep_matrix(mu_cluster_effects_sd, num_clusters);
      
      cluster_mu_rep = cluster_mu_rep .* exp(mu_cluster_effects);
    }
    
    if (use_mu_county_effects) {
      matrix[num_counties, num_treatments] mu_county_effects =  mu_county_effects_raw .* rep_matrix(mu_county_effects_sd, num_counties);
      
      cluster_mu_rep = cluster_mu_rep .* exp(mu_county_effects[cluster_county_id]);
    } 
  }
  
  if (use_salience_effect) {
    structural_treatment_effect += beta_salience * mu_rep';
  }
  
  if (use_cost_model == COST_MODEL_TYPE_PARAM_KAPPA) { 
    dist_cost_k = dist_cost_k_raw;
   
    if (use_cost_k_restrictions) { 
      dist_cost_k[2] += dist_cost_k[1];
      dist_cost_k[3] += dist_cost_k[1];
      dist_cost_k[4] += dist_cost_k[1];
    }
    
    cluster_dist_cost = param_kappa_dist_cost(cluster_standard_dist, dist_cost_k[cluster_assigned_treatment]);
  } else if (use_cost_model != COST_MODEL_TYPE_DISCRETE) {
    if (in_array(use_cost_model, { COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE, COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE, COST_MODEL_TYPE_SEMIPARAM_SALIENCE })) {
      linear_dist_cost = rep_vector(dist_beta_v[1], num_treatments) + dist_beta_salience * mu_rep';
      
      if (use_param_dist_cluster_effects) {
        cluster_linear_dist_cost += rep_matrix(dist_beta_cluster_raw[, 1] * dist_beta_cluster_sd[1], num_treatments);
      }
      
      if (use_param_dist_county_effects) {
        cluster_linear_dist_cost += rep_matrix((dist_beta_county_raw[, 1] * dist_beta_county_sd[1])[cluster_county_id], num_treatments);
      }
      
      if (use_cost_model == COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE) {
        quadratic_dist_cost = rep_vector(dist_quadratic_beta_v[1], num_treatments) + dist_quadratic_beta_salience * mu_rep';
      }
    } else {
      if (use_single_cost_model) {
        linear_dist_cost = rep_vector(dist_beta_v[1], num_treatments);
        
        if (use_param_dist_cluster_effects) {
          cluster_linear_dist_cost += rep_matrix(dist_beta_cluster_raw[, 1] * dist_beta_cluster_sd[1], num_treatments);
        }
        
        if (use_param_dist_county_effects) {
          cluster_linear_dist_cost += rep_matrix((dist_beta_county_raw[, 1] * dist_beta_county_sd[1])[cluster_county_id], num_treatments);
        }
        
        if (num_treatments_param_quadratic > 0) {
          quadratic_dist_cost = rep_vector(dist_quadratic_beta_v[1], num_treatments);
        }
      } else {
        linear_dist_cost = treatment_map_design_matrix * dist_beta_v;
        
        if (use_param_dist_cluster_effects) {
          cluster_linear_dist_cost += dist_beta_cluster_raw .* rep_matrix(dist_beta_cluster_sd, num_clusters);
        }
        
        if (use_param_dist_county_effects) {
          cluster_linear_dist_cost += (dist_beta_county_raw .* rep_matrix(dist_beta_county_sd, num_counties))[cluster_county_id];
        }
        
        if (num_treatments_param_quadratic > 0) {
          quadratic_dist_cost = treatment_map_design_matrix * dist_quadratic_beta_v;
        }
      }
      
      cluster_linear_dist_cost += rep_matrix(linear_dist_cost', num_clusters);
    } 
    
    if (use_semiparam) {
      for (treatment_index in 1:num_treatments_semiparam) {
        if (use_dist_salience) {
          u_splines_v[treatment_index] = u_splines_v_raw[1] * u_splines_v_sigma * mu_rep[treatment_index];
        } else {
          u_splines_v[treatment_index] = u_splines_v_raw[treatment_index] * u_splines_v_sigma;
        }
      }        
    } 
      
    cluster_dist_cost = param_dist_cost(cluster_standard_dist, 
                                        to_vector(cluster_linear_dist_cost)[long_cluster_by_treatment_index],
                                        quadratic_dist_cost[cluster_assigned_treatment],
                                        u_splines_v[cluster_assigned_treatment],
                                        Z_splines_v[cluster_assigned_treatment]);
  }
  
  if (use_cost_model != COST_MODEL_TYPE_DISCRETE) {
    structural_cluster_benefit_cost = structural_treatment_effect[cluster_assigned_treatment] - cluster_dist_cost;
  } else {
    structural_cluster_benefit_cost = structural_treatment_effect[cluster_assigned_dist_group_treatment];
  }
  
  if (use_cluster_effects) {
    vector[num_clusters] cluster_effects;
    structural_beta_cluster = structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_clusters);
    
    cluster_effects = rows_dot_product(cluster_treatment_design_matrix, structural_beta_cluster); 
    structural_cluster_benefit_cost += cluster_effects;
  }
  
  if (use_county_effects) {
    vector[num_clusters] county_effects;
    structural_beta_county = structural_beta_county_raw .* rep_matrix(structural_beta_county_sd, num_counties);
    
    county_effects = rows_dot_product(cluster_treatment_design_matrix, structural_beta_cluster[cluster_county_id]); 
    structural_cluster_benefit_cost += county_effects;
  }

  if (suppress_reputation) {
    structural_cluster_obs_v = - structural_cluster_benefit_cost;
  } else {
    for (cluster_index in 1:num_clusters) {
      structural_cluster_obs_v[cluster_index] = algebra_solver(v_fixedpoint_solution_normal,
                                                               [ - structural_cluster_benefit_cost[cluster_index] ]',
                                                               prepare_solver_theta(structural_cluster_benefit_cost[cluster_index], 
                                                                                    cluster_mu_rep[cluster_index, cluster_assigned_treatment[cluster_index]],
                                                                                    use_shifting_v_dist ? v_mu : 0,
                                                                                    lambda_v_mix,
                                                                                    v_mix_mean,
                                                                                    1,
                                                                                    v_mix_sd),
                                                               { 0.0 },
                                                               { num_v_mix }, 
                                                               1e-10,
                                                               1e-5,
                                                               1e6)[1];
                                                               
    }
  }
  
  if (use_name_matched_obs) {
    vector[num_clusters] cluster_effects_nm = rep_vector(0, num_clusters);
    vector[num_clusters] county_effects_nm = rep_vector(0, num_clusters);
    
    if (use_cluster_effects) {
      structural_beta_cluster_nm = beta_nm_effect_cluster_raw .* rep_matrix(beta_nm_effect_cluster_sd, num_clusters);
      cluster_effects_nm = rows_dot_product(cluster_treatment_design_matrix, structural_beta_cluster_nm); 
    }
    
    if (use_county_effects) {
      structural_beta_county_nm = beta_nm_effect_county_raw .* rep_matrix(beta_nm_effect_county_sd, num_counties);
      county_effects_nm = rows_dot_product(cluster_treatment_design_matrix, structural_beta_county_nm[cluster_county_id]); 
    }
    
    structural_cluster_obs_v_nm = structural_cluster_obs_v 
                                - (structural_treatment_effect_nm[cluster_assigned_treatment] + cluster_effects_nm + county_effects_nm); 
  }
  
  if (!suppress_shocks) {
    vector[num_shocks] all_u_sd = append_row(1, ub_ur_sd);
    
    L_all_u_vcov = diag_pre_multiply(all_u_sd, L_all_u_corr);
    total_error_sd = sqrt(sum(L_all_u_vcov * L_all_u_vcov'));
  } else {
    total_error_sd = 1;
  }
  
  for (mix_index in 1:num_v_mix) {
    if (mix_index == 1) {
      structural_cluster_takeup_prob[1] = Phi(- structural_cluster_obs_v / total_error_sd)';
      structural_cluster_takeup_prob_nm[1] = Phi(- structural_cluster_obs_v_nm / total_error_sd)';
      
      // structural_obs_takeup_prob[1] = structural_cluster_takeup_prob[1, obs_cluster_id]; // BUGBUG no name matching yet 
    } else {
      structural_cluster_takeup_prob[mix_index] = Phi(- (structural_cluster_obs_v - v_mix_mean[mix_index - 1]) / v_mix_sd[mix_index - 1])';
      structural_cluster_takeup_prob_nm[mix_index] = Phi(- (structural_cluster_obs_v_nm - v_mix_mean[mix_index - 1]) / v_mix_sd[mix_index - 1])';
      
      // structural_obs_takeup_prob[mix_index] = structural_cluster_takeup_prob[mix_index, obs_cluster_id]; // BUGBUG no name matching yet 
    }
  }
}

model {
  beta_control[1] ~ normal(0, 5);
  
  if (use_cost_model == COST_MODEL_TYPE_DISCRETE) {
    beta_control[2:num_discrete_dist] ~ normal(0, 1);
  }
  
  beta_ink_effect ~ normal(0, 1);
  beta_calendar_effect ~ normal(0, 1);
  beta_bracelet_effect ~ normal(0, 1);
  
  beta_salience ~ normal(0, 1);
  dist_beta_salience ~ normal(0, 1);
  dist_quadratic_beta_salience ~ normal(0, 1);
  
  v_mu ~ normal(0, 1);
  
  // beta_salience_nm_effect ~ normal(0, 1);
  
  if (num_treatments_name_matched > 0) {
    beta_nm_effect ~ normal(0, 1);
  }

  if (!suppress_reputation) { 
    mu_rep_raw ~ normal(0, mu_rep_sd);
    
    if (use_mu_cluster_effects) {
      to_vector(mu_cluster_effects_raw) ~ normal(0, 1);
      mu_cluster_effects_sd ~ normal(0, 0.05);
    }
    
    if (use_mu_county_effects) {
      to_vector(mu_county_effects_raw) ~ normal(0, 1);
      mu_county_effects_sd ~ normal(0, 0.1);
    }
  }

  if (num_treatments_param_kappa > 0) {
    dist_cost_k_raw ~ normal(0, 1);
  } 
  
  u_splines_v_sigma ~ normal(0, u_splines_v_sigma_sd);
  
  if (num_treatments_param > 0) {
    dist_beta_v ~ normal(0, 1);
    
    if (use_param_dist_cluster_effects) {
      to_vector(dist_beta_cluster_raw) ~ normal(0, 1);
      dist_beta_cluster_sd ~ normal(0, 0.25);
    }
    
    if (use_param_dist_county_effects) {
      to_vector(dist_beta_county_raw) ~ normal(0, 1);
      dist_beta_county_sd ~ normal(0, 0.25);
    }
  }
  
  if (num_treatments_semiparam > 0) {
    to_vector(u_splines_v_raw) ~ normal(0, 1);
  }
  
  if (num_treatments_param_quadratic > 0) {
    dist_quadratic_beta_v ~ normal(0, 1);
  }
  
  if (use_cluster_effects) {
    to_vector(structural_beta_cluster_raw) ~ normal(0, 1);
    structural_beta_cluster_sd ~ normal(0, 0.25);
    
    if (use_name_matched_obs) {
      to_vector(beta_nm_effect_cluster_raw) ~ normal(0, 1);
      beta_nm_effect_cluster_sd ~ normal(0, 0.25);
    }
  }
  
  if (use_county_effects) {
    to_vector(structural_beta_county_raw) ~ normal(0, 1);
    structural_beta_county_sd ~ normal(0, 0.25);
    
    if (use_name_matched_obs) {
      to_vector(beta_nm_effect_county_raw) ~ normal(0, 1);
      beta_nm_effect_county_sd ~ normal(0, 0.25);
    }
  }
  
  if (!suppress_shocks) {
    ub_ur_sd ~ normal(0, 1);
    L_all_u_corr ~ lkj_corr_cholesky(2);
  }
  
  if (num_v_mix > 1) {
    v_mix_mean_diff ~ normal(0, 1);
    v_mix_sd ~ normal(0, 1);
  }
  
  if (num_v_mix == 1) {
    if (use_binomial) {
      if (use_name_matched_obs) {
        reject("Cannot use binomial model if using name-matched observations.");
      }
      
      cluster_takeup_count[included_clusters] ~ binomial(cluster_size[included_clusters], structural_cluster_takeup_prob[1, included_clusters]);
    } else {
      takeup[included_monitored_obs] ~ bernoulli(structural_cluster_takeup_prob[1, obs_cluster_id[included_monitored_obs]]);
      
      if (use_name_matched_obs) {
        takeup[included_name_matched_obs] ~ bernoulli(structural_cluster_takeup_prob_nm[1, obs_cluster_id[included_name_matched_obs]]);
      }
    }
  } else {
    if (use_binomial) {
      if (use_name_matched_obs) {
        reject("Cannot use binomial model if using name-matched observations.");
      }
      
      cluster_takeup_count[included_clusters] ~ mixed_binomial(lambda_v_mix, cluster_size[included_clusters], structural_cluster_takeup_prob[, included_clusters]);
    } else {
      takeup[included_monitored_obs] ~ mixed_bernoulli(lambda_v_mix, structural_cluster_takeup_prob[, obs_cluster_id[included_monitored_obs]]);
      
      if (use_name_matched_obs) {
        takeup[included_name_matched_obs] ~ mixed_bernoulli(lambda_v_mix, structural_cluster_takeup_prob_nm[, obs_cluster_id[included_name_matched_obs]]);
      }
    }
  }
  
  group_dist_mean ~ normal(0, 1);
  group_dist_sd ~ normal(0, 1);
  
  for (cluster_index in 1:num_clusters) {
    int dist_group_pos = 1;
    
    cluster_standard_dist[cluster_index] ~ normal(group_dist_mean[cluster_assigned_dist_group[cluster_index]], group_dist_sd[cluster_assigned_dist_group[cluster_index]]) T[0,];
    
    for (dist_group_index in 1:num_discrete_dist) {
      if (dist_group_index != cluster_assigned_dist_group[cluster_index]) {
        missing_cluster_standard_dist[cluster_index, dist_group_pos] ~ normal(group_dist_mean[dist_group_index], group_dist_sd[dist_group_index]) T[0,];
        dist_group_pos += 1;
      }
    }
  }
}

generated quantities { 
  matrix[num_shocks, num_shocks] all_u_corr;
 
  vector[num_clusters] cluster_cf_benefit_cost[num_actual_treatments]; 
  
  vector[generate_rep ? num_clusters : 0] cluster_rep_benefit_cost; 
  
  // Cross Validation
  vector[use_binomial || cluster_log_lik ? num_included_clusters : num_included_obs] log_lik = 
    rep_vector(negative_infinity(), use_binomial || cluster_log_lik ? num_included_clusters : num_included_obs); 
  vector[use_binomial || cluster_log_lik ? num_excluded_clusters : num_excluded_obs] log_lik_heldout = 
    rep_vector(negative_infinity(), use_binomial || cluster_log_lik ? num_excluded_clusters : num_excluded_obs); 
    
  if (num_shocks > 0) {
    all_u_corr = L_all_u_corr * L_all_u_corr';
  }
  
  {
    int treatment_cluster_pos = 1;
    int cluster_treatment_cf_pos = 1;
    
    for (treatment_index in 1:num_actual_treatments) {
      vector[num_clusters] treatment_dist_cost;
     
      if (use_cost_model == COST_MODEL_TYPE_DISCRETE) {
        treatment_dist_cost = rep_vector(0, num_clusters);
      } else {                                                            
        treatment_dist_cost = param_dist_cost(cluster_standard_dist, 
                                              cluster_linear_dist_cost[, treatment_index],
                                              rep_vector(quadratic_dist_cost[treatment_index], num_clusters),
                                              u_splines_v[rep_array(treatment_index, num_clusters)],
                                              Z_splines_v[rep_array(treatment_index, num_clusters)]);
      }
                                                                   
      cluster_cf_benefit_cost[treatment_index] =
        // Not using "restricted" design matrix because restrictions are only on the top-level parameters not village and county level effects
        (structural_beta_cluster + structural_beta_cluster_nm + (structural_beta_county + structural_beta_county_nm)[cluster_county_id]) * treatment_map_design_matrix[treatment_index]'
        + rep_vector(structural_treatment_effect[treatment_index], num_clusters) 
        - treatment_dist_cost;
    }
  }
  
  // Cross Validation 
  
  if (use_binomial) {
    for (cluster_index_index in 1:num_included_clusters) {
      int cluster_index = included_clusters[cluster_index_index];
      
      if (num_v_mix == 1) {
        log_lik[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], structural_cluster_takeup_prob[1, cluster_index]);
      } else {
        log_lik[cluster_index_index] = mixed_binomial_lpmf({ cluster_takeup_count[cluster_index] } | lambda_v_mix, { cluster_size[cluster_index] }, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
    }
    
    for (cluster_index_index in 1:num_excluded_clusters) {
      int cluster_index = excluded_clusters[cluster_index_index];
      
      if (num_v_mix == 1) {
        log_lik_heldout[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], structural_cluster_takeup_prob[1, cluster_index]);
      } else {
        log_lik_heldout[cluster_index_index] = mixed_binomial_lpmf({ cluster_takeup_count[cluster_index] } | lambda_v_mix, { cluster_size[cluster_index] }, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
    }
  } else {
    vector[num_clusters] temp_log_lik = rep_vector(0, num_clusters);
    
    for (obs_index_index in 1:num_included_monitored_obs) {
      int obs_index = included_monitored_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      real curr_log_lik;
      
      if (num_v_mix == 1) {
        curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[1, cluster_index:cluster_index]);
      } else {
        curr_log_lik = mixed_bernoulli_lpmf({ takeup[obs_index] } | lambda_v_mix, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
      
      if (cluster_log_lik) {
        temp_log_lik[cluster_index] += curr_log_lik;
      } else {
        log_lik[obs_index_index] = curr_log_lik;
      }
    }
    
    for (obs_index_index in 1:num_included_name_matched_obs) {
      int obs_index = included_name_matched_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      real curr_log_lik;
      
      if (num_v_mix == 1) {
        curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob_nm[1, cluster_index:cluster_index]);
      } else {
        curr_log_lik = mixed_bernoulli_lpmf({ takeup[obs_index] } | lambda_v_mix, structural_cluster_takeup_prob_nm[, cluster_index:cluster_index]);
      }
      
      if (cluster_log_lik) {
        temp_log_lik[cluster_index] += curr_log_lik;
      } else {
        log_lik[obs_index_index] = curr_log_lik;
      }
    }
    
    for (obs_index_index in 1:num_excluded_monitored_obs) {
      int obs_index = excluded_monitored_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      real curr_log_lik;
      
      if (num_v_mix == 1) {
        curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[1, cluster_index:cluster_index]);
      } else {
        curr_log_lik = mixed_bernoulli_lpmf({ takeup[obs_index] } | lambda_v_mix, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
      
      if (cluster_log_lik) {
        temp_log_lik[cluster_index] += curr_log_lik;
      } else {
        log_lik_heldout[obs_index_index] = curr_log_lik;
      }
    }
    
    for (obs_index_index in 1:num_excluded_name_matched_obs) {
      int obs_index = excluded_name_matched_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      real curr_log_lik;
      
      if (num_v_mix == 1) {
        curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob_nm[1, cluster_index:cluster_index]);
      } else {
        curr_log_lik = mixed_bernoulli_lpmf({ takeup[obs_index] } | lambda_v_mix, structural_cluster_takeup_prob_nm[, cluster_index:cluster_index]);
      }
      
      if (cluster_log_lik) {
        temp_log_lik[cluster_index] += curr_log_lik;
      } else {
        log_lik_heldout[obs_index_index] = curr_log_lik;
      }
    }
    
    if (cluster_log_lik) {
      log_lik = temp_log_lik[included_clusters];
      log_lik_heldout = temp_log_lik[excluded_clusters];
    }
  }
 
  // Simulation 
  
  if (generate_rep) {
    for (cluster_index in 1:num_clusters) {
      vector[num_actual_treatments] rep_beta_cluster = 
        use_cluster_effects ? to_vector(normal_rng(rep_array(0, num_actual_treatments), structural_beta_cluster_sd')) : rep_vector(0, num_actual_treatments);
        
      vector[num_actual_treatments] rep_beta_county = 
        use_county_effects ? to_vector(normal_rng(rep_array(0, num_actual_treatments), structural_beta_county_sd')) : rep_vector(0, num_actual_treatments);
      
      if (use_cost_model != COST_MODEL_TYPE_DISCRETE) {
        real rep_dist_cost;
        
        if (use_cost_model == COST_MODEL_TYPE_PARAM_KAPPA) { 
          rep_dist_cost = param_kappa_dist_cost([ cluster_standard_dist[cluster_index] ]', [ dist_cost_k[cluster_assigned_treatment[cluster_index]] ]')[1];
        } else {
         real rep_linear_dist_cost = linear_dist_cost[cluster_assigned_treatment[cluster_index]];
          
          if (in_array(use_cost_model, { COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE, COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE, COST_MODEL_TYPE_SEMIPARAM_SALIENCE })) {
            if (use_param_dist_cluster_effects) {
              rep_linear_dist_cost += normal_rng(0, dist_beta_cluster_sd[1]); 
            }
            
            if (use_param_dist_county_effects) {
              rep_linear_dist_cost += normal_rng(0, dist_beta_county_sd[1]); 
            }
          } else {
            if (use_single_cost_model) {
              if (use_param_dist_cluster_effects) {
                rep_linear_dist_cost += normal_rng(0, dist_beta_cluster_sd[1]); 
              }
              
              if (use_param_dist_county_effects) {
                rep_linear_dist_cost += normal_rng(0, dist_beta_county_sd[1]); 
              }
            } else {
              vector[num_treatments] rep_dist_beta = rep_vector(0, num_treatments);
              
              if (use_param_dist_cluster_effects) {
                rep_dist_beta += to_vector(normal_rng(rep_array(0, num_actual_treatments), dist_beta_cluster_sd'));
              }
              
              if (use_param_dist_county_effects) {
                rep_dist_beta += to_vector(normal_rng(rep_array(0, num_actual_treatments), dist_beta_county_sd'));
              }
              
              rep_linear_dist_cost += treatment_map_design_matrix[cluster_assigned_treatment[cluster_index]] * rep_dist_beta;
            }
          } 
            
          rep_dist_cost = param_dist_cost([ cluster_standard_dist[cluster_index] ]', 
                                          [ rep_linear_dist_cost ]',
                                          [ quadratic_dist_cost[cluster_assigned_treatment[cluster_index]] ]',
                                          to_matrix(u_splines_v[cluster_assigned_treatment[cluster_index]]),
                                          to_matrix(Z_splines_v[cluster_assigned_treatment[cluster_index]]))[1];
        }
        
        cluster_rep_benefit_cost[cluster_index] = structural_treatment_effect[cluster_assigned_treatment[cluster_index]] 
          + treatment_map_design_matrix[cluster_assigned_treatment[cluster_index]] * (rep_beta_cluster + rep_beta_county) - rep_dist_cost;
      } else {
        cluster_rep_benefit_cost[cluster_index] = structural_treatment_effect[cluster_assigned_dist_group_treatment[cluster_index]]
          + treatment_map_design_matrix[cluster_assigned_dist_group_treatment[cluster_index]] * (rep_beta_cluster + rep_beta_county);
      }
    }
  }
}
