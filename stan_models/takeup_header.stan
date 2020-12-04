functions {
#include /../multilvlr/util.stan
  
  real reputational_returns_normal(real v, vector lambda, vector mix_mean, vector mix_sd) {
    int num_mix = num_elements(lambda);
    real rep = 0; 
    
    for (mix_index in 1:num_mix) {
      real mix_Phi_v = Phi_approx((v - mix_mean[mix_index]) / mix_sd[mix_index]);
      
      rep += lambda[mix_index] * exp(normal_lpdf(v | mix_mean[mix_index], mix_sd[mix_index])) / (mix_Phi_v * (1 - mix_Phi_v));
    }
    
    return rep;
  }
  
  
  vector rep_normal_1std(vector v) {
    int num_v = num_elements(v);
    vector[num_v] result;
   
    for (v_index in 1:num_v) {
      real phi_v = exp(std_normal_lpdf(v[v_index])); 
      real Phi_v = Phi_approx(v[v_index]); 
      real Phi_1mPhi_v = Phi_v * (1- Phi_v);
      
      result[v_index] = (phi_v / Phi_1mPhi_v^2) * (-(v[v_index] * Phi_1mPhi_v) + (phi_v * (2 * Phi_v - 1)));
    } 
    
    return result;
  }

  vector social_multiplier(vector delta_1st, real mu) {
    return - 1.0 ./ (1 + mu * delta_1st);
  }

  vector expect_y_partial_bbar(vector v, vector sm) {
    int num_v = num_elements(v);
    vector[num_v] result = sm;
   
    for (v_index in 1:num_v) {
      real phi_v = exp(std_normal_lpdf(v[v_index])); 
      
      result[v_index] *= - phi_v; 
    } 
    
    return result;
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
  
  real normal_density(real x,          // Function argument
                    real xc,         // Complement of function argument
                                     //  on the domain (defined later)
                    real[] theta,    // parameters
                    real[] x_r,      // data (real)
                    int[] x_i) {     // data (integer)
    real mu = theta[1];
    real sigma = theta[2];
  
    return 1 / (sqrt(2 * pi()) * sigma) * exp(-0.5 * ((x - mu) / sigma)^2);
  } 
  
  real expected_delta(real u, real xc, real[] theta, data real[] x_r, data int[] x_i) {
    real cutoff = theta[1] - u;
    real u_sd = theta[2];

    // return reputational_returns_normal(v - u, [ 1 ]', [ 0 ]', [ 1 ]') * exp(normal_lpdf(u | 0, u_sd));
    return exp(normal_lpdf(cutoff | 0, 1) - normal_lcdf(cutoff| 0, 1) - normal_lccdf(cutoff | 0, 1) + normal_lpdf(u | 0, u_sd));
  }
  
  vector v_fixedpoint_solution_normal(vector model_param, vector theta, data real[] x_r, data int[] x_i) {
    real v_cutoff = model_param[1];
    
    real benefit_cost = theta[1];
    real mu = theta[2];
    
    int num_v_mix = x_i[1];
    int use_u_in_delta = x_i[2];
    
    vector[num_v_mix] lambda = theta[3:(3 + num_v_mix - 1)];
    vector[num_v_mix] mix_mean = theta[(3 + num_v_mix):(3 + 2 * num_v_mix - 1)];
    vector[num_v_mix] mix_sd = theta[(3 + 2 * num_v_mix):(3 + 3 * num_v_mix - 1)];
    real u_sd = theta[3 + 3 * num_v_mix]; 
    
    if (use_u_in_delta && u_sd > 0) {
      real delta = integrate_1d(expected_delta, negative_infinity(), positive_infinity(), { v_cutoff, u_sd }, x_r, x_i, 0.01);

      return [ v_cutoff + benefit_cost + mu * delta ]';
    } else {
      return [ v_cutoff + benefit_cost + mu * reputational_returns_normal(v_cutoff, lambda, mix_mean, mix_sd) ]';
    }
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
  
  vector prepare_solver_theta(real benefit_cost, real mu_rep, real v_mu, vector lambda, vector mix_mean, real v_sd, vector mix_sd, real u_sd) { //, vector u_shocks) {
    int num_v_mix = num_elements(lambda);
    vector[2 + 3 * num_v_mix + 1] solver_theta;
    
    solver_theta[1:2] = [ benefit_cost, mu_rep ]';
    solver_theta[3:(3 + num_v_mix - 1)] = lambda;
    solver_theta[(3 + num_v_mix):(3 + num_v_mix + num_v_mix - 1)] = append_row(v_mu, mix_mean);
    solver_theta[(3 + num_v_mix + num_v_mix):(3 + 3 * num_v_mix - 1)] = append_row(v_sd, mix_sd);
    solver_theta[3 + 3 * num_v_mix] = u_sd; 
 
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
  int<lower = 0, upper = 1> use_u_in_delta;
  int<lower = 0, upper = 1> generate_rep;
  int<lower = 0, upper = 1> generate_sim;
  int<lower = 0, upper = 1> predict_prior;

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