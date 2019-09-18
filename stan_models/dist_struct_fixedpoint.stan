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
    vector[num_v_mix] mix_mean = append_row(0, theta[(3 + num_v_mix):(3 + num_v_mix + num_v_mix - 2)]);
    vector[num_v_mix] mix_sd = theta[(3 + num_v_mix + num_v_mix - 1):(3 + 3 * num_v_mix - 2)];
    
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
  
  vector prepare_solver_theta(real benefit_cost, real mu_rep, vector lambda, vector mix_mean, real v_sd, vector mix_sd) {
    int num_v_mix = num_elements(lambda);
    vector[2 + 2 * num_v_mix + (num_v_mix - 1)] solver_theta;
    
    solver_theta[1:2] = [ benefit_cost, mu_rep ]';
    solver_theta[3:(3 + num_v_mix - 1)] = lambda;
    solver_theta[(3 + num_v_mix):(3 + num_v_mix + num_v_mix - 2)] = mix_mean;
    solver_theta[(3 + num_v_mix + num_v_mix - 1):(3 + 3 * num_v_mix - 2)] = append_row(v_sd, mix_sd);
 
    return solver_theta; 
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
  
  int<lower = 1> num_v_mix;
  
  int<lower = 0, upper = 1> use_binomial;
  int<lower = 0, upper = 1> use_cluster_effects;
  int<lower = 0, upper = 1> use_mu_cluster_effects;
  int<lower = 0, upper = 1> use_cost_k_restrictions;
  int<lower = 0, upper = 1> use_private_incentive_restrictions;
  int<lower = 0, upper = 1> use_salience_effect;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> use_cost_model;
  int<lower = 0, upper = 1> suppress_reputation;
  int<lower = 0, upper = 1> suppress_shocks; 
  int<lower = 0, upper = 1> simulate_new_data; 

  int<lower = 1> num_obs; // Actual observations
  int<lower = 1> num_clusters;
  int<lower = 1> num_treatments;
  
  int<lower = 1, upper = num_clusters> obs_cluster_id[num_obs];
 
  int<lower = 0, upper = 1> takeup[num_obs]; // Observed outcome variable
  vector[num_clusters] cluster_standard_dist; // Standardized distance to treatment
  
  int<lower = 1, upper = num_treatments> cluster_assigned_treatment[num_clusters]; // Actual assigned treatments 
  
  // Reputation
  
  real<lower = 0> mu_rep_sd;
  
  // Semiparametric Cost Model (Splines)
  
  int<lower = 3> num_knots_v;
  matrix[num_clusters, num_knots_v] Z_splines_v; 
  
  real<lower = 0> u_splines_v_sigma_sd;
  
  // K-fold CV 
  
  int<lower = 0> num_excluded_clusters;
  int<lower = 1, upper = num_clusters> excluded_clusters[num_excluded_clusters];
  
  // Simulation
  
  // int<lower = 1> num_grid_obs; // Simulation observations
  // vector[num_grid_obs] grid_dist2[num_treatments]; // Simulation distances
  // matrix[num_grid_obs, num_knots_v] Z_grid_v2[num_treatments];
  // vector[num_grid_obs] grid_dist; // Simulation distances
  // matrix[num_grid_obs, num_knots_v] Z_grid_v;
}

transformed data {
  matrix<lower = 0, upper = 1>[num_treatments, num_treatments] treatment_map_design_matrix = rep_matrix(0, num_treatments, num_treatments);
  matrix<lower = 0, upper = 1>[num_treatments, num_treatments] restricted_treatment_map_design_matrix;
  matrix<lower = 0, upper = 1>[num_clusters, num_treatments] cluster_treatment_design_matrix = rep_matrix(0, num_clusters, num_treatments);
  
  int<lower = 1> cluster_size[num_clusters] = count(num_clusters, obs_cluster_id);
  int<lower = 0> cluster_takeup_count[num_clusters] = count_by_group_test(takeup, obs_cluster_id, { 1 }, 1);
  
  int<lower = 1, upper = num_treatments> assigned_treatment[num_obs] = cluster_assigned_treatment[obs_cluster_id]; 
 
  int<lower = 0, upper = num_clusters> num_included_clusters = num_clusters - num_excluded_clusters;
  int<lower = 1, upper = num_clusters> included_clusters[num_included_clusters] = num_excluded_clusters > 0 ? which(seq(1, num_clusters, 1), excluded_clusters, 0) : seq(1, num_clusters, 1);
  int<lower = 0, upper = num_obs> num_included_obs = sum(cluster_size[included_clusters]);
  int<lower = 0, upper = num_obs> num_excluded_obs = num_obs - num_included_obs;
  int<lower = 0, upper = num_obs> included_obs[num_included_obs] = num_included_obs != num_obs ? which(obs_cluster_id, included_clusters, 1) : seq(1, num_obs, 1);
  int<lower = 0, upper = num_obs> excluded_obs[num_excluded_obs]; 
  
  int num_treatments_param_kappa = use_cost_model == COST_MODEL_TYPE_PARAM_KAPPA ? num_treatments : 0;
  int num_treatments_param;
  int num_treatments_param_quadratic = in_array(use_cost_model, { COST_MODEL_TYPE_PARAM_QUADRATIC }) ? num_treatments : 0;
  
  if (in_array(use_cost_model, { COST_MODEL_TYPE_SEMIPARAM, COST_MODEL_TYPE_PARAM_LINEAR, COST_MODEL_TYPE_PARAM_QUADRATIC })) {
    num_treatments_param = num_treatments;
  } else if (use_cost_model == COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE) {
    num_treatments_param = 1;
  }
  
  if (num_excluded_obs > 0) {
    excluded_obs = which(obs_cluster_id, included_clusters, 0);
  }
  
  treatment_map_design_matrix[, 1] = rep_vector(1, num_treatments);

  for(treatment_index in 1:num_treatments) {
    treatment_map_design_matrix[treatment_index, treatment_index] = 1;
  }
  
  restricted_treatment_map_design_matrix = treatment_map_design_matrix;
  
  if (use_private_incentive_restrictions) {
    restricted_treatment_map_design_matrix[3,4] = 1; // Calendar effect is > then bracelet
  }
  
  cluster_treatment_design_matrix = treatment_map_design_matrix[cluster_assigned_treatment];
}

parameters {
  // Levels: control ink calendar bracelet
  real beta_control;
  real beta_ink_effect;
  real <lower = (use_private_incentive_restrictions ? 0 : negative_infinity())> beta_calendar_effect;
  real <lower = (use_private_incentive_restrictions ? 0 : negative_infinity())> beta_bracelet_effect;
  
  matrix[use_cluster_effects ? num_clusters : 0, num_treatments] structural_beta_cluster_raw;
  row_vector<lower = 0>[use_cluster_effects ? num_treatments : 0] structural_beta_cluster_sd;
  
  // Salience
  
  real<lower = 0> beta_salience;
  real<lower = 0> dist_beta_salience;
  
  // V Mixture
  
  vector<lower = 0.5, upper = 1>[num_v_mix - 1] recursive_lambda_v_mix;
  
  vector<lower = 0>[num_v_mix - 1] v_mix_mean_diff;
  vector<lower = 0>[num_v_mix - 1] v_mix_sd;
  
  // Reputational Returns
  
  row_vector<lower = 0>[suppress_reputation ? 0 : num_treatments] mu_rep_raw;
  // real<lower = 0, upper = 1> v_sd;
  vector<lower = 0>[suppress_shocks ? 0 : 2] ub_ur_sd;
  cholesky_factor_corr[suppress_shocks ? 0 : 3] L_all_u_corr;
  
  matrix[use_mu_cluster_effects && !suppress_reputation ? num_clusters : 0, num_treatments] mu_cluster_effects_raw;
  row_vector<lower = 0>[use_mu_cluster_effects && !suppress_reputation ? num_treatments : 0] mu_cluster_effects_sd;
  
  // Parameteric Cost Model
  
  vector<lower = 0>[num_treatments_param_kappa] dist_cost_k_raw;
  
  // Semiparameteric Cost Model
  
  vector[num_treatments_param] dist_beta_v; // Linear distance*treatment effects
  vector[num_treatments_param_quadratic] dist_quadratic_beta_v; // Quadratic distance*treatment effects
  
  matrix[use_cost_model == COST_MODEL_TYPE_SEMIPARAM ? num_treatments_param : 0, num_knots_v] u_splines_v_raw;
  real<lower = 0> u_splines_v_sigma;
  
  // Simulation
  
  // matrix[use_cluster_effects && simulate_new_data ? num_grid_obs : 0, use_private_incentive_restrictions ? num_treatments - 1 : num_treatments] sim_structural_beta_cluster_raw;
  
}

transformed parameters {
  // vector[num_treatments] restricted_structural_beta;
  vector[num_treatments] beta;
  vector[num_treatments] structural_treatment_effect;
  vector[num_clusters] structural_cluster_benefit_cost;
  vector[num_clusters] cluster_effects = rep_vector(0, num_clusters);
  
  row_vector<lower = 0>[suppress_reputation ? 0 : num_treatments] mu_rep = mu_rep_raw;
  matrix<lower = 0>[!suppress_reputation ? num_clusters : 0, num_treatments] cluster_mu_rep;
  
  vector[num_clusters] structural_cluster_obs_v = rep_vector(0, num_clusters);
  matrix<lower = 0, upper = 1>[num_v_mix, num_clusters] structural_cluster_takeup_prob;
  
  simplex[num_v_mix] lambda_v_mix;
  ordered[num_v_mix - 1] v_mix_mean = cumulative_sum(v_mix_mean_diff);
  
  vector<lower = 0>[num_treatments_param_kappa] dist_cost_k;
  vector[num_treatments] linear_dist_cost = rep_vector(0, num_treatments);
  vector[num_treatments] quadratic_dist_cost = rep_vector(0, num_treatments);
  matrix[num_treatments, num_knots_v] u_splines_v = rep_matrix(0, num_treatments, num_knots_v);
  vector[num_clusters] cluster_dist_cost = rep_vector(0, num_clusters);
  
  real total_error_sd;
  
  // Simulation
  // matrix[num_grid_obs, num_treatments] sim_benefit_cost = rep_matrix(0, num_grid_obs, num_treatments); 
  // matrix[num_grid_obs, num_treatments] sim_cluster_effects = rep_matrix(0, num_grid_obs, num_treatments);
  // matrix[num_grid_obs, num_treatments] sim_v = rep_matrix(0, num_grid_obs, num_treatments);
  // matrix<lower = 0, upper = 1>[num_grid_obs, num_treatments] sim_takeup_prob = rep_matrix(0, num_grid_obs, num_treatments);
  
  // if (use_private_incentive_restrictions) {
  //   restricted_structural_beta = append_row(structural_beta, structural_beta[3] + onesided_structural_beta[1]);
  // } else {
  //   restricted_structural_beta = structural_beta;
  // }
  
  beta = [ beta_control, beta_ink_effect, beta_calendar_effect, beta_bracelet_effect ]';
  
  structural_treatment_effect = restricted_treatment_map_design_matrix * beta;
  
  
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
    
    if (use_mu_cluster_effects) {
      matrix[num_clusters, num_treatments] mu_cluster_effects =  mu_cluster_effects_raw .* rep_matrix(mu_cluster_effects_sd, num_clusters);
      
      cluster_mu_rep = rep_matrix(mu_rep, num_clusters) .* exp(mu_cluster_effects);
    } else {
      cluster_mu_rep = rep_matrix(mu_rep, num_clusters);
    }
  }
  
  if (use_salience_effect) {
    structural_treatment_effect += use_salience_effect * mu_rep';
  }
  
  if (use_cost_model == COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE) {
    linear_dist_cost = rep_vector(dist_beta_v[1], num_treatments) + dist_beta_salience * mu_rep';
    
    cluster_dist_cost = param_dist_cost(cluster_standard_dist, 
                                        linear_dist_cost[cluster_assigned_treatment],
                                        quadratic_dist_cost[cluster_assigned_treatment],
                                        u_splines_v[cluster_assigned_treatment],
                                        Z_splines_v[cluster_assigned_treatment]);
  } else if (num_treatments_param_kappa > 0) { 
    dist_cost_k = dist_cost_k_raw;
   
    if (use_cost_k_restrictions) { 
      dist_cost_k[2] += dist_cost_k[1];
      dist_cost_k[3] += dist_cost_k[1];
      dist_cost_k[4] += dist_cost_k[1];
    }
    
    cluster_dist_cost = param_kappa_dist_cost(cluster_standard_dist, dist_cost_k[cluster_assigned_treatment]);
  } else {
    if (num_treatments_param > 0) {
      linear_dist_cost = treatment_map_design_matrix * dist_beta_v;
    } 
    
    if (num_treatments_param_quadratic > 0) {
      quadratic_dist_cost = treatment_map_design_matrix * dist_quadratic_beta_v;
    }
     
    if (use_cost_model == COST_MODEL_TYPE_SEMIPARAM) {
      for (treatment_index in 1:num_treatments_param) {
        u_splines_v[treatment_index] = u_splines_v_raw[treatment_index] * u_splines_v_sigma;
      }        
    }
      
    cluster_dist_cost = param_dist_cost(cluster_standard_dist, 
                                        linear_dist_cost[cluster_assigned_treatment],
                                        quadratic_dist_cost[cluster_assigned_treatment],
                                        u_splines_v[cluster_assigned_treatment],
                                        Z_splines_v[cluster_assigned_treatment]);
  }
  
  structural_cluster_benefit_cost = structural_treatment_effect[cluster_assigned_treatment] - cluster_dist_cost;
  
  if (use_cluster_effects) {
    matrix[num_clusters, num_treatments] structural_beta_cluster = structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_clusters);
    
    cluster_effects = rows_dot_product(cluster_treatment_design_matrix, structural_beta_cluster); 
    structural_cluster_benefit_cost += cluster_effects;
  }

  if (suppress_reputation) {
    structural_cluster_obs_v = - structural_cluster_benefit_cost;
  } else {
    for (cluster_index in 1:num_clusters) {
      structural_cluster_obs_v[cluster_index] = algebra_solver(v_fixedpoint_solution_normal,
                                                               [ - structural_cluster_benefit_cost[cluster_index] ]',
                                                               prepare_solver_theta(structural_cluster_benefit_cost[cluster_index], 
                                                                                    cluster_mu_rep[cluster_index, cluster_assigned_treatment[cluster_index]],
                                                                                    lambda_v_mix,
                                                                                    v_mix_mean,
                                                                                    // suppress_shocks ? 1 : v_sd,
                                                                                    1,
                                                                                    v_mix_sd),
                                                               { 0.0 },
                                                               { num_v_mix },
                                                               1e-10,
                                                               1e-5,
                                                               1e6)[1];
    }
  }
  
  if (!suppress_shocks) {
    vector[3] all_u_sd = append_row(1, ub_ur_sd);
    matrix[3, 3] L_all_u_vcov = diag_pre_multiply(all_u_sd, L_all_u_corr);
    
    total_error_sd = sqrt(sum(L_all_u_vcov * L_all_u_vcov'));
  } else {
    total_error_sd = 1;
  }
  
  for (mix_index in 1:num_v_mix) {
    if (mix_index == 1) {
      structural_cluster_takeup_prob[1] = Phi(- structural_cluster_obs_v / total_error_sd)';
    } else {
      structural_cluster_takeup_prob[mix_index] = Phi(- (structural_cluster_obs_v - v_mix_mean[mix_index - 1]) / v_mix_sd[mix_index - 1])';
    }
  }
  
  /*if (simulate_new_data) { 
    if (use_cluster_effects) {
      matrix[num_grid_obs, num_treatments] sim_structural_beta_cluster;
      
      if (use_private_incentive_restrictions) {
        sim_structural_beta_cluster[, 1:3] = sim_structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_grid_obs);
        sim_structural_beta_cluster[, 4] = sim_structural_beta_cluster[, 3];
      } else {
        sim_structural_beta_cluster = sim_structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_grid_obs);
      }
      
      for (treatment_index in 1:num_treatments) {
        sim_cluster_effects[, treatment_index] = sim_structural_beta_cluster * treatment_map_design_matrix[treatment_index]';
      }
    }
    
    sim_benefit_cost = rep_matrix(structural_treatment_effect', num_grid_obs) + sim_cluster_effects;  
    
    for (treatment_index in 1:num_treatments) {
      if (use_cost_model == COST_MODEL_TYPE_SEMIPARAM) {
        sim_benefit_cost[, treatment_index] -= param_dist_cost(grid_dist2[treatment_index], [ linear_dist_cost[treatment_index] ]', rep_matrix(u_splines_v[treatment_index], num_grid_obs), Z_grid_v2[treatment_index]);
      } else {
        sim_benefit_cost[, treatment_index] -= param_kappa_dist_cost(grid_dist2[treatment_index], [ dist_cost_k[treatment_index] ]');
      }
      
      if (suppress_reputation) {
        sim_v[, treatment_index] = - sim_benefit_cost[, treatment_index];
      } else {
        for (sim_obs_index in 1:num_grid_obs) {
          sim_v[sim_obs_index, treatment_index] = algebra_solver(v_fixedpoint_solution_normal,
                                                                 [ -sim_benefit_cost[sim_obs_index, treatment_index] ]',
                                                                 prepare_solver_theta(sim_benefit_cost[sim_obs_index, treatment_index],
                                                                                      mu_rep[treatment_index],
                                                                                      lambda_v_mix,
                                                                                      v_mix_mean,
                                                                                      suppress_shocks ? 1 : v_sd,
                                                                                      v_mix_sd),
                                                                 { 0.0 },
                                                                 { num_v_mix },
                                                                 1e-10,
                                                                 1e-5,
                                                                 1e6)[1];
        }
      }
      
      sim_takeup_prob[, treatment_index] = Phi(- sim_v[, treatment_index]);
    }
  }*/
}

model {
  beta_control ~ normal(0, 5);
  beta_ink_effect ~ normal(0, 1);
  beta_calendar_effect ~ normal(0, 1);
  beta_bracelet_effect ~ normal(0, 1);
  
  beta_salience ~ normal(0, 1);
  dist_beta_salience ~ normal(0, 1);

  if (!suppress_reputation) { 
    mu_rep_raw ~ normal(0, mu_rep_sd);
    
    if (use_mu_cluster_effects) {
      to_vector(mu_cluster_effects_raw) ~ normal(0, 1);
      mu_cluster_effects_sd ~ normal(0, 1);
    }
  }

  if (num_treatments_param_kappa > 0) {
    dist_cost_k_raw ~ normal(0, 1);
  } 
  
  u_splines_v_sigma ~ normal(0, u_splines_v_sigma_sd);
  
  if (num_treatments_param > 0) {
    dist_beta_v ~ normal(0, 1);
    
    if (use_cost_model == COST_MODEL_TYPE_SEMIPARAM) {
      to_vector(u_splines_v_raw) ~ normal(0, 1);
    }
  }
  
  if (num_treatments_param_quadratic > 0) {
    dist_quadratic_beta_v ~ normal(0, 1);
  }
  
  if (use_cluster_effects) {
    to_vector(structural_beta_cluster_raw) ~ normal(0, 1);
    structural_beta_cluster_sd ~ normal(0, 0.25);
    
    // if (simulate_new_data) {
    //   to_vector(sim_structural_beta_cluster_raw) ~ normal(0, 1);
    // }
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
      cluster_takeup_count[included_clusters] ~ binomial(cluster_size[included_clusters], structural_cluster_takeup_prob[1, included_clusters]);
    } else {
      takeup[included_obs] ~ bernoulli(structural_cluster_takeup_prob[1, obs_cluster_id[included_obs]]);
    }
  } else {
    if (use_binomial) {
      cluster_takeup_count[included_clusters] ~ mixed_binomial(lambda_v_mix, cluster_size[included_clusters], structural_cluster_takeup_prob[, included_clusters]);
    } else {
      takeup[included_obs] ~ mixed_bernoulli(lambda_v_mix, structural_cluster_takeup_prob[, obs_cluster_id[included_obs]]);
    }
  }
}

generated quantities { 
  // Cross Validation
  vector[use_binomial ? num_included_clusters : num_included_obs] log_lik = rep_vector(negative_infinity(), use_binomial ? num_included_clusters : num_included_obs); 
  vector[use_binomial ? num_excluded_clusters : num_excluded_obs] log_lik_heldout = rep_vector(negative_infinity(), use_binomial ? num_excluded_clusters : num_excluded_obs); 
  
  corr_matrix[3] all_u_corr = L_all_u_corr * L_all_u_corr';
  
  // Simulation
  // matrix[num_grid_obs, num_treatments] sim_benefit_cost = rep_matrix((treatment_map_design_matrix * restricted_structural_beta)', num_grid_obs); 
  // matrix[num_grid_obs, num_treatments] sim_v;
  // matrix<lower = 0, upper = 1>[num_grid_obs, num_treatments] sim_takeup_prob;
  // matrix[num_grid_obs, num_treatments] sim_cluster_effect = rep_matrix(0, num_grid_obs, num_treatments);
  
  // Cross Validation 
  
  if (use_binomial) {
    for (cluster_index_index in 1:num_included_clusters) {
      int cluster_index = included_clusters[cluster_index_index];
      
      if (num_v_mix == 1) {
        log_lik[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], structural_cluster_takeup_prob[cluster_index]);
      } else {
        log_lik[cluster_index_index] = mixed_binomial_lpmf({ cluster_takeup_count[cluster_index] } | lambda_v_mix, { cluster_size[cluster_index] }, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
    }
    
    for (cluster_index_index in 1:num_excluded_clusters) {
      int cluster_index = excluded_clusters[cluster_index_index];
      
      if (num_v_mix == 1) {
        log_lik_heldout[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], structural_cluster_takeup_prob[cluster_index]);
      } else {
        log_lik_heldout[cluster_index_index] = mixed_binomial_lpmf({ cluster_takeup_count[cluster_index] } | lambda_v_mix, { cluster_size[cluster_index] }, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
    }
  } else {
    for (obs_index_index in 1:num_included_obs) {
      int obs_index = included_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      if (num_v_mix == 1) {
        log_lik[obs_index_index] = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[1, cluster_index:cluster_index]);
      } else {
        log_lik[obs_index_index] = mixed_bernoulli_lpmf({ takeup[obs_index] } | lambda_v_mix, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
    }
    
    for (obs_index_index in 1:num_excluded_obs) {
      int obs_index = excluded_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      if (num_v_mix == 1) {
        log_lik_heldout[obs_index_index] = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[1, cluster_index:cluster_index]);
      } else {
        log_lik_heldout[obs_index_index] = mixed_bernoulli_lpmf({ takeup[obs_index] } | lambda_v_mix, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
    }
  }
 
  // Simulation 
  
  // for (treatment_index in 1:num_treatments) {
  //   if (use_semiparam_cost) {
  //     sim_benefit_cost[, treatment_index] -= param_dist_cost(grid_dist2[treatment_index], [ dist_beta_v[treatment_index] ]', rep_matrix(u_splines_v[treatment_index], num_grid_obs), Z_grid_v2[treatment_index]);
  //   } else {
  //     sim_benefit_cost[, treatment_index] -= param_kappa_dist_cost(grid_dist2[treatment_index], [ dist_cost_k[treatment_index] ]');
  //   }
  //   
  //   for (sim_obs_index in 1:num_grid_obs) {
  //     // sim_cluster_effect[sim_obs_index, treatment_index] = normal_rng(0, structural_beta_cluster_sd[min(treatment_index, 3)]);
  //     // sim_benefit_cost[sim_obs_index, treatment_index] += sim_cluster_effect[sim_obs_index, treatment_index];
  //     
  //     sim_v[sim_obs_index, treatment_index] = algebra_solver(v_fixedpoint_solution_normal,
  //                                                            [ -sim_benefit_cost[sim_obs_index, treatment_index] ]',
  //                                                            prepare_solver_theta(sim_benefit_cost[sim_obs_index, treatment_index],
  //                                                                                 mu_rep[treatment_index],
  //                                                                                 lambda_v_mix,
  //                                                                                 v_mix_mean,
  //                                                                                 v_mix_sd),
  //                                                            { 0.0 },
  //                                                            { num_v_mix },
  //                                                            1e-10,
  //                                                            1e-4,
  //                                                            1e10)[1];
  //                                                              
  //     sim_takeup_prob[sim_obs_index, treatment_index] = Phi(- sim_v[sim_obs_index, treatment_index]);
  //   }
  // }
}
