functions {
#include /../multilvlr/util.stan
  
  real reputational_returns_normal(real v, vector lambda, vector mix_mean, vector mix_sd) {
    int num_mix = num_elements(lambda);
    real Phi_v = Phi(v);
    real rep = exp(std_normal_lpdf(v)) / (Phi_v * (1 - Phi_v));
     
    if (num_mix > 1) { 
      rep = lambda[1] * rep;
      
      for (mix_index in 2:num_mix) {
        real mix_Phi_v = Phi((v - mix_mean[mix_index - 1]) / mix_sd[mix_index - 1]);
        
        rep += lambda[mix_index] * exp(normal_lpdf(v | mix_mean[mix_index - 1], mix_sd[mix_index - 1])) / (mix_Phi_v * (1 - mix_Phi_v));
      }
    }
    
    return rep;
  }
  
  vector param_dist_cost(vector dist, vector k) {
    if (num_elements(k) == 1) {
      return (k[1] * square(dist)) / 2;
    } else {
      return (k .* square(dist)) / 2;
    }
  }
  
  vector semiparam_dist_cost(vector dist, vector linear_dist_cost, matrix u_splines, matrix Z_splines) {
    if (num_elements(linear_dist_cost)) {
      return dist * linear_dist_cost[1] + rows_dot_product(u_splines, Z_splines);
    } else {
      return dist .* linear_dist_cost + rows_dot_product(u_splines, Z_splines);
    }
  }
  
  vector v_fixedpoint_solution_normal(vector model_param, vector theta, real[] x_r, int[] x_i) {
    real v_cutoff = model_param[1];
    
    real benefit_cost = theta[1];
    real mu = theta[2];
    
    int num_v_mix = x_i[1];
    
    vector[num_v_mix] lambda = theta[3:(3 + num_v_mix - 1)];
    vector[num_v_mix - 1] mix_mean = theta[(3 + num_v_mix):(3 + num_v_mix + num_v_mix - 2)];
    vector[num_v_mix - 1] mix_sd = theta[(3 + num_v_mix + num_v_mix - 1):(3 + 3 * num_v_mix - 3)];
    
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
  
  vector prepare_solver_theta(real benefit_cost, real mu_rep, vector lambda, vector mix_mean, vector mix_sd) {
    int num_v_mix = num_elements(lambda);
    vector[2 + num_v_mix + 2 * (num_v_mix - 1)] solver_theta;
    
    solver_theta[1:2] = [ benefit_cost, mu_rep ]';
    solver_theta[3:(3 + num_v_mix - 1)] = lambda;
    solver_theta[(3 + num_v_mix):(3 + num_v_mix + num_v_mix - 2)] = mix_mean;
    solver_theta[(3 + num_v_mix + num_v_mix - 1):(3 + 3 * num_v_mix - 3)] = mix_sd;
 
    return solver_theta; 
  }
}

data {
  int<lower = 1> num_v_mix;
  
  int<lower = 0, upper = 1> use_binomial;
  int<lower = 0, upper = 1> use_cluster_effects;
  int<lower = 0, upper = 1> use_cost_k_restrictions;
  int<lower = 0, upper = 1> use_private_incentive_restrictions;
  int<lower = 0, upper = 1> use_semiparam_cost;
  int<lower = 0, upper = 1> suppress_reputation;
  int<lower = 0, upper = 1> simulate_new_data; 

  int<lower = 1> num_obs; // Actual observations
  int<lower = 1> num_clusters;
  int<lower = 1> num_treatments;
  
  int<lower = 1, upper = num_clusters> obs_cluster_id[num_obs];
  
  int<lower = 0, upper = 1> takeup[num_obs]; // Observed outcome variable
  vector[num_clusters] cluster_standard_dist; // Standardized distance to treatment
  
  int<lower = 1, upper = num_treatments> cluster_assigned_treatment[num_clusters]; // Actual assigned treatments 
  
  // Semiparametric Cost Model (Splines)
  
  int<lower = 3> num_knots_v;
  matrix[num_clusters, num_knots_v] Z_splines_v; 
  
  real<lower = 0> u_splines_v_sigma_sd;
  
  // K-fold CV 
  
  int<lower = 0> num_excluded_clusters;
  int<lower = 1, upper = num_clusters> excluded_clusters[num_excluded_clusters];
  
  // Simulation
  
  int<lower = 1> num_grid_obs; // Simulation observations
  vector[num_grid_obs] grid_dist2[num_treatments]; // Simulation distances
  matrix[num_grid_obs, num_knots_v] Z_grid_v2[num_treatments];
  // vector[num_grid_obs] grid_dist; // Simulation distances
  // matrix[num_grid_obs, num_knots_v] Z_grid_v;
}

transformed data {
  matrix<lower = 0, upper = 1>[num_treatments, num_treatments] treatment_map_design_matrix = rep_matrix(0, num_treatments, num_treatments);
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
  
  int num_treatments_param = use_semiparam_cost ? 0 : num_treatments;
  int num_treatments_semiparam = use_semiparam_cost ? num_treatments : 0;
  
  if (num_excluded_obs > 0) {
    excluded_obs = which(obs_cluster_id, included_clusters, 0);
  }
  
  treatment_map_design_matrix[, 1] = rep_vector(1, num_treatments);

  for(treatment_index in 1:num_treatments) {
    treatment_map_design_matrix[treatment_index, treatment_index] = 1;
  }
  
  cluster_treatment_design_matrix = treatment_map_design_matrix[cluster_assigned_treatment];
}

parameters {
  // Levels: control ink calendar bracelet
  vector[use_private_incentive_restrictions ? num_treatments - 1 : num_treatments] structural_beta;
  matrix[use_cluster_effects ? num_clusters : 0, num_treatments] structural_beta_cluster_raw;
  // matrix[use_cluster_effects ? num_clusters : 0, use_private_incentive_restrictions ? num_treatments - 1 : num_treatments] structural_beta_cluster_raw;
  row_vector<lower = 0>[use_cluster_effects ? num_treatments : 0] structural_beta_cluster_sd;
  // row_vector<lower = 0>[use_cluster_effects ? (use_private_incentive_restrictions ? num_treatments - 1 : num_treatments) : 0] structural_beta_cluster_sd;
  
  vector<lower = 0>[suppress_reputation ? 0 : num_treatments] mu_rep_raw;
  
  // Parameteric Cost Model
  
  vector<lower = 0>[num_treatments_param] dist_cost_k_raw;
  
  // Semiparameteric Cost Model
  
  vector[num_treatments_semiparam] dist_beta_v; // Linear distance*treatment effects
  matrix[num_treatments_semiparam, num_knots_v] u_splines_v_raw;
  real<lower = 0> u_splines_v_sigma;
  
  // Simulation
  
  matrix[use_cluster_effects && simulate_new_data ? num_grid_obs : 0, use_private_incentive_restrictions ? num_treatments - 1 : num_treatments] sim_structural_beta_cluster_raw;
  
  // V Mixture
  
  vector<lower = 0.5, upper = 1>[num_v_mix - 1] recursive_lambda_v_mix;
  
  vector<lower = 0>[num_v_mix - 1] v_mix_mean_diff;
  vector<lower = 0>[num_v_mix - 1] v_mix_sd;
}

transformed parameters {
  vector[num_treatments] restricted_structural_beta;
  vector[num_treatments] structural_treatment_effect;
  vector[num_clusters] structural_cluster_benefit_cost;
  vector[num_clusters] cluster_effects = rep_vector(0, num_clusters);
  vector<lower = 0>[suppress_reputation ? 0 : num_treatments] mu_rep = mu_rep_raw;
  vector<lower = 0>[num_treatments_param] dist_cost_k;
  vector[num_clusters] structural_cluster_obs_v = rep_vector(0, num_clusters);
  matrix<lower = 0, upper = 1>[num_v_mix, num_clusters] structural_cluster_takeup_prob;
  
  simplex[num_v_mix] lambda_v_mix;
  ordered[num_v_mix - 1] v_mix_mean = cumulative_sum(v_mix_mean_diff);
  
  vector[num_treatments_semiparam] linear_dist_cost;
  matrix[num_treatments_semiparam, num_knots_v] u_splines_v;
  
  // Simulation
  matrix[num_grid_obs, num_treatments] sim_benefit_cost = rep_matrix(0, num_grid_obs, num_treatments); 
  matrix[num_grid_obs, num_treatments] sim_cluster_effects = rep_matrix(0, num_grid_obs, num_treatments);
  matrix[num_grid_obs, num_treatments] sim_v = rep_matrix(0, num_grid_obs, num_treatments);
  matrix<lower = 0, upper = 1>[num_grid_obs, num_treatments] sim_takeup_prob = rep_matrix(0, num_grid_obs, num_treatments);
  
  if (use_private_incentive_restrictions) {
    restricted_structural_beta = append_row(structural_beta, structural_beta[3]);
  } else {
    restricted_structural_beta = structural_beta;
  }
  
  structural_treatment_effect = treatment_map_design_matrix * restricted_structural_beta;
  
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
  }
  
  if (num_treatments_param > 0) { 
    dist_cost_k = dist_cost_k_raw;
   
    if (use_cost_k_restrictions) { 
      dist_cost_k[2] += dist_cost_k[1];
      dist_cost_k[3] += dist_cost_k[1];
      dist_cost_k[4] += dist_cost_k[1];
    }
    
    structural_cluster_benefit_cost = structural_treatment_effect[cluster_assigned_treatment] - param_dist_cost(cluster_standard_dist, dist_cost_k[cluster_assigned_treatment]);
  }
  
  if (num_treatments_semiparam > 0) {
    linear_dist_cost = treatment_map_design_matrix * dist_beta_v;
    
    for (treatment_index in 1:num_treatments_semiparam) {
      u_splines_v[treatment_index] = u_splines_v_raw[treatment_index] * u_splines_v_sigma;
    }
    
    structural_cluster_benefit_cost = structural_treatment_effect[cluster_assigned_treatment] - semiparam_dist_cost(cluster_standard_dist, 
                                                                                                                    linear_dist_cost[cluster_assigned_treatment],
                                                                                                                    u_splines_v[cluster_assigned_treatment],
                                                                                                                    Z_splines_v[cluster_assigned_treatment]);
  }
  
  if (use_cluster_effects) {
    matrix[num_clusters, num_treatments] structural_beta_cluster;
    
    if (use_private_incentive_restrictions) {
      structural_beta_cluster = structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_clusters);
      // structural_beta_cluster[, 1:3] = structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_clusters);
      // structural_beta_cluster[, 4] = structural_beta_cluster[, 3];
    } else {
      structural_beta_cluster = structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_clusters);
    }
    
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
                                                                                    mu_rep[cluster_assigned_treatment[cluster_index]],
                                                                                    lambda_v_mix,
                                                                                    v_mix_mean,
                                                                                    v_mix_sd),
                                                               { 0.0 },
                                                               { num_v_mix },
                                                               1e-10,
                                                               1e-5,
                                                               1e6)[1];
    }
  }
  
  for (mix_index in 1:num_v_mix) {
    if (mix_index == 1) {
      structural_cluster_takeup_prob[1] = Phi(- structural_cluster_obs_v)';
    } else {
      structural_cluster_takeup_prob[mix_index] = Phi(- (structural_cluster_obs_v - v_mix_mean[mix_index - 1]) / v_mix_sd[mix_index - 1])';
    }
  }
  
  if (simulate_new_data) { 
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
      if (use_semiparam_cost) {
        sim_benefit_cost[, treatment_index] -= semiparam_dist_cost(grid_dist2[treatment_index], [ linear_dist_cost[treatment_index] ]', rep_matrix(u_splines_v[treatment_index], num_grid_obs), Z_grid_v2[treatment_index]);
      } else {
        sim_benefit_cost[, treatment_index] -= param_dist_cost(grid_dist2[treatment_index], [ dist_cost_k[treatment_index] ]');
      }
      
      if (suppress_reputation) {
        sim_v[, treatment_index] = - sim_benefit_cost[, treatment_index];
      } else {
        for (sim_obs_index in 1:num_grid_obs) {
          // sim_cluster_effect[sim_obs_index, treatment_index] = normal_rng(0, structural_beta_cluster_sd[min(treatment_index, 3)]);
          // sim_benefit_cost[sim_obs_index, treatment_index] += sim_cluster_effect[sim_obs_index, treatment_index];
          
          sim_v[sim_obs_index, treatment_index] = algebra_solver(v_fixedpoint_solution_normal,
                                                                 [ -sim_benefit_cost[sim_obs_index, treatment_index] ]',
                                                                 prepare_solver_theta(sim_benefit_cost[sim_obs_index, treatment_index],
                                                                                      mu_rep[treatment_index],
                                                                                      lambda_v_mix,
                                                                                      v_mix_mean,
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
  }
}

model {
  u_splines_v_sigma ~ normal(0, u_splines_v_sigma_sd);
  
  structural_beta ~ normal(0, 1);

  if (!suppress_reputation) { 
    mu_rep_raw ~ normal(0, 0.1);
  }

  if (num_treatments_param > 0) {
    dist_cost_k_raw ~ normal(0, 1);
  } 
  
  if (num_treatments_semiparam > 0) {
    dist_beta_v ~ normal(0, 1);
    u_splines_v_sigma ~ normal(0, u_splines_v_sigma_sd);
    to_vector(u_splines_v_raw) ~ normal(0, 1);
  }
  
  if (use_cluster_effects) {
    to_vector(structural_beta_cluster_raw) ~ normal(0, 1);
    structural_beta_cluster_sd ~ normal(0, 0.25);
    
    if (simulate_new_data) {
      to_vector(sim_structural_beta_cluster_raw) ~ normal(0, 1);
    }
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
  // Cross validation
  vector[use_binomial ? num_included_clusters : num_included_obs] log_lik = rep_vector(negative_infinity(), use_binomial ? num_included_clusters : num_included_obs); 
  vector[use_binomial ? num_excluded_clusters : num_excluded_obs] log_lik_heldout = rep_vector(negative_infinity(), use_binomial ? num_excluded_clusters : num_excluded_obs); 
  
  // Simulation
  // matrix[num_grid_obs, num_treatments] sim_benefit_cost = rep_matrix((treatment_map_design_matrix * restricted_structural_beta)', num_grid_obs); 
  // matrix[num_grid_obs, num_treatments] sim_v;
  // matrix<lower = 0, upper = 1>[num_grid_obs, num_treatments] sim_takeup_prob;
  // matrix[num_grid_obs, num_treatments] sim_cluster_effect = rep_matrix(0, num_grid_obs, num_treatments);
  
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
  
  // for (treatment_index in 1:num_treatments) {
  //   if (use_semiparam_cost) {
  //     sim_benefit_cost[, treatment_index] -= semiparam_dist_cost(grid_dist2[treatment_index], [ dist_beta_v[treatment_index] ]', rep_matrix(u_splines_v[treatment_index], num_grid_obs), Z_grid_v2[treatment_index]);
  //   } else {
  //     sim_benefit_cost[, treatment_index] -= param_dist_cost(grid_dist2[treatment_index], [ dist_cost_k[treatment_index] ]');
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
