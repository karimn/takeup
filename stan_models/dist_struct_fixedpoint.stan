functions {
#include /../multilvlr/util.stan
  
  real reputational_returns_normal(real v, vector lambda, vector mix_mean) {
    int num_mix = num_elements(lambda);
    real Phi_v = Phi(v);
    real rep = exp(std_normal_lpdf(v)) / (Phi_v * (1 - Phi_v));
     
    if (num_mix > 1) { 
      rep = lambda[1] * rep;
      
      for (mix_index in 2:num_mix) {
        real mix_Phi_v = Phi(v - mix_mean[mix_index - 1]);
        
        rep += lambda[mix_index] * exp(normal_lpdf(v | mix_mean[mix_index - 1], 1)) / (mix_Phi_v * (1 - mix_Phi_v));
      }
    }
    
    return rep;
  }
  
  // vector reputation_returns_normal_vec(vector v) {
  //   int num_v = num_elements(v);
  //   vector[num_v] rep;
  //   
  //   for (v_index in 1:num_v) {
  //     rep[v_index] = reputational_returns_normal(v[v_index]);
  //   }
  //   
  //   return rep;
  // }
  
  vector param_dist_cost(vector dist, vector k) {
  // vector param_dist_cost(vector dist, real k) {
    // return (k * square(dist)) / 2;
    return (k .* square(dist)) / 2;
  }
  
  vector v_fixedpoint_solution_normal(vector model_param, vector theta, real[] x_r, int[] x_i) {
    real v_cutoff = model_param[1];
    
    real benefit_cost = theta[1];
    real mu = theta[2];
    
    int num_v_mix = x_i[1];
    
    vector[num_v_mix] lambda = theta[3:(3 + num_v_mix - 1)];
    vector[num_v_mix - 1] mix_mean = theta[(3 + num_v_mix):(3 + num_v_mix + num_v_mix - 2)];
    
    return [ v_cutoff + benefit_cost + mu * reputational_returns_normal(v_cutoff, lambda, mix_mean) ]';
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
}

data {
  int v_dist_type;
  int<lower = 1> num_v_mix;
  
  int<lower = 0, upper = 1> use_binomial;
  int<lower = 0, upper = 1> use_cluster_effects;
  int<lower = 0, upper = 1> use_cost_k_restrictions;

  int<lower = 1> num_obs; // Actual observations
  int<lower = 1> num_clusters;
  int<lower = 1> num_treatments;
  
  int<lower = 1, upper = num_clusters> obs_cluster_id[num_obs];
  
  int<lower = 0, upper = 1> takeup[num_obs]; // Observed outcome variable
  vector[num_clusters] cluster_standard_dist; // Standardized distance to treatment
  
  int<lower = 1, upper = num_treatments> cluster_assigned_treatment[num_clusters]; // Actual assigned treatments 
  
  // K-fold CV 
  
  int<lower = 0> num_excluded_clusters;
  int<lower = 1, upper = num_clusters> excluded_clusters[num_excluded_clusters];
}

transformed data {
  int V_DIST_TYPE_NORMAL = 1;
  
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
  vector[num_treatments] structural_beta;
  matrix[use_cluster_effects ? num_clusters : 0, num_treatments] structural_beta_cluster_raw;
  row_vector<lower = 0>[use_cluster_effects ? num_treatments : 0] structural_beta_cluster_sd;
  // 
  vector<lower = 0>[num_treatments] mu_rep_raw;
  // real<lower = 0> dist_cost_k;
  vector<lower = 0>[num_treatments] dist_cost_k_raw;
  
  // V Mixture
  
  simplex[num_v_mix] lambda_v_mix;
  ordered[num_v_mix - 1] v_mix_mean;
}

transformed parameters {
  vector[num_treatments] structural_treatment_effect = treatment_map_design_matrix * structural_beta;
  vector[num_clusters] structural_cluster_benefit_cost;
  vector<lower = 0>[num_treatments] mu_rep = mu_rep_raw;
  vector<lower = 0>[num_treatments] dist_cost_k = dist_cost_k_raw;
  vector[num_clusters] structural_cluster_obs_v = rep_vector(0, num_clusters);
  matrix<lower = 0, upper = 1>[num_v_mix, num_clusters] structural_cluster_takeup_prob;

  // Levels: control ink calendar bracelet
  
  mu_rep[2] += mu_rep[1];
  mu_rep[3] += mu_rep[1];
  mu_rep[4] += mu_rep[1];
  
  if (use_cost_k_restrictions) {
    dist_cost_k[2] += dist_cost_k[1];
    dist_cost_k[3] += dist_cost_k[1];
    dist_cost_k[4] += dist_cost_k[1];
  }
  
  structural_cluster_benefit_cost = structural_treatment_effect[cluster_assigned_treatment] - param_dist_cost(cluster_standard_dist, dist_cost_k[cluster_assigned_treatment]);
  
  {
    vector[num_clusters] cluster_effect = rep_vector(0, num_clusters);
    
    if (use_cluster_effects) {
      matrix[num_clusters, num_treatments] structural_beta_cluster = structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_clusters);
      
      cluster_effect = rows_dot_product(cluster_treatment_design_matrix, structural_beta_cluster);
      structural_cluster_benefit_cost += cluster_effect;
    }
  
    for (cluster_index in 1:num_clusters) {
      // print(cluster_index, ": ", structural_cluster_benefit_cost[cluster_index]);
     
      if (v_dist_type == V_DIST_TYPE_NORMAL) { 
        vector[2 + num_v_mix + (num_v_mix - 1)] solver_theta;
        
        solver_theta[1:2] = [ structural_cluster_benefit_cost[cluster_index], mu_rep[cluster_assigned_treatment[cluster_index]] ]';
        solver_theta[3:(3 + num_v_mix - 1)] = lambda_v_mix;
        solver_theta[(3 + num_v_mix):(3 + num_v_mix + num_v_mix - 2)] = v_mix_mean;
        
        // structural_cluster_obs_v[cluster_index] = 0;
        
        structural_cluster_obs_v[cluster_index] = algebra_solver(v_fixedpoint_solution_normal,
                                                                 [ -structural_cluster_benefit_cost[cluster_index] ]',
                                                                 solver_theta,
                                                                 { 0.0 },
                                                                 { num_v_mix },
                                                                 1e-10,
                                                                 1e-5,
                                                                 1e6)[1];
      }
    }
  }
  
  if (v_dist_type == V_DIST_TYPE_NORMAL) { 
    for (mix_index in 1:num_v_mix) {
      if (mix_index == 1) {
        structural_cluster_takeup_prob[1] = Phi(- structural_cluster_obs_v)';
      } else {
        structural_cluster_takeup_prob[mix_index] = Phi(- (structural_cluster_obs_v - v_mix_mean[mix_index - 1]))';
      }
    }
  }
}

model {
  // print("structural_beta = ", structural_beta);
  structural_beta ~ normal(0, 1);
 
  // print("mu_rep_raw = ", mu_rep_raw);
  mu_rep_raw ~ normal(0, 0.1);
 
  // print("dist_cost_k = ", dist_cost_k_raw);
  dist_cost_k_raw ~ normal(0, 1);
  
  if (use_cluster_effects) {
    to_vector(structural_beta_cluster_raw) ~ normal(0, 1);
    structural_beta_cluster_sd ~ normal(0, 0.25);
  }
  
  if (num_v_mix > 1) {
    v_mix_mean ~ normal(0, 1);
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
  vector[use_binomial ? num_included_clusters : num_included_obs] log_lik; 
  vector[use_binomial ? num_excluded_clusters : num_excluded_obs] log_lik_heldout; 
  
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
}
