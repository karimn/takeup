functions {
#include /../multilvlr/util.stan
  
  real reputational_returns_normal(real v) {
    real Phi_v = Phi(v);
      
    return exp(std_normal_lpdf(v)) / (Phi_v * (1 - Phi_v));
  }
  
  vector reputation_returns_normal_vec(vector v) {
    int num_v = num_elements(v);
    vector[num_v] rep;
    
    for (v_index in 1:num_v) {
      rep[v_index] = reputational_returns_normal(v[v_index]);
    }
    
    return rep;
  }
  
  vector param_dist_cost(vector dist, vector k) {
    return (k .* square(dist)) / 2;
  }
  
  // vector cutoff_solution(vector v_cutoff, vector theta, vector x_r, int[] x_i) {
  //   int num_clusters = x_i[1];
  //   int num_treatments = x_i[2];
  //   
  //   real dist_cost_k = theta[1];
  //   vector[num_treatments] benefit_param = theta[2:(2 + num_treatments - 1)];
  //   vector[num_treatments] rep_mu = theta[(num_treatments + 1):(num_treatments + num_treatments)];
  //   
  //   vector[num_clusters] dist = x_r[1:num_clusters];
  //   matrix[num_clusters, num_treatments] cluster_treatment = to_matrix(x_r[(num_clusters + 1):(num_clusters + (num_clusters * num_treatments))]);
  //   
  //   vector[num_clusters] dist_cost = param_dist_cost(dist, dist_cost_k);
  //   
  //   return v_cutoff - dist_cost + cluster_treatment * benefit_param + (cluster_treatment * rep_mu) .* reputation_returns_normal_vec(v_cutoff);
  // }
  
  vector cutoff_solution(vector v_cutoff, vector theta, real[] x_r, int[] x_i) {
    int num_clusters = x_i[1];
    
    real dist_cost_k = theta[1];
    real benefit = theta[2];
    real mu = theta[3];
    
    vector[num_clusters] dist = to_vector(x_r);
    vector[num_clusters] dist_cost = param_dist_cost(dist, rep_vector(dist_cost_k, num_clusters));
    
    return v_cutoff - dist_cost + benefit + mu * reputation_returns_normal_vec(v_cutoff);
  }
  
  vector v_fixedpoint_solution(vector model_param, vector theta, real[] x_r, int[] x_i) {
    real mu = model_param[1];
    
    int num_clusters = x_i[1];
    
    real dist_cost_k = theta[1];
    real benefit = theta[2];
    vector[num_clusters] v_cutoff = theta[3:(3 + num_clusters - 1)];
    
    vector[num_clusters] dist = to_vector(x_r);
    vector[num_clusters] dist_cost = param_dist_cost(dist, rep_vector(dist_cost_k, num_clusters));
    
    return v_cutoff - dist_cost + benefit + mu * reputation_returns_normal_vec(v_cutoff);
  }
}

data {
  // Constants for selecting model to fit
  int MODEL_TYPE_NO_DIST;
  int MODEL_TYPE_DISCRETE_DIST;
  int MODEL_TYPE_LINEAR_DIST;
  int MODEL_TYPE_QUADRATIC_NONLINEAR_DIST;
  int MODEL_TYPE_CUBIC_NONLINEAR_DIST;
  int MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_OSULLIVAN;
  int MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_BSPLINE; // https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html
  int MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_ISPLINE;
  
  int model_type; // Model to fit
  int<lower = 0, upper = 1> use_binomial;
  int<lower = 0, upper = 1> use_cluster_effects;
  int<lower = 0, upper = 1> use_dist_cluster_effects;

    /* 
       0: not structural 
       1: structural with semiparam v* used to calculate structural v*
       2: Fixed point estimation (not working yet)
    */
  int<lower = 0, upper = 2> is_structural;
  
  int<lower = 1> num_obs; // Actual observations
  int<lower = 1> num_clusters;
  int<lower = 1> num_grid_obs; // Simulation observations
  int<lower = 1> num_treatments;
  int<lower = 1> num_dist_group_treatments; // Number of treatments for a dist*treatment fully saturated model (for discrete distance grouping)
  
  int<lower = 1, upper = num_clusters> obs_cluster_id[num_obs];
  
  int<lower = 0, upper = 1> takeup[num_obs]; // Observed outcome variable
  vector[num_clusters] cluster_standard_dist; // Standardized distance to treatment
 
  // Actual assigned treatments 
  int<lower = 1, upper = num_treatments> cluster_assigned_treatment[num_clusters]; 
  int<lower = 1, upper = num_dist_group_treatments> cluster_assigned_dist_group_treatment[num_clusters];
  
  vector[num_grid_obs] grid_dist; // Simulation distances
 
  // Splines 
  int<lower = 3> num_knots_v;
  matrix[num_clusters, num_knots_v] Z_splines_v; 
  matrix[num_grid_obs, num_knots_v] Z_grid_v;
  
  // Hyper parameters
  
  real<lower = 0> u_splines_v_sigma_sd;
  
  // K-fold CV 
  
  int<lower = 0> num_excluded_clusters;
  int<lower = 1, upper = num_clusters> excluded_clusters[num_excluded_clusters];
}

transformed data {
  // These are used to turn on/off some parameters depending on model
  int num_actual_treatments = model_type == MODEL_TYPE_DISCRETE_DIST ? num_dist_group_treatments : num_treatments;
  
  int num_treatments_linear = in_array(model_type, { MODEL_TYPE_LINEAR_DIST,
                                                     MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                                                     MODEL_TYPE_CUBIC_NONLINEAR_DIST,
                                                     MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_OSULLIVAN, 
                                                     MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_BSPLINE, 
                                                     MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_ISPLINE
                                                     }) ? num_treatments : 0;
  
  int num_treatments_quadratic = in_array(model_type, { MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                                                        MODEL_TYPE_CUBIC_NONLINEAR_DIST }) ? num_treatments : 0;
                                                        
  int num_treatments_cubic = model_type == MODEL_TYPE_CUBIC_NONLINEAR_DIST ? num_treatments : 0;
  
  int num_treatments_semiparam = in_array(model_type, { MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_OSULLIVAN, 
                                                        MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_BSPLINE, 
                                                        MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_ISPLINE }) ? num_treatments : 0;
 
  // A lower bound for splines parameters. Normally there is no bound, but for a i-splines model I need the lower to be 0. 
  real u_splines_lb = model_type == MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_ISPLINE ? 0 : negative_infinity();
  
  matrix<lower = 0, upper = 1>[num_actual_treatments, num_actual_treatments] treatment_map_design_matrix = rep_matrix(0, num_actual_treatments, num_actual_treatments);
  matrix<lower = 0, upper = 1>[num_clusters, num_actual_treatments] cluster_treatment_design_matrix = rep_matrix(0, num_clusters, num_actual_treatments);
  
  vector[num_clusters] cluster_standard_dist_square = square(cluster_standard_dist);
  vector[num_clusters] cluster_standard_dist_cube; 
  
  vector[num_grid_obs] grid_dist_square = square(grid_dist);
  vector[num_grid_obs] grid_dist_cube; 
  
  int<lower = 1> cluster_size[num_clusters] = count(num_clusters, obs_cluster_id);
  int<lower = 0> cluster_takeup_count[num_clusters] = count_by_group_test(takeup, obs_cluster_id, { 1 }, 1);
  
  int<lower = 1, upper = num_treatments> assigned_treatment[num_obs] = cluster_assigned_treatment[obs_cluster_id]; 
  int<lower = 1, upper = num_dist_group_treatments> assigned_dist_group_treatment[num_obs] = cluster_assigned_dist_group_treatment[obs_cluster_id];
  
  real<lower = 0> dist_param_sd = 1;
 
  int<lower = 0, upper = num_clusters> num_included_clusters = num_clusters - num_excluded_clusters;
  int<lower = 1, upper = num_clusters> included_clusters[num_included_clusters] = num_excluded_clusters > 0 ? which(seq(1, num_clusters, 1), excluded_clusters, 0) : seq(1, num_clusters, 1);
  int<lower = 0, upper = num_obs> num_included_obs = sum(cluster_size[included_clusters]);
  int<lower = 0, upper = num_obs> num_excluded_obs = num_obs - num_included_obs;
  int<lower = 0, upper = num_obs> included_obs[num_included_obs] = num_included_obs != num_obs ? which(obs_cluster_id, included_clusters, 1) : seq(1, num_obs, 1);
  int<lower = 0, upper = num_obs> excluded_obs[num_excluded_obs]; 
  
  int<lower = 1, upper = num_clusters> treatment_sorted_clusters[num_clusters] = sort_indices_asc(cluster_assigned_treatment);
  int<lower = 0, upper = num_clusters> treatment_cluster_size[num_treatments] = count(num_treatments, cluster_assigned_treatment);
  
  if (num_excluded_obs > 0) {
    excluded_obs = which(obs_cluster_id, included_clusters, 0);
  }
  
  treatment_map_design_matrix[, 1] = rep_vector(1, num_actual_treatments);

  for(treatment_index in 1:num_actual_treatments) {
    treatment_map_design_matrix[treatment_index, treatment_index] = 1;
  }
  
  if (model_type == MODEL_TYPE_DISCRETE_DIST) { 
    cluster_treatment_design_matrix = treatment_map_design_matrix[cluster_assigned_dist_group_treatment];
  } else {
    cluster_treatment_design_matrix = treatment_map_design_matrix[cluster_assigned_treatment];
  }
  
  for (cluster_index in 1:num_clusters) {
    cluster_standard_dist_cube[cluster_index] = pow(cluster_standard_dist[cluster_index], 3);
  }
  
  for (obs_index in 1:num_grid_obs) {
    grid_dist_cube[obs_index] = pow(grid_dist[obs_index], 3);
  }
}

parameters {
  vector[is_structural < 2 ? num_actual_treatments : 0] beta; // Threatment intercepts
  vector[is_structural ? num_actual_treatments : 0] structural_beta_raw; 
  
  matrix[use_cluster_effects && is_structural < 2 ? num_clusters : 0, num_actual_treatments] beta_cluster_raw;
  row_vector<lower = 0>[use_cluster_effects && is_structural < 2 ? num_actual_treatments : 0] beta_cluster_sd;
  
  matrix[use_cluster_effects && is_structural ? num_clusters : 0, num_treatments] structural_beta_cluster_raw;
  row_vector<lower = 0>[use_cluster_effects && is_structural ? num_treatments : 0] structural_beta_cluster_sd;
  
  vector[num_treatments_linear] dist_beta_v; // Linear distance*treatment effects
  vector[num_treatments_quadratic] dist_beta_v_quadratic;
  vector[num_treatments_cubic] dist_beta_v_cubic;
  
  matrix[use_dist_cluster_effects ? num_clusters : 0, num_treatments_linear] dist_beta_cluster_v_raw; // Linear distance*treatment effects
  matrix[use_dist_cluster_effects ? num_clusters : 0, num_treatments_quadratic] dist_beta_cluster_v_quadratic_raw;
  matrix[use_dist_cluster_effects ? num_clusters : 0, num_treatments_cubic] dist_beta_cluster_v_cubic_raw;
  
  row_vector<lower = 0>[use_dist_cluster_effects ? num_treatments_linear : 0] dist_beta_cluster_linear_sd;
  row_vector<lower = 0>[use_dist_cluster_effects ? num_treatments_quadratic : 0] dist_beta_cluster_quadratic_sd;
  row_vector<lower = 0>[use_dist_cluster_effects ? num_treatments_cubic : 0] dist_beta_cluster_cubic_sd;
 
  // Spline parameters 
  matrix<lower = u_splines_lb>[num_treatments_semiparam, num_knots_v] u_splines_v_raw;
  real<lower = 0> u_splines_v_sigma;
  
  matrix<lower = u_splines_lb>[num_treatments_semiparam, num_knots_v] u_splines_cluster_v_raw[use_dist_cluster_effects ? num_clusters : 0];
  vector<lower = 0>[use_dist_cluster_effects ? num_clusters : 0] u_splines_cluster_v_sigma;
  
  vector<lower = 0>[is_structural ? num_treatments : 0] mu_rep_raw;
  vector<lower = 0>[is_structural ? num_treatments : 0] dist_cost_k;
}

transformed parameters {
  vector[is_structural ? num_actual_treatments : 0] structural_beta;
  vector<lower = 0>[is_structural ? num_treatments : 0] mu_rep;
  
  matrix<lower = u_splines_lb>[num_treatments_semiparam, num_knots_v] u_splines_v;
  matrix<lower = u_splines_lb>[num_treatments_semiparam, num_knots_v] u_splines_cluster_v[use_dist_cluster_effects ? num_clusters : 0];
  
  vector[is_structural < 2 ? num_clusters : 0] cluster_obs_treatment_effect = cluster_treatment_design_matrix * beta; // v^* in the social signaling model
  vector[is_structural < 2 ? num_clusters : 0] cluster_obs_v = cluster_obs_treatment_effect; // v^* in the social signaling model
  vector<lower = 0, upper = 1>[is_structural < 2 ? num_clusters : 0] cluster_takeup_prob;
  
  vector[is_structural ? num_clusters : 0] structural_cluster_obs_v; 
  vector<lower = 0, upper = 1>[is_structural ? num_clusters : 0] structural_cluster_takeup_prob;
  
  if (use_cluster_effects && is_structural < 2) {
    matrix[num_clusters, num_actual_treatments] beta_cluster = beta_cluster_raw .* rep_matrix(beta_cluster_sd, num_clusters);
    
    cluster_obs_v += rows_dot_product(cluster_treatment_design_matrix, beta_cluster);
  }
 
  if (is_structural < 2) { 
    if (in_array(model_type, { MODEL_TYPE_LINEAR_DIST,
                               MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                               MODEL_TYPE_CUBIC_NONLINEAR_DIST })) {
      cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v) .* cluster_standard_dist; // Add a linear component to all models using a continuous distance effect
      
      if (use_dist_cluster_effects) {
        matrix[num_clusters, num_actual_treatments] dist_beta_cluster_linear = dist_beta_cluster_v_raw .* rep_matrix(dist_beta_cluster_linear_sd, num_clusters);
        
        cluster_obs_v += rows_dot_product(cluster_treatment_design_matrix, dist_beta_cluster_linear) .* cluster_standard_dist;
      }
      
      if (in_array(model_type, { MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                                 MODEL_TYPE_CUBIC_NONLINEAR_DIST })) {
        cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v_quadratic) .* cluster_standard_dist_square; 
        
        if (use_dist_cluster_effects) {
          matrix[num_clusters, num_actual_treatments] dist_beta_cluster_quadratic = dist_beta_cluster_v_quadratic_raw .* rep_matrix(dist_beta_cluster_quadratic_sd, num_clusters);
          
          cluster_obs_v += rows_dot_product(cluster_treatment_design_matrix, dist_beta_cluster_quadratic) .* cluster_standard_dist_square;
        }
        
        if (model_type == MODEL_TYPE_CUBIC_NONLINEAR_DIST) {
          cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v_cubic) .* cluster_standard_dist_cube;
          
          if (use_dist_cluster_effects) {
            matrix[num_clusters, num_actual_treatments] dist_beta_cluster_cubic = dist_beta_cluster_v_cubic_raw .* rep_matrix(dist_beta_cluster_cubic_sd, num_clusters);
            
            cluster_obs_v += rows_dot_product(cluster_treatment_design_matrix, dist_beta_cluster_cubic) .* cluster_standard_dist_cube;
          }
        }
      }
    } else if (model_type == MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_OSULLIVAN) { // O'Sullivan Splines 
      cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v) .* cluster_standard_dist; // Add a linear component to all models using a continuous distance effect
    
      for (treatment_index in 1:num_treatments) {
        u_splines_v[treatment_index] = u_splines_v_raw[treatment_index] * u_splines_v_sigma;
        
        if (use_dist_cluster_effects) {
          for (cluster_index in 1:num_clusters) {
            u_splines_cluster_v[cluster_index, treatment_index] = u_splines_cluster_v_raw[cluster_index, treatment_index] * u_splines_cluster_v_sigma[cluster_index];
          }
        }
      }
      
      if (use_dist_cluster_effects) {
        for (cluster_index in 1:num_clusters) {
          cluster_obs_v[cluster_index] += u_splines_cluster_v[cluster_index, cluster_assigned_treatment[cluster_index]] * Z_splines_v[cluster_index]';
        }
      } else {
        cluster_obs_v += rows_dot_product(u_splines_v[cluster_assigned_treatment], Z_splines_v);
      }
    } else if (model_type == MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_BSPLINE) {
      cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v) .* cluster_standard_dist; // Add a linear component to all models using a continuous distance effect
      
      for (treatment_index in 1:num_treatments) {
        u_splines_v[treatment_index] = cumulative_sum(append_col(1, rep_row_vector(u_splines_v_sigma, num_knots_v - 1)) .* u_splines_v_raw[treatment_index]); // Random walk
        
        if (use_dist_cluster_effects) {
          for (cluster_index in 1:num_clusters) {
            u_splines_cluster_v[cluster_index, treatment_index] = cumulative_sum(append_col(1, rep_row_vector(u_splines_cluster_v_sigma[cluster_index], num_knots_v - 1)) .* u_splines_cluster_v_raw[cluster_index, treatment_index]); // Random walk
          }
        }
      }
      
      if (use_dist_cluster_effects) {
        for (cluster_index in 1:num_clusters) {
          cluster_obs_v[cluster_index] += u_splines_cluster_v[cluster_index, cluster_assigned_treatment[cluster_index]] * Z_splines_v[cluster_index]';
        }
      } else {
        cluster_obs_v += rows_dot_product(u_splines_v[cluster_assigned_treatment], Z_splines_v);
      }
    } else if (model_type == MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_ISPLINE) {
      cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v) .* cluster_standard_dist; // Add a linear component to all models using a continuous distance effect
      
      for (treatment_index in 1:num_treatments) {
        u_splines_v[treatment_index, 1] = u_splines_v_raw[treatment_index, 1];
        
        for (knot_index in 2:num_knots_v) {
          u_splines_v[treatment_index, knot_index] = u_splines_v[treatment_index, knot_index - 1] * u_splines_v_raw[treatment_index, knot_index];
        }
      }
      
      cluster_obs_v += rows_dot_product(u_splines_v[cluster_assigned_treatment], Z_splines_v);
    }
  
    cluster_takeup_prob = Phi(- cluster_obs_v);
  }
  
  if (is_structural) {
    vector[num_treatments] structural_treatment_effect;
    vector[num_clusters] structural_cluster_treatment_effect;
    
    structural_beta = structural_beta_raw;
    structural_treatment_effect = treatment_map_design_matrix * structural_beta;
    structural_cluster_treatment_effect = structural_treatment_effect[cluster_assigned_treatment];
    
    if (use_cluster_effects) {
      matrix[num_clusters, num_treatments] structural_beta_cluster = structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_clusters);
      structural_cluster_treatment_effect += rows_dot_product(cluster_treatment_design_matrix, structural_beta_cluster);   
    } 
    
    mu_rep = mu_rep_raw;
    
    // Levels: control ink calendar bracelet
    // structural_beta[1] += structural_beta[2];
    // structural_beta[3] += structural_beta[1];
    // structural_beta[4] += structural_beta[1];
    
    // mu_rep[3] += mu_rep[1];
    // mu_rep[2] += mu_rep[3];
    // mu_rep[4] += mu_rep[3];
  
    if (is_structural == 1) {
      structural_cluster_obs_v = param_dist_cost(cluster_standard_dist, dist_cost_k[cluster_assigned_treatment]) 
        - structural_cluster_treatment_effect 
        - mu_rep[cluster_assigned_treatment] .* reputation_returns_normal_vec(cluster_obs_v);
    } else if (is_structural == 2) {
      int cluster_pos = 1;
      
      for (treatment_index in 1:num_treatments) {
        int num_treatment_clusters = treatment_cluster_size[treatment_index];
        int cluster_end = cluster_pos + num_treatment_clusters - 1;
        int treatment_cluster_ids[num_treatment_clusters] = treatment_sorted_clusters[cluster_pos:cluster_end];
        
        structural_cluster_obs_v[treatment_cluster_ids] = algebra_solver(cutoff_solution, 
                                                                         cluster_obs_v[treatment_cluster_ids],
                                                                         [ dist_cost_k[treatment_index], structural_treatment_effect[treatment_index], mu_rep[treatment_index] ]',
                                                                         to_array_1d(cluster_standard_dist[treatment_cluster_ids]),
                                                                         { num_treatment_clusters });
        
        // mu_rep[treatment_index] = algebra_solver(v_fixedpoint_solution, 
        //                                          [ 0 ]', 
        //                                          append_row([ dist_cost_k, structural_treatment_effect[treatment_index] ]', cluster_obs_v[treatment_cluster_ids]),
        //                                          to_array_1d(cluster_standard_dist[treatment_cluster_ids]),
        //                                          { num_treatment_clusters })[1];
        
        cluster_pos = cluster_end + 1;
      }
    }
    
    structural_cluster_takeup_prob = Phi(- structural_cluster_obs_v);
  }
}

model {
  u_splines_v_sigma ~ normal(0, u_splines_v_sigma_sd);
  
  if (use_cluster_effects && is_structural) {
    structural_beta_cluster_sd ~ normal(0, 0.25);
    to_vector(structural_beta_cluster_raw) ~ normal(0, 1);
  }
  
  if (is_structural < 2) {
    beta[1] ~ normal(0, 1);
    beta[2:num_actual_treatments] ~ normal(0, 1);
    
    if (use_cluster_effects) {
      beta_cluster_sd ~ normal(0, 0.25);
      to_vector(beta_cluster_raw) ~ normal(0, 1);
    }
  
    if (in_array(model_type, { MODEL_TYPE_LINEAR_DIST,
                               MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                               MODEL_TYPE_CUBIC_NONLINEAR_DIST })) {
      dist_beta_v ~ normal(0, dist_param_sd);
      
      if (use_dist_cluster_effects) {
        dist_beta_cluster_linear_sd ~ normal(0, 0.25);
        to_vector(dist_beta_cluster_v_raw) ~ normal(0, 1);
      }
      
      if (in_array(model_type, { MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                                 MODEL_TYPE_CUBIC_NONLINEAR_DIST })) {
        dist_beta_v_quadratic ~ normal(0, dist_param_sd);
        
        if (use_dist_cluster_effects) {
          dist_beta_cluster_quadratic_sd ~ normal(0, 0.25);
          to_vector(dist_beta_cluster_v_quadratic_raw) ~ normal(0, 1);
        }
        
        if (model_type == MODEL_TYPE_CUBIC_NONLINEAR_DIST) {
          dist_beta_v_cubic ~ normal(0, dist_param_sd);
          
          if (use_dist_cluster_effects) {
            dist_beta_cluster_cubic_sd ~ normal(0, 0.25);
            to_vector(dist_beta_cluster_v_cubic_raw) ~ normal(0, 1);
          }
        }
      }
    } else if (in_array(model_type, { MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_OSULLIVAN, 
                                      MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_BSPLINE })) { 
      dist_beta_v ~ normal(0, dist_param_sd);
      
      if (use_dist_cluster_effects) {
        dist_beta_cluster_linear_sd ~ normal(0, 0.25);
        to_vector(dist_beta_cluster_v_raw) ~ normal(0, 1);
      }
                                        
      to_vector(u_splines_v_raw) ~ std_normal();
      
      if (use_dist_cluster_effects) {
        u_splines_cluster_v_sigma ~ normal(0, u_splines_v_sigma_sd); 
        
        for (cluster_index in 1:num_clusters) {
          to_vector(u_splines_cluster_v_raw[cluster_index]) ~ std_normal();
        }
      }
    } else if (model_type == MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_ISPLINE) { 
      to_vector(u_splines_v_raw) ~ exponential(3);
    } 
    
    if (use_binomial) { 
      cluster_takeup_count[included_clusters] ~ binomial(cluster_size[included_clusters], cluster_takeup_prob[included_clusters]);
    } else {
      takeup[included_obs] ~ bernoulli(cluster_takeup_prob[obs_cluster_id[included_obs]]);
    }
  }
  
  if (is_structural) {
    structural_beta_raw ~ normal(0, 1);
    mu_rep_raw ~ std_normal();
    dist_cost_k ~ std_normal();
    
    if (use_binomial) { 
      cluster_takeup_count[included_clusters] ~ binomial(cluster_size[included_clusters], structural_cluster_takeup_prob[included_clusters]);
    } else {
      takeup[included_obs] ~ bernoulli(structural_cluster_takeup_prob[obs_cluster_id[included_obs]]);
    }
  }
}

generated quantities {
  matrix[num_grid_obs, num_actual_treatments] sim_v = rep_matrix((treatment_map_design_matrix * beta)', num_grid_obs); 
  matrix[num_grid_obs, num_actual_treatments] sim_takeup_prob;
  
  // matrix[num_grid_obs, num_treatments] sim_reputation;
  
  vector[use_binomial ? num_included_clusters : num_included_obs] log_lik; // loo
  vector[use_binomial ? num_excluded_clusters : num_excluded_obs] log_lik_heldout; // loo
  
  if (in_array(model_type, { MODEL_TYPE_LINEAR_DIST,
                             MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                             MODEL_TYPE_CUBIC_NONLINEAR_DIST })) {
    sim_v += rep_matrix((treatment_map_design_matrix * dist_beta_v)', num_grid_obs) .* rep_matrix(grid_dist, num_treatments); 
    
    if (in_array(model_type, { MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                               MODEL_TYPE_CUBIC_NONLINEAR_DIST })) {
      sim_v += rep_matrix((treatment_map_design_matrix * dist_beta_v_quadratic)', num_grid_obs) .* rep_matrix(grid_dist_square, num_treatments); 
      
      if (model_type == MODEL_TYPE_CUBIC_NONLINEAR_DIST) {
        sim_v += rep_matrix((treatment_map_design_matrix * dist_beta_v_cubic)', num_grid_obs) .* rep_matrix(grid_dist_cube, num_treatments); 
      }
    }
  } else if (in_array(model_type, { MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_OSULLIVAN, 
                                    MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_BSPLINE, 
                                    MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_ISPLINE })) {
                                      
    sim_v += rep_matrix((treatment_map_design_matrix * dist_beta_v)', num_grid_obs) .* rep_matrix(grid_dist, num_treatments); 
    
    for (treatment_index in 1:num_treatments) {
      sim_v[, treatment_index] += rows_dot_product(rep_matrix(u_splines_v[treatment_index], num_grid_obs), Z_grid_v);
       
      // for (grid_obs_index in 1:num_grid_obs) {
      //   sim_reputation[grid_obs_index, treatment_index] = reputational_returns_normal(sim_v[grid_obs_index, treatment_index]);
      // }
    }
  }

  sim_takeup_prob = Phi(- sim_v);
 
  if (use_binomial) {
    for (cluster_index_index in 1:num_included_clusters) {
      int cluster_index = included_clusters[cluster_index_index];
      
      log_lik[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], cluster_takeup_prob[cluster_index]);
    }
    
    for (cluster_index_index in 1:num_excluded_clusters) {
      int cluster_index = excluded_clusters[cluster_index_index];
      
      log_lik_heldout[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], cluster_takeup_prob[cluster_index]);
    }
  } else {
    for (obs_index_index in 1:num_included_obs) {
      int obs_index = included_obs[obs_index_index];
      
      log_lik[obs_index_index] = bernoulli_lpmf(takeup[obs_index] | cluster_takeup_prob[obs_cluster_id[obs_index]]);
    }
    
    for (obs_index_index in 1:num_excluded_obs) {
      int obs_index = excluded_obs[obs_index_index];
      
      log_lik_heldout[obs_index_index] = bernoulli_lpmf(takeup[obs_index] | cluster_takeup_prob[obs_cluster_id[obs_index]]);
    }
  }
}
