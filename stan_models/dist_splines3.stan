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
  
  real param_dist_cost(real dist, real k) {
    return (k * dist^2) / 2;
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
  int use_binomial;
  int use_cluster_effects;
  
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
  // vector[use_cluster_effects ? num_clusters : 0] cluster_effects;
  // real<lower = 0> cluster_effects_sd;
  
  vector[num_actual_treatments] beta; // Threatment intercepts
  
  matrix[use_cluster_effects ? num_clusters : 0, num_actual_treatments] beta_cluster_raw;
  row_vector<lower = 0>[use_cluster_effects ? num_actual_treatments : 0] beta_cluster_sd;
  
  // vector[model_type >= MODEL_TYPE_LINEAR_DIST ? num_treatments : 0] dist_beta_v; // Linear distance*treatment effects
  vector[num_treatments_linear] dist_beta_v; // Linear distance*treatment effects
  vector[num_treatments_quadratic] dist_beta_v_quadratic;
  vector[num_treatments_cubic] dist_beta_v_cubic;
 
  // Spline parameters 
  matrix<lower = u_splines_lb>[num_treatments_semiparam, num_knots_v] u_splines_v_raw;
  real<lower = 0> u_splines_v_sigma;
  
  real<lower = 0> beta_link_kappa;
  
  // vector<lower =  0>[num_treatments] mu_rep;
  // real<lower = 0> dist_cost_k;
}

transformed parameters {
  matrix<lower = u_splines_lb>[num_treatments_semiparam, num_knots_v] u_splines_v;
  vector[num_clusters] cluster_obs_v = cluster_treatment_design_matrix * beta; // v^* in the social signaling model
  // vector[num_obs] obs_v = cluster_obs_v[obs_cluster_id]; // cluster_treatment_design_matrix[obs_cluster_id] * beta + cluster_effects[obs_cluster_id]; // v^* in the social signaling model
  // vector[num_obs] obs_takeup_prob;
  vector[num_clusters] cluster_takeup_prob;
  
  if (use_cluster_effects) {
    matrix[num_clusters, num_actual_treatments] beta_cluster = beta_cluster_raw .* rep_matrix(beta_cluster_sd, num_clusters);
    
    cluster_obs_v += rows_dot_product(cluster_treatment_design_matrix, beta_cluster);
     
     // cluster_obs_v += cluster_effects;
  }
  
  if (in_array(model_type, { MODEL_TYPE_LINEAR_DIST,
                             MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                             MODEL_TYPE_CUBIC_NONLINEAR_DIST })) {
    // obs_v += (cluster_treatment_design_matrix[obs_cluster_id] * dist_beta_v) .* standard_dist; // Add a linear component to all models using a continuous distance effect
    cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v) .* cluster_standard_dist; // Add a linear component to all models using a continuous distance effect
    
    if (in_array(model_type, { MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                               MODEL_TYPE_CUBIC_NONLINEAR_DIST })) {
      // obs_v += (cluster_treatment_design_matrix[obs_cluster_id] * dist_beta_v_quadratic) .* standard_dist_square; 
      cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v_quadratic) .* cluster_standard_dist_square; 
      
      if (model_type == MODEL_TYPE_CUBIC_NONLINEAR_DIST) {
        // obs_v += (cluster_treatment_design_matrix[obs_cluster_id] * dist_beta_v_cubic) .* standard_dist_cube;
        cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v_cubic) .* cluster_standard_dist_cube;
      }
    }
  } else if (model_type == MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_OSULLIVAN) { // O'Sullivan Splines 
    cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v) .* cluster_standard_dist; // Add a linear component to all models using a continuous distance effect
  
    for (treatment_index in 1:num_treatments) {
      u_splines_v[treatment_index] = u_splines_v_raw[treatment_index] * u_splines_v_sigma;
    }
    
    // obs_v += rows_dot_product(u_splines_v[assigned_treatment], Z_splines_v);
    cluster_obs_v += rows_dot_product(u_splines_v[cluster_assigned_treatment], Z_splines_v);
  } else if (model_type == MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_BSPLINE) {
    cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v) .* cluster_standard_dist; // Add a linear component to all models using a continuous distance effect
    
    for (treatment_index in 1:num_treatments) {
      u_splines_v[treatment_index] = cumulative_sum(append_col(1, rep_row_vector(u_splines_v_sigma, num_knots_v - 1)) .* u_splines_v_raw[treatment_index]); // Random walk
    }
    
    // obs_v += rows_dot_product(u_splines_v[assigned_treatment], Z_splines_v);
    cluster_obs_v += rows_dot_product(u_splines_v[cluster_assigned_treatment], Z_splines_v);
  } else if (model_type == MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_ISPLINE) {
    cluster_obs_v += (cluster_treatment_design_matrix * dist_beta_v) .* cluster_standard_dist; // Add a linear component to all models using a continuous distance effect
    
    for (treatment_index in 1:num_treatments) {
      u_splines_v[treatment_index, 1] = u_splines_v_raw[treatment_index, 1];
      
      for (knot_index in 2:num_knots_v) {
        u_splines_v[treatment_index, knot_index] = u_splines_v[treatment_index, knot_index - 1] * u_splines_v_raw[treatment_index, knot_index];
      }
    }
    
    // obs_v += rows_dot_product(u_splines_v[assigned_treatment], Z_splines_v);
    cluster_obs_v += rows_dot_product(u_splines_v[cluster_assigned_treatment], Z_splines_v);
  }
  
  // obs_takeup_prob = Phi(- obs_v);
  cluster_takeup_prob = Phi(- cluster_obs_v);
}

model {
  beta[1] ~ normal(0, 1);
  beta[2:num_actual_treatments] ~ normal(0, 1);
  
  u_splines_v_sigma ~ normal(0, u_splines_v_sigma_sd);
  // cluster_effects_sd ~ normal(0, 1);
  beta_link_kappa ~ normal(0, 1);
  
  
  if (use_cluster_effects) {
    beta_cluster_sd ~ normal(0, 1);
    to_vector(beta_cluster_raw) ~ normal(0, 1);
    
    // for (cluster_index in 1:num_clusters) {
    //   beta_cluster[cluster_index] ~ normal(0, beta_cluster_sd);
    // }
    
    // cluster_effects ~ normal(0, cluster_effects_sd);
  }
 
  if (in_array(model_type, { MODEL_TYPE_LINEAR_DIST,
                             MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                             MODEL_TYPE_CUBIC_NONLINEAR_DIST })) {
    dist_beta_v ~ normal(0, dist_param_sd);
    
    if (in_array(model_type, { MODEL_TYPE_QUADRATIC_NONLINEAR_DIST, 
                               MODEL_TYPE_CUBIC_NONLINEAR_DIST })) {
      dist_beta_v_quadratic ~ normal(0, dist_param_sd);
      
      if (model_type == MODEL_TYPE_CUBIC_NONLINEAR_DIST) {
        dist_beta_v_cubic ~ normal(0, dist_param_sd);
      }
    }
  } else if (in_array(model_type, { MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_OSULLIVAN, 
                                    MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_BSPLINE })) { 
    to_vector(u_splines_v_raw) ~ std_normal();
  } else if (model_type == MODEL_TYPE_SEMIPARAM_NONLINEAR_DIST_ISPLINE) { 
    to_vector(u_splines_v_raw) ~ exponential(3);
  }
  
  // dist_cost_k ~ exponential(1);
  // mu_rep ~ exponential(1); 
 
  if (use_binomial) { 
    cluster_takeup_count ~ binomial(cluster_size, cluster_takeup_prob);
  } else {
    takeup ~ bernoulli(cluster_takeup_prob[obs_cluster_id]);
  }
}

generated quantities {
  matrix[num_grid_obs, num_actual_treatments] sim_v = rep_matrix((treatment_map_design_matrix * beta)', num_grid_obs); 
  matrix[num_grid_obs, num_actual_treatments] sim_takeup_prob;
  
  // matrix[num_grid_obs, num_treatments] sim_reputation;
  
  vector[use_binomial ? num_clusters : num_obs] log_lik; // loo
  
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
    for (treatment_index in 1:num_treatments) {
      sim_v[, treatment_index] += rows_dot_product(rep_matrix(u_splines_v[treatment_index], num_grid_obs), Z_grid_v);
       
      // for (grid_obs_index in 1:num_grid_obs) {
      //   sim_reputation[grid_obs_index, treatment_index] = reputational_returns_normal(sim_v[grid_obs_index, treatment_index]);
      // }
    }
  }

  sim_takeup_prob = Phi(- sim_v);
 
  if (use_binomial) {
    for (cluster_index in 1:num_clusters) {
      log_lik[cluster_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], cluster_takeup_prob[cluster_index]);
    }
  } else {
    for (obs_index in 1:num_obs) {
      log_lik[obs_index] = bernoulli_lpmf(takeup[obs_index] | cluster_takeup_prob[obs_cluster_id[obs_index]]);
    }
  }
}
