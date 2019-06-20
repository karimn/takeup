data {
  int<lower = 1> num_obs;
  int<lower = 3> num_knots_B;
  int<lower = 3> num_knots_R;
  int<lower = 1> num_grid_obs;
  int<lower = 1> num_treatments;
  
  int<lower = 0, upper = 1> takeup[num_obs];
  vector[num_obs] standard_dist;
  int<lower = 1, upper = num_treatments> assigned_treatment[num_obs];
 
  matrix[num_obs, num_knots_B] Z_splines_B; // Needs to be monotonically increasing
  matrix[num_obs, num_knots_R] Z_splines_R; 
  
  vector[num_grid_obs] grid_dist;
  matrix[num_grid_obs, num_knots_B] Z_grid_B;
  matrix[num_grid_obs, num_knots_R] Z_grid_R;
}

transformed data {
  matrix<lower = 0, upper = 1>[num_treatments, num_treatments] treatment_map_design_matrix = rep_matrix(0, num_treatments, num_treatments);
  matrix<lower = 0, upper = 1>[num_obs, num_treatments] treatment_design_matrix = rep_matrix(0, num_obs, num_treatments);
  
  treatment_map_design_matrix[, 1] = rep_vector(1, num_treatments);
  
  for(treatment_index in 1:num_treatments) {
    treatment_map_design_matrix[treatment_index, treatment_index] = 1;
  }
  
  treatment_design_matrix = treatment_map_design_matrix[assigned_treatment];
}

parameters {
  vector[num_treatments] beta;
  
  vector[num_treatments] dist_beta;
  
  // simplex[num_knots_B] u_splines_B_simplex[num_treatments];  
  // vector<lower = 0>[num_treatments] u_splines_B_scale;  
  // simplex[num_knots_B] u_splines_B_simplex;  
  // real<lower = 0> u_splines_B_scale;  
  
  matrix[num_treatments, num_knots_R] u_splines_R_raw;
  real<lower = 0> u_splines_R_sigma;
}

transformed parameters {
  // matrix[num_treatments, num_knots_B] u_splines_B;
  
  // vector[num_knots_B] u_splines_B;
  
  vector[num_obs] B;
  
  matrix[num_treatments, num_knots_R] u_splines_R;
  vector[num_obs] R;

  for (treatment_index in 1:num_treatments) {
    // u_splines_B[treatment_index] = u_splines_B_scale[treatment_index] * u_splines_B_simplex[treatment_index]';
    
    u_splines_R[treatment_index] = u_splines_R_raw[treatment_index] * u_splines_R_sigma;
    
    // u_splines_R[treatment_index, 1] = u_splines_R_raw[treatment_index, 1];
    // 
    // for (knots_R_index in 2:num_knots_R) {
    //   u_splines_R[treatment_index, knots_R_index] = u_splines_R[treatment_index, knots_R_index - 1] + u_splines_R_raw[treatment_index, knots_R_index] * u_splines_R_sigma;
    // }
  }
  
  // u_splines_B = u_splines_B_scale * u_splines_B_simplex;
  
  // B = treatment_design_matrix * beta - rows_dot_product(u_splines_B[assigned_treatment], Z_splines_B);
  // B = treatment_design_matrix * beta - Z_splines_B * u_splines_B; 
  B = treatment_design_matrix * beta + (treatment_design_matrix * dist_beta) .* standard_dist;
  
  R = rows_dot_product(u_splines_R[assigned_treatment], Z_splines_R);
}

model {
  beta[1] ~ normal(0, 2);
  beta[2:num_treatments] ~ normal(0, 1);
  
  dist_beta ~ normal(0, 0.125);
  
  // u_splines_B_scale ~ exponential(1);
  
  u_splines_R_sigma ~ normal(0, 0.125);
  
  to_vector(u_splines_R_raw) ~ normal(0, 1);
  
  // to_vector(u_splines_R) ~ normal(0, u_splines_R_sigma);
  
  takeup ~ bernoulli_logit(B + R);
}

generated quantities {
  matrix[num_grid_obs, num_treatments] sim_B = rep_matrix((treatment_map_design_matrix * beta)', num_grid_obs) + 
    rep_matrix((treatment_map_design_matrix * dist_beta)', num_grid_obs) .* rep_matrix(grid_dist, num_treatments); 
  matrix[num_grid_obs, num_treatments] sim_R; 
  matrix[num_grid_obs, num_treatments] sim_latent;
  matrix[num_grid_obs, num_treatments] sim_takeup;

  for (treatment_index in 1:num_treatments) {
    // sim_B[, treatment_index] -= rows_dot_product(rep_matrix(u_splines_B[treatment_index], num_grid_obs), Z_grid_B); 
    // sim_B[, treatment_index] -= Z_grid_B * u_splines_B;
    sim_R[, treatment_index] = rows_dot_product(rep_matrix(u_splines_R[treatment_index], num_grid_obs), Z_grid_R);
  }

  sim_latent = sim_B + sim_R;
  sim_takeup = inv_logit(sim_latent);
}
