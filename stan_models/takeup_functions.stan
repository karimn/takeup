real reputational_returns_normal(real v) {
  real mix_Phi_v = Phi_approx(v);
  
  return exp(normal_lpdf(v | 0, 1)) / (mix_Phi_v * (1 - mix_Phi_v));
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

vector param_dist_cost(vector dist, vector linear_dist_cost, vector quadratic_dist_cost) {
  int num_cost = num_elements(dist);
  vector[num_cost] cost = rep_vector(0, num_cost); 
  
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

vector param_dist_cost(vector dist, vector linear_dist_cost, vector quadratic_dist_cost, matrix u_splines, matrix Z_splines) {
  return rows_dot_product(u_splines, Z_splines) + param_dist_cost(dist, linear_dist_cost, quadratic_dist_cost);
}

matrix param_dist_cost(real dist, matrix linear_dist_cost, matrix quadratic_dist_cost) {
  // int num_cost = num_elements(dist);
  // vector[num_cost] cost = rep_vector(0, num_cost); 
  // int num_clusters = rows(linear_dist_cost);
  // int num_treatments = cols(linear_dist_cost);
  // matrix[num_clusters, num_treatments] cost = rep_matrix(0, num_clusters, num_treatments); 
  
  return dist * linear_dist_cost + square(dist) * quadratic_dist_cost;
}


// vector calculate_mu_rep_log(array[] int treatment_ids, vector dist,
//                         matrix design_matrix,
//                         matrix beta, matrix dist_beta) {
// }

// matrix calculate_mu_rep_log_deriv(int treatment_id, vector dist,
//                               matrix design_matrix,
//                               matrix beta, matrix dist_beta) {
//   matrix[rows(beta), 2] mu_rep;
  
//   mu_rep[, 1] = calculate_mu_rep_log({ treatment_id }, dist, design_matrix, beta, dist_beta);

//   if (rows(design_matrix) == 1) {
//     mu_rep_lin_pred = (beta * design_matrix[1]') + ((dist_beta * design_matrix[1]') .* dist);
//   } else {
//     mu_rep_lin_pred = rows_dot_product(design_matrix, beta) + (rows_dot_product(design_matrix, dist_beta) .* dist);
//   }
//   mu_rep[, 2] = rows_dot_product(design_matrix, dist_beta) / mu_rep_lin_pred;
//   return mu_rep;
// }

vector calculate_mu_rep(array[] int treatment_ids, vector dist,
                        real base_mu_rep, real mu_beliefs_effect,
                        matrix design_matrix,
                        matrix beta, matrix dist_beta, int mu_rep_type) {
    vector[rows(beta)] beliefs_latent = calculate_beliefs_latent_predictor(design_matrix[treatment_ids], beta, dist_beta, dist);
    if (mu_rep_type == 1) { // log
      return log(beliefs_latent) ;
    } else if (mu_rep_type == 2) { // linear
      return beliefs_latent;
    } else if (mu_rep_type == 4) { // mu = x \lambda, x = \hat{p}, \lambda = base_mu_rep
      return base_mu_rep * inv_logit(beliefs_latent);
    } else { // exp
      return base_mu_rep * exp(mu_beliefs_effect * (beliefs_latent - beta[, 1])); // Remove intercept 
    }

}




matrix calculate_mu_rep_deriv(int treatment_id, vector dist,
                              real base_mu_rep, real mu_beliefs_effect,
                              matrix design_matrix,
                              matrix beta, matrix dist_beta, int mu_rep_type) {
  
  matrix[rows(beta), 2] mu_rep;
  mu_rep[, 1] = calculate_mu_rep({ treatment_id }, dist, base_mu_rep, mu_beliefs_effect, design_matrix, beta, dist_beta, mu_rep_type);

  if (mu_rep_type == 1) { // log
    mu_rep[, 2] = (dist_beta * design_matrix[treatment_id]')  ./ exp(mu_rep[, 1]);
    return mu_rep;
  } else if (mu_rep_type == 2) { // linear
    mu_rep[, 2] = dist_beta * design_matrix[treatment_id]';
  } else if (mu_rep_type == 4) { // mu = x \lambda, x = \hat{p}, \lambda = base_mu_rep
    vector[rows(beta)] beliefs_latent = calculate_beliefs_latent_predictor(design_matrix[{treatment_id}], beta, dist_beta, dist);
    mu_rep[, 2] = base_mu_rep .* (dist_beta * design_matrix[treatment_id]') .* exp(-beliefs_latent) ./ (1 + exp(-beliefs_latent))^2;
  } else { // exp
    mu_rep[, 2] = mu_rep[, 1] .* (mu_beliefs_effect * (dist_beta * design_matrix[treatment_id]')); 
  } 
  if (mu_rep_type == 3) {
    reject("mu_rep_type = 3 not yet implemented.");
  }
  return mu_rep;
}

real expected_delta_part(real v, real xc, array[] real theta, data array[] real x_r, data array[] int x_i) {
  real w = theta[1];
  real u_sd = theta[2];
 
  real std_wmv = (w - v) / u_sd; 
  real v_lpdf = normal_lpdf(v | 0, 1);
  real wmv_lcdf = normal_lcdf(w - v | 0, u_sd);
  
  return v * exp(v_lpdf + wmv_lcdf);
}

real expected_delta(real w, real total_error_sd, real u_sd, data array[] real x_r, data array[] int x_i) {
  real F_w = Phi_approx(w / total_error_sd); 
  real r;

  r = (-1/u_sd) * exp(-0.5 * (w^2)/(1 + u_sd^2)) * (1/sqrt(2*pi())) * sqrt((u_sd^2)/(1 + u_sd^2));
 
  // real delta_part = integrate_1d(expected_delta_part, negative_infinity(), positive_infinity(), { w, u_sd }, x_r, x_i, 0.00001);

  return - r / (F_w * (1 - F_w));
}

real expected_delta_deriv_part(real v, real xc, array[] real theta, data array[] real x_r, data array[] int x_i) {
  real w = theta[1];
  real u_sd = theta[2];
 
  real v_lpdf = normal_lpdf(v | 0, 1);
  real wmv_lpdf = normal_lpdf(w - v | 0, u_sd);
  
  return v * exp(v_lpdf + wmv_lpdf);
}

// Calculate Delta'(w) (also returns Delta(w) as first argument)
vector expected_delta_deriv(real w, real total_error_sd, real u_sd, data array[] real x_r, data array[] int x_i) {
  real F_w = Phi_approx(w / total_error_sd); 
 
  real delta = expected_delta(w, total_error_sd, u_sd, x_r, x_i);
  // real delta_deriv_part = integrate_1d(expected_delta_deriv_part, negative_infinity(), positive_infinity(), { w, u_sd }, x_r, x_i, 0.00001);
  real Sigma = sqrt((u_sd^2)/(1 + u_sd^2));
  real mu = w/(u_sd^2 + 1);
  real H = (1/u_sd) * (1/sqrt(2*pi())) * exp(-0.5 *((w^2)/(u_sd^2 + 1))) * Sigma;
   
  return [delta, - (H*mu + exp(normal_lpdf(w | 0, total_error_sd)) * delta * (1 - 2 * F_w)) / (F_w * (1 - F_w)) ]';
}

vector v_fixedpoint_solution_normal(vector model_param, vector theta, data array[] real x_r, data array[] int x_i) {
  real cutoff = model_param[1];
  
  int use_u_in_delta = x_i[1];
  
  real benefit_cost = theta[1];
  real mu = theta[2];
  real total_error_sd = theta[3];
  real u_sd = theta[4];
  
  real delta;
  
  if (use_u_in_delta && u_sd > 0) {
    delta = expected_delta(cutoff, total_error_sd, u_sd, x_r, x_i);
  } else {
    delta = reputational_returns_normal(cutoff);
  }
  
  return [ cutoff + benefit_cost + mu * delta ]';
}

real mixed_binomial_lpmf(array[] int outcomes, vector lambda, array[] int N, matrix prob) {
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

real mixed_bernoulli_lpmf(array[] int outcomes, vector lambda, matrix prob) {
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

real find_fixedpoint_solution(real benefit_cost, real mu_rep, 
                              real total_error_sd, real u_sd,
                              data int use_u_in_delta, data real alg_sol_rel_tol, data real alg_sol_f_tol, data real alg_sol_max_steps) {
  vector[4] solver_theta = [ benefit_cost, mu_rep, total_error_sd, u_sd ]';
  array[1] int x_i = { use_u_in_delta };
  
  return algebra_solver(v_fixedpoint_solution_normal, [ - benefit_cost ]', solver_theta, { 0.0 }, x_i, alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps)[1];
}

vector find_fixedpoint_solution_rect(vector phi, vector theta, data array[] real x_r, data array[] int x_i) {
  int use_u_in_delta = x_i[1];
  
  real benefit_cost = theta[1];
  real mu_rep = theta[2];
  real total_error_sd = theta[3];
  real u_sd = theta[4];
  
  return [ find_fixedpoint_solution(benefit_cost, mu_rep, total_error_sd, u_sd, use_u_in_delta, x_r[1], x_r[2], x_r[3]) ]'; 
}

vector map_find_fixedpoint_solution(vector benefit_cost, vector mu_rep, 
                                     vector total_error_sd, vector u_sd,
                                     data int use_u_in_delta, data real alg_sol_rel_tol, data real alg_sol_f_tol, data real alg_sol_max_steps) {
  int num_clusters = num_elements(benefit_cost);
  vector[1] phi = total_error_sd[1:1];
  array[num_clusters] vector[4] thetas;
  array[num_clusters, 1] int x_is = rep_array({ use_u_in_delta }, num_clusters);
  
  for (cluster_index in 1:num_clusters) {
    thetas[cluster_index] = [ benefit_cost[cluster_index], mu_rep[cluster_index], total_error_sd[cluster_index], u_sd[cluster_index] ]';
  }
  
  return map_rect(find_fixedpoint_solution_rect, phi, thetas, rep_array({ alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps }, num_clusters) , x_is);
}

real calculate_roc(real w, real total_error_sd, real dist_beta, 
                     real mu_rep, real mu_rep_deriv, 
                     real delta, real delta_deriv) {
  real f_w = exp(normal_lpdf(w | 0, total_error_sd));
  return (- f_w * (dist_beta - mu_rep_deriv * delta)) / (1 + mu_rep * delta_deriv);
} 

vector calculate_roc_rect(vector phi, vector theta, data array[] real x_r, data array[] int x_i) {
  real benefit_cost = theta[1];
  real benefit_cost_control = theta[2];
  real total_error_sd = theta[3];
  real total_error_sd_control = theta[4];
  real u_sd = theta[5];
  real u_sd_control = theta[6];
  real dist_beta = theta[7];
  real mu_rep = theta[8];
  real mu_rep_control = theta[9];
  real mu_rep_deriv = theta[10];
  
  real w = find_fixedpoint_solution(benefit_cost, mu_rep, total_error_sd, u_sd, x_i[1], x_r[1], x_r[2], x_r[3]);
  real w_control = find_fixedpoint_solution(benefit_cost_control, mu_rep_control, total_error_sd_control, u_sd_control, x_i[1], x_r[1], x_r[2], x_r[3]);
   
  vector[2] delta = expected_delta_deriv(w_control, total_error_sd, u_sd, x_r, x_i);

  real roc_no_vis = -exp(normal_lpdf(w_control | 0, total_error_sd))*dist_beta;

  return append_row(
  append_row(
    append_row([w, w_control]', delta),
    calculate_roc(w_control, total_error_sd, dist_beta, mu_rep, mu_rep_deriv, delta[1], delta[2])),
    [roc_no_vis]'
  ); 
}

matrix map_calculate_roc(
    vector benefit_cost, 
    vector benefit_cost_control, 
    vector total_error_sd, vector total_error_sd_control, vector u_sd, vector u_sd_control, vector dist_beta, 
    vector mu, vector mu_control, vector mu_deriv,
    data int use_u_in_delta, data real alg_sol_rel_tol, data real alg_sol_f_tol, data real alg_sol_max_steps) {
  int num_clusters = num_elements(benefit_cost_control);
  vector[1] phi = total_error_sd[1:1];
  array[num_clusters] vector[10] thetas;
  
  array[num_clusters, 1] int x_is = rep_array({ use_u_in_delta }, num_clusters);
  
  for (cluster_index in 1:num_clusters) {
    thetas[cluster_index] = [ benefit_cost[cluster_index], benefit_cost_control[cluster_index],
                              total_error_sd[cluster_index], total_error_sd_control[cluster_index], u_sd[cluster_index], u_sd_control[cluster_index], dist_beta[cluster_index], 
                              mu[cluster_index], mu_control[cluster_index], mu_deriv[cluster_index] ]';
  }
  
  return to_matrix(map_rect(calculate_roc_rect, phi, thetas, rep_array({ alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps }, num_clusters), x_is), num_clusters, 6, 0);
}


array[] int prepare_cluster_assigned_dist_group_treatment(array[] int cluster_assigned_treatment, array[] int cluster_assigned_dist_group) {
  int num_clusters = num_elements(cluster_assigned_treatment);
  int num_treatments = max(cluster_assigned_treatment);
  int num_dist = max(cluster_assigned_dist_group);
  
  array[num_clusters] int treatment_id;
 
  array[num_treatments, num_dist] int treatment_id_map;
  
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
