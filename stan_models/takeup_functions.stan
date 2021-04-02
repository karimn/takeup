#include /../multilvlr/util.stan

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

vector param_kappa_dist_cost(vector dist, vector k) {
  if (num_elements(k) == 1) {
    return (k[1] * square(dist)) / 2;
  } else {
    return (k .* square(dist)) / 2;
  }
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

vector param_dist_cost_with_splines(vector dist, vector linear_dist_cost, vector quadratic_dist_cost, matrix u_splines, matrix Z_splines) {
  return rows_dot_product(u_splines, Z_splines) + param_dist_cost(dist, linear_dist_cost, quadratic_dist_cost);
}

vector calculate_mu_rep(int[] treatment_ids, vector dist,
                        real base_mu_rep, real mu_beliefs_effect,
                        matrix design_matrix,
                        matrix beta, matrix dist_beta) {
  vector[rows(beta)] beliefs_latent = calculate_beliefs_latent_predictor(design_matrix[treatment_ids], beta, dist_beta, dist);
  
  return base_mu_rep * exp(mu_beliefs_effect * (beliefs_latent - beta[, 1])); // Remove intercept 
}

matrix calculate_mu_rep_deriv(int treatment_id, vector dist,
                              real base_mu_rep, real mu_beliefs_effect,
                              matrix design_matrix,
                              matrix beta, matrix dist_beta) {
  matrix[rows(beta), 2] mu_rep;
  
  mu_rep[, 1] = calculate_mu_rep({ treatment_id }, dist, base_mu_rep, mu_beliefs_effect, design_matrix, beta, dist_beta);
  mu_rep[, 2] = mu_rep[, 1] .* (mu_beliefs_effect * (dist_beta * design_matrix[treatment_id]')); 
  
  return mu_rep;
}
  

real expected_delta_part(real v, real xc, real[] theta, data real[] x_r, data int[] x_i) {
  real w = theta[1];
  real u_sd = theta[2];
 
  real std_wmv = (w - v) / u_sd; 
  real v_lpdf = normal_lpdf(v | 0, 1);
  real wmv_lcdf = normal_lcdf(w - v | 0, u_sd);
  
  // real wmv_cdf;
  // if (std_wmv > 5) { 
  //   wmv_cdf = 1; 
  // } else if (std_wmv < -5) {
  //   wmv_cdf = 0; 
  // } else {
    // wmv_cdf = Phi_approx(std_wmv);
  // }
  
  // return v * exp(v_lpdf) * wmv_cdf; 
  return v * exp(v_lpdf + wmv_lcdf);
}

real expected_delta(real w, real total_error_sd, real u_sd, data real[] x_r, data int[] x_i) {
  real F_w = Phi_approx(w / total_error_sd); 
 
  real delta_part = integrate_1d(expected_delta_part, negative_infinity(), positive_infinity(), { w, u_sd }, x_r, x_i, 0.00001);
  
  return - delta_part / (F_w * (1 - F_w));
}

real expected_delta_deriv_part(real v, real xc, real[] theta, data real[] x_r, data int[] x_i) {
  real w = theta[1];
  real u_sd = theta[2];
 
  real v_lpdf = normal_lpdf(v | 0, 1);
  real wmv_lpdf = normal_lpdf(w - v | 0, u_sd);
  
  return v * exp(v_lpdf + wmv_lpdf);
}

vector expected_delta_deriv(real w, real total_error_sd, real u_sd, data real[] x_r, data int[] x_i) {
  real F_w = Phi_approx(w / total_error_sd); 
 
  real delta = expected_delta(w, total_error_sd, u_sd, x_r, x_i);
  real delta_deriv_part = integrate_1d(expected_delta_deriv_part, negative_infinity(), positive_infinity(), { w, u_sd }, x_r, x_i, 0.00001);
  
  return [delta, - (delta_deriv_part + delta * (1 - 2 * F_w)) / (F_w * (1 - F_w)) ]';
}

vector v_fixedpoint_solution_normal(vector model_param, vector theta, data real[] x_r, data int[] x_i) {
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

real find_fixedpoint_solution(real benefit_cost, real mu_rep, 
                              real total_error_sd, real u_sd,
                              data int use_u_in_delta, data real alg_sol_rel_tol, data real alg_sol_f_tol, data real alg_sol_max_steps) {
  vector[4] solver_theta = [ benefit_cost, mu_rep, total_error_sd, u_sd ]';
  int x_i[1] = { use_u_in_delta };
  
  return algebra_solver(v_fixedpoint_solution_normal, [ - benefit_cost ]', solver_theta, { 0.0 }, x_i, alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps)[1];
}

vector find_fixedpoint_solution_rect(vector phi, vector theta, data real[] x_r, data int[] x_i) {
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
  vector[4] thetas[num_clusters];
  // real x_rs[num_clusters, 3] = rep_array({ alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps }, num_clusters);
  int x_is[num_clusters, 1] = rep_array({ use_u_in_delta }, num_clusters);
  
  for (cluster_index in 1:num_clusters) {
    thetas[cluster_index] = [ benefit_cost[cluster_index], mu_rep[cluster_index], total_error_sd[cluster_index], u_sd[cluster_index] ]';
  }
  
  return map_rect(find_fixedpoint_solution_rect, phi, thetas, rep_array({ alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps }, num_clusters) , x_is);
}

vector calculate_roc_diff(int treatment_id_left, int treatment_id_right, 
                          real w, real total_error_sd, real u_sd, real dist_beta, 
                          real mu_rep_left, real mu_rep_right, real mu_rep_deriv_left, real mu_rep_deriv_right, 
                          real delta, real delta_deriv) {
  real f_w = exp(normal_lpdf(w | 0, total_error_sd));
  real left = (- f_w * (dist_beta - mu_rep_deriv_left * delta)) / (1 + mu_rep_left * delta_deriv);
  real right = (- f_w * (dist_beta - mu_rep_deriv_right * delta)) / (1 + mu_rep_right * delta_deriv);
  
  return [ left, right, f_w ]';
} 
  
vector calculate_roc_diff_rect(vector phi, vector theta, data real[] x_r, data int[] x_i) {
  real benefit_cost_left = theta[1];
  real benefit_cost_right = theta[2];
  real total_error_sd = theta[3];
  real u_sd = theta[4];
  real dist_beta = theta[5];
  real mu_rep_left = theta[6];
  real mu_rep_right = theta[7];
  real mu_rep_deriv_left = theta[8];
  real mu_rep_deriv_right = theta[9];
  
  real w_left = find_fixedpoint_solution(benefit_cost_left, mu_rep_left, total_error_sd, u_sd, x_i[1], x_r[1], x_r[2], x_r[3]);
  real w_right = find_fixedpoint_solution(benefit_cost_right, mu_rep_right, total_error_sd, u_sd, x_i[1], x_r[1], x_r[2], x_r[3]);
  vector[2] delta_left= expected_delta_deriv(w_left, total_error_sd, u_sd, x_r, x_i);
  vector[2] delta_right = expected_delta_deriv(w_right, total_error_sd, u_sd, x_r, x_i);
   
  return append_row(
    append_row([w_left, w_right]', append_row(delta_left, delta_right)),
    calculate_roc_diff(x_i[2], x_i[3], w_right, total_error_sd, u_sd, dist_beta, mu_rep_left, mu_rep_right, mu_rep_deriv_left, mu_rep_deriv_right, delta_right[1], delta_right[2])
  ); 
}

matrix map_calculate_roc_diff(
    data int treatment_id_left, data int treatment_id_right, 
    vector benefit_cost_left, vector benefit_cost_right,
    vector total_error_sd, vector u_sd, vector dist_beta, 
    vector mu_left, vector mu_right, vector mu_deriv_left, vector mu_deriv_right, 
    data int use_u_in_delta, data real alg_sol_rel_tol, data real alg_sol_f_tol, data real alg_sol_max_steps) {
  int num_clusters = num_elements(benefit_cost_left);
  vector[1] phi = total_error_sd[1:1];
  vector[9] thetas[num_clusters];
  
  // real x_rs[num_clusters, 1] = rep_array(x_r, num_clusters); // There seems to be a bug in STAN where this doesn't work. I have to call rep_array() below
  int x_is[num_clusters, 3] = rep_array({ use_u_in_delta, treatment_id_left, treatment_id_right }, num_clusters);
  
  for (cluster_index in 1:num_clusters) {
    // thetas[cluster_index] = [ w[cluster_index], total_error_sd[cluster_index], u_sd[cluster_index], dist_beta[cluster_index], 
    thetas[cluster_index] = [ benefit_cost_left[cluster_index], benefit_cost_right[cluster_index], 
                              total_error_sd[cluster_index], u_sd[cluster_index], dist_beta[cluster_index], 
                              mu_left[cluster_index], mu_right[cluster_index], mu_deriv_left[cluster_index], mu_deriv_right[cluster_index]  ]';
  }
  
  return to_matrix(map_rect(calculate_roc_diff_rect, phi, thetas, rep_array({ alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps }, num_clusters), x_is), num_clusters, 9, 0);
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