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
  
  
  real total_error_sd = phi[1];
  real u_sd = phi[2];
  
  return [ find_fixedpoint_solution(benefit_cost, mu_rep, total_error_sd, u_sd, use_u_in_delta, x_r[1], x_r[2], x_r[3]) ]'; 
}

vector map_find_fixedpoint_solution(vector benefit_cost, vector mu_rep, 
                                     real total_error_sd, real u_sd,
                                     data int use_u_in_delta, data real alg_sol_rel_tol, data real alg_sol_f_tol, data real alg_sol_max_steps) {
  int num_clusters = num_elements(benefit_cost);
  vector[2] phi = [total_error_sd, u_sd]';
  vector[2] thetas[num_clusters];
  real x_rs[num_clusters, 3] = rep_array({ alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps }, num_clusters);
  int x_is[num_clusters, 1] = rep_array({ use_u_in_delta }, num_clusters);
  
  for (cluster_index in 1:num_clusters) {
    thetas[cluster_index] = [ benefit_cost[cluster_index], mu_rep[cluster_index] ]';
  }
  
  return map_rect(find_fixedpoint_solution_rect, phi, thetas, rep_array({ alg_sol_rel_tol, alg_sol_f_tol, alg_sol_max_steps }, num_clusters) , x_is);
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