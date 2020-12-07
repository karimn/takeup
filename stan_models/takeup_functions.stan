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

real expected_delta_part(real v, real xc, real[] theta, data real[] x_r, data int[] x_i) {
  real w = theta[1];
  real u_sd = theta[2];

  return v * exp(normal_lcdf(w - v | 0, u_sd) + normal_lpdf(v | 0, 1));
}

real expected_delta(real w, real u_sd) {
  real delta_part = integrate_1d(expected_delta_part, negative_infinity(), positive_infinity(), { w, u_sd }, { 0.0 }, { 0 }, 0.01);
  real F_w = Phi_approx(w / u_sd);

  return - delta_part / (F_w * (1 - F_w));
}

vector v_fixedpoint_solution_normal(vector model_param, vector theta, data real[] x_r, data int[] x_i) {
  real cutoff = model_param[1];
  
  real benefit_cost = theta[1];
  real mu = theta[2];
  
  int num_v_mix = x_i[1];
  int use_u_in_delta = x_i[2];
  
  vector[num_v_mix] lambda = theta[3:(3 + num_v_mix - 1)];
  vector[num_v_mix] mix_mean = theta[(3 + num_v_mix):(3 + 2 * num_v_mix - 1)];
  vector[num_v_mix] mix_sd = theta[(3 + 2 * num_v_mix):(3 + 3 * num_v_mix - 1)];
  real u_sd = theta[3 + 3 * num_v_mix]; 
  
  if (use_u_in_delta && u_sd > 0) {
    // return [ cutoff + benefit_cost + mu * expected_delta(cutoff, u_sd) ]';
    reject("Not supported yet."); // BUG This is temporary because we still can't compile the above line.
    return [ cutoff + benefit_cost + mu * reputational_returns_normal(cutoff, lambda, mix_mean, mix_sd) ]';
  } else {
    return [ cutoff + benefit_cost + mu * reputational_returns_normal(cutoff, lambda, mix_mean, mix_sd) ]';
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