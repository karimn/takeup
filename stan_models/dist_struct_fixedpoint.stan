#include /takeup_header.stan

parameters {
  // Levels: control ink calendar bracelet
  vector[use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_discrete_dist : 1] beta_control;
  vector[use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_discrete_dist : 1] beta_ink_effect;
  vector<lower = (use_private_incentive_restrictions ? 0 : negative_infinity())>[use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_discrete_dist : 1] beta_calendar_effect;
  vector<lower = (use_private_incentive_restrictions ? 0 : negative_infinity())>[use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_discrete_dist : 1] beta_bracelet_effect;
  
  matrix[use_cluster_effects ? num_clusters : 0, use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_treatments * num_discrete_dist : num_treatments] structural_beta_cluster_raw;
  row_vector<lower = 0>[use_cluster_effects ? (use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_treatments * num_discrete_dist : num_treatments) : 0] structural_beta_cluster_sd;
  
  matrix[use_county_effects ? num_counties : 0, (use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_dist_group_treatments : num_treatments)] structural_beta_county_raw;
  row_vector<lower = 0>[use_county_effects ? (use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_dist_group_treatments : num_treatments) : 0] structural_beta_county_sd;
  
  // Salience
  
  real<lower = 0> beta_salience;
  real<lower = 0> dist_beta_salience;
  real<lower = 0> dist_quadratic_beta_salience;
  
  // Name Matched
  
  // vector[num_treatments_name_matched] beta_nm_effect;
  // 
  // matrix[use_cluster_effects && use_name_matched_obs ? num_clusters : 0, (use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_treatments * num_discrete_dist : num_treatments)] beta_nm_effect_cluster_raw;
  // row_vector<lower = 0>[use_cluster_effects && use_name_matched_obs ? (use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_treatments * num_discrete_dist : num_treatments) : 0] beta_nm_effect_cluster_sd;
  // 
  // matrix[use_county_effects && use_name_matched_obs ? num_clusters : 0, (use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_treatments * num_discrete_dist : num_treatments)] beta_nm_effect_county_raw;
  // row_vector<lower = 0>[use_county_effects && use_name_matched_obs ? (use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_treatments * num_discrete_dist : num_treatments) : 0] beta_nm_effect_county_sd;
  
  // V Mixture
  
  vector<lower = 0.5, upper = 1>[num_v_mix - 1] recursive_lambda_v_mix;
  
  vector<lower = 0>[num_v_mix - 1] v_mix_mean_diff;
  vector<lower = 0>[num_v_mix - 1] v_mix_sd;
  
  // Reputational Returns
  
  real v_mu;
  
  row_vector<lower = 0>[suppress_reputation ? 0 : num_treatments] mu_rep_raw;
  // real<lower = 0, upper = 1> v_sd;
  vector<lower = 0>[suppress_shocks ? 0 : num_shocks - 1] ub_ur_sd;
  cholesky_factor_corr[suppress_shocks ? 0 : num_shocks] L_all_u_corr;
  
  matrix[use_mu_cluster_effects && !suppress_reputation ? num_clusters : 0, num_treatments] mu_cluster_effects_raw;
  row_vector<lower = 0>[use_mu_cluster_effects && !suppress_reputation ? num_treatments : 0] mu_cluster_effects_sd;
  
  matrix[use_mu_county_effects && !suppress_reputation ? num_counties : 0, num_treatments] mu_county_effects_raw;
  row_vector<lower = 0>[use_mu_county_effects && !suppress_reputation ? num_treatments : 0] mu_county_effects_sd;
  
  // Parameteric Cost Model
  
  vector<lower = 0>[num_treatments_param_kappa] dist_cost_k_raw;
  
  // Linear Parametric Cost
  
  vector<lower = (use_dist_salience ? 0 : negative_infinity())>[num_treatments_param] dist_beta_v; // Linear distance*treatment effects
  
  matrix[use_param_dist_cluster_effects ? num_clusters : 0, num_treatments_param] dist_beta_cluster_raw;
  row_vector<lower = 0>[use_param_dist_cluster_effects ? num_treatments_param : 0] dist_beta_cluster_sd;
  
  matrix[use_param_dist_county_effects ? num_counties : 0, num_treatments_param] dist_beta_county_raw;
  row_vector<lower = 0>[use_param_dist_county_effects ? num_treatments_param : 0] dist_beta_county_sd;
  
  // Quadratic Cost Model
  
  vector<lower = 0>[num_treatments_param_quadratic] dist_quadratic_beta_v; // Quadratic distance*treatment effects
  
  // Semiparameteric Cost Model
  
  matrix[num_treatments_semiparam, num_knots_v] u_splines_v_raw;
  real<lower = 0> u_splines_v_sigma;
  
  // Assigned Distance Model
  
  simplex[num_dist_group_mix] group_dist_mix[num_discrete_dist];
  
  // vector<lower = 0>[num_discrete_dist] group_dist_mean;
  // vector<lower = 0>[num_discrete_dist] group_dist_sd;
  ordered[num_dist_group_mix] group_dist_mean[num_discrete_dist];
  vector<lower = 0>[num_dist_group_mix] group_dist_sd[num_discrete_dist];
  
  matrix<lower = 0>[num_clusters, num_discrete_dist - 1] missing_cluster_standard_dist; 
}

transformed parameters {
  vector[num_dist_group_treatments] beta; 
  vector[num_dist_group_treatments] structural_treatment_effect;
  vector[num_clusters] structural_cluster_benefit_cost;
  matrix[num_clusters, num_dist_group_treatments] structural_beta_cluster = rep_matrix(0, num_clusters, num_dist_group_treatments);
  // matrix[num_clusters, num_actual_treatments] structural_beta_cluster = rep_matrix(0, num_clusters, num_actual_treatments);
  // matrix[num_clusters, num_actual_treatments] structural_beta_cluster_nm = rep_matrix(0, num_clusters, num_actual_treatments);
  matrix[num_counties, num_dist_group_treatments] structural_beta_county = rep_matrix(0, num_counties, num_dist_group_treatments);
  // matrix[num_counties, num_actual_treatments] structural_beta_county = rep_matrix(0, num_counties, num_actual_treatments);
  // matrix[num_counties, num_actual_treatments] structural_beta_county_nm = rep_matrix(0, num_counties, num_actual_treatments);
  
  row_vector<lower = 0>[suppress_reputation ? 0 : num_treatments] mu_rep = mu_rep_raw;
  matrix<lower = 0>[!suppress_reputation ? num_clusters : 0, num_treatments] cluster_mu_rep;
  
  vector[num_clusters] structural_cluster_obs_v = rep_vector(0, num_clusters);
  vector[use_name_matched_obs ? num_clusters : 0] structural_cluster_obs_v_nm;
  matrix<lower = 0, upper = 1>[num_v_mix, num_clusters] structural_cluster_takeup_prob;
  matrix<lower = 0, upper = 1>[num_v_mix, use_name_matched_obs ? num_clusters : 0] structural_cluster_takeup_prob_nm;
  
  simplex[num_v_mix] lambda_v_mix;
  ordered[num_v_mix - 1] v_mix_mean = cumulative_sum(v_mix_mean_diff);
  
  vector<lower = 0>[num_treatments_param_kappa] dist_cost_k;
  // vector[num_treatments] linear_dist_cost = rep_vector(0, num_treatments);
  // vector[num_treatments] quadratic_dist_cost = rep_vector(0, num_treatments);
  vector[num_dist_group_treatments] linear_dist_cost = rep_vector(0, num_dist_group_treatments);
  vector[num_dist_group_treatments] quadratic_dist_cost = rep_vector(0, num_dist_group_treatments);
  // matrix[num_clusters, num_treatments] cluster_linear_dist_cost = rep_matrix(0, num_clusters, num_treatments);
  matrix[num_clusters, num_dist_group_treatments] cluster_linear_dist_cost = rep_matrix(0, num_clusters, num_dist_group_treatments);
  matrix[num_treatments, num_knots_v] u_splines_v = rep_matrix(0, num_treatments, num_knots_v);
  vector[num_clusters] cluster_dist_cost = rep_vector(0, num_clusters);
  
  matrix<lower = 0>[num_clusters, num_discrete_dist] cf_cluster_standard_dist; 
  
  matrix[num_shocks, num_shocks] L_all_u_vcov;
  real total_error_sd;
  
  for (cluster_index in 1:num_clusters) {
    int dist_group_pos = 1;
    
    for (dist_group_index in 1:num_discrete_dist) {
      if (cluster_treatment_map[cluster_assigned_dist_group_treatment[cluster_index], 2] == dist_group_index) {
        cf_cluster_standard_dist[cluster_index, dist_group_index] = cluster_standard_dist[cluster_index];
      } else {
        cf_cluster_standard_dist[cluster_index, dist_group_index] = missing_cluster_standard_dist[cluster_index, dist_group_pos];
        
        dist_group_pos += 1;
      }
    }
  }
  
  for (dist_index in 1:num_discrete_dist) {
    if (use_cost_model == COST_MODEL_TYPE_DISCRETE) {
      beta[(num_treatments * (dist_index - 1) + 1):(num_treatments * dist_index)] = 
        [ beta_control[dist_index], beta_ink_effect[dist_index], beta_calendar_effect[dist_index], beta_bracelet_effect[dist_index] ]';
    } else {
      beta[(num_treatments * (dist_index - 1) + 1):(num_treatments * dist_index)] = 
        [ beta_control[1] * (dist_index == 1), beta_ink_effect[1], beta_calendar_effect[1], beta_bracelet_effect[1] ]';
    }
  }
 
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
    
    cluster_mu_rep = rep_matrix(mu_rep, num_clusters);
    
    if (use_mu_cluster_effects) {
      matrix[num_clusters, num_treatments] mu_cluster_effects =  mu_cluster_effects_raw .* rep_matrix(mu_cluster_effects_sd, num_clusters);
      
      cluster_mu_rep = cluster_mu_rep .* exp(mu_cluster_effects);
    }
    
    if (use_mu_county_effects) {
      matrix[num_counties, num_treatments] mu_county_effects =  mu_county_effects_raw .* rep_matrix(mu_county_effects_sd, num_counties);
      
      cluster_mu_rep = cluster_mu_rep .* exp(mu_county_effects[cluster_county_id]);
    } 
  }
  
  if (use_salience_effect) {
    // structural_treatment_effect += beta_salience * mu_rep[cluster_treatment_map[, 1]]';
  }
  
  if (use_cost_model == COST_MODEL_TYPE_PARAM_KAPPA) { 
    dist_cost_k = dist_cost_k_raw;
   
    if (use_cost_k_restrictions) { 
      dist_cost_k[2] += dist_cost_k[1];
      dist_cost_k[3] += dist_cost_k[1];
      dist_cost_k[4] += dist_cost_k[1];
    }
    
    cluster_dist_cost = param_kappa_dist_cost(cluster_standard_dist, dist_cost_k[cluster_treatment_map[, 1]]);
  } else if (use_cost_model != COST_MODEL_TYPE_DISCRETE) {
    if (in_array(use_cost_model, { COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE, COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE, COST_MODEL_TYPE_SEMIPARAM_SALIENCE })) {
      linear_dist_cost = rep_vector(dist_beta_v[1], num_dist_group_treatments) + dist_beta_salience * append_col(mu_rep, mu_rep)';
      
      if (use_param_dist_cluster_effects) {
        cluster_linear_dist_cost += rep_matrix(dist_beta_cluster_raw[, 1] * dist_beta_cluster_sd[1], num_dist_group_treatments);
      }
      
      if (use_param_dist_county_effects) {
        cluster_linear_dist_cost += rep_matrix((dist_beta_county_raw[, 1] * dist_beta_county_sd[1])[cluster_county_id], num_dist_group_treatments);
      }
      
      if (use_cost_model == COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE) {
        quadratic_dist_cost = rep_vector(dist_quadratic_beta_v[1], num_dist_group_treatments) + dist_quadratic_beta_salience * append_col(mu_rep, mu_rep)';
      }
    } else {
      if (use_single_cost_model) {
        linear_dist_cost = rep_vector(dist_beta_v[1], num_dist_group_treatments);
        
        if (use_param_dist_cluster_effects) {
          cluster_linear_dist_cost += rep_matrix(dist_beta_cluster_raw[, 1] * dist_beta_cluster_sd[1], num_dist_group_treatments);
        }
        
        if (use_param_dist_county_effects) {
          cluster_linear_dist_cost += rep_matrix((dist_beta_county_raw[, 1] * dist_beta_county_sd[1])[cluster_county_id], num_dist_group_treatments);
        }
        
        if (num_treatments_param_quadratic > 0) {
          quadratic_dist_cost = rep_vector(dist_quadratic_beta_v[1], num_dist_group_treatments);
        }
      } else {
        linear_dist_cost = treatment_map_design_matrix * append_row(append_row(dist_beta_v, 0), dist_beta_v[2:]);
        
        if (use_param_dist_cluster_effects) {
          cluster_linear_dist_cost += dist_beta_cluster_raw .* rep_matrix(dist_beta_cluster_sd, num_clusters);
        }
        
        if (use_param_dist_county_effects) {
          cluster_linear_dist_cost += (dist_beta_county_raw .* rep_matrix(dist_beta_county_sd, num_counties))[cluster_county_id];
        }
        
        if (num_treatments_param_quadratic > 0) {
          quadratic_dist_cost = treatment_map_design_matrix * append_row(append_row(dist_quadratic_beta_v, 0), dist_quadratic_beta_v[2:]);
        }
      }
      
      cluster_linear_dist_cost += rep_matrix(linear_dist_cost', num_clusters);
    } 
    
    if (use_semiparam) {
      for (treatment_index in 1:num_treatments_semiparam) {
        if (use_dist_salience) {
          u_splines_v[treatment_index] = u_splines_v_raw[1] * u_splines_v_sigma * mu_rep[treatment_index];
        } else {
          u_splines_v[treatment_index] = u_splines_v_raw[treatment_index] * u_splines_v_sigma;
        }
      }        
    } 
      
    cluster_dist_cost = param_dist_cost(cluster_standard_dist, 
                                        to_vector(cluster_linear_dist_cost)[long_cluster_by_treatment_index],
                                        quadratic_dist_cost[cluster_treatment_map[cluster_assigned_dist_group_treatment, 1]],
                                        u_splines_v[cluster_treatment_map[cluster_assigned_dist_group_treatment, 1]],
                                        Z_splines_v[cluster_treatment_map[cluster_assigned_dist_group_treatment, 1]]);
  }
  
  structural_cluster_benefit_cost = structural_treatment_effect[cluster_assigned_dist_group_treatment] - cluster_dist_cost;
  
  if (use_cluster_effects) {
    vector[num_clusters] cluster_effects;
    
    if (use_cost_model != COST_MODEL_TYPE_DISCRETE) {
      structural_beta_cluster[, 1:num_treatments] = structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_clusters);
      structural_beta_cluster[, (num_treatments + 2):num_dist_group_treatments] = structural_beta_cluster[, 2:num_treatments];
    } else {
      structural_beta_cluster = structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_clusters);
    }
    
    cluster_effects = rows_dot_product(cluster_treatment_design_matrix, structural_beta_cluster); 
    structural_cluster_benefit_cost += cluster_effects;
  }
  
  if (use_county_effects) {
    vector[num_clusters] county_effects;
    
    if (use_cost_model != COST_MODEL_TYPE_DISCRETE) {
      structural_beta_county[, 1:num_treatments] = structural_beta_county_raw .* rep_matrix(structural_beta_county_sd, num_counties);
      structural_beta_county[, (num_treatments + 2):num_dist_group_treatments] = structural_beta_county[, 2:num_treatments];
    } else {
      structural_beta_county = structural_beta_county_raw .* rep_matrix(structural_beta_county_sd, num_counties);
    }
    
    county_effects = rows_dot_product(cluster_treatment_design_matrix, structural_beta_county[cluster_county_id]); 
    structural_cluster_benefit_cost += county_effects;
  }

  if (suppress_reputation) {
    structural_cluster_obs_v = - structural_cluster_benefit_cost;
  } else {
    for (cluster_index in 1:num_clusters) {
      structural_cluster_obs_v[cluster_index] = algebra_solver(v_fixedpoint_solution_normal,
                                                               [ - structural_cluster_benefit_cost[cluster_index] ]',
                                                               prepare_solver_theta(structural_cluster_benefit_cost[cluster_index], 
                                                                                    cluster_mu_rep[cluster_index, cluster_treatment_map[cluster_assigned_dist_group_treatment[cluster_index], 1]],
                                                                                    use_shifting_v_dist ? v_mu : 0,
                                                                                    lambda_v_mix,
                                                                                    v_mix_mean,
                                                                                    1,
                                                                                    v_mix_sd),
                                                               { 0.0 },
                                                               { num_v_mix }, 
                                                               1e-10,
                                                               1e-5,
                                                               1e6)[1];
                                                               
    }
  }
  
  // if (use_name_matched_obs) {
  //   vector[num_clusters] cluster_effects_nm = rep_vector(0, num_clusters);
  //   vector[num_clusters] county_effects_nm = rep_vector(0, num_clusters);
  //   
  //   if (use_cluster_effects) {
  //     structural_beta_cluster_nm = beta_nm_effect_cluster_raw .* rep_matrix(beta_nm_effect_cluster_sd, num_clusters);
  //     cluster_effects_nm = rows_dot_product(cluster_treatment_design_matrix, structural_beta_cluster_nm); 
  //   }
  //   
  //   if (use_county_effects) {
  //     structural_beta_county_nm = beta_nm_effect_county_raw .* rep_matrix(beta_nm_effect_county_sd, num_counties);
  //     county_effects_nm = rows_dot_product(cluster_treatment_design_matrix, structural_beta_county_nm[cluster_county_id]); 
  //   }
  //   
  //   structural_cluster_obs_v_nm = structural_cluster_obs_v 
  //                               - (structural_treatment_effect_nm[cluster_assigned_treatment] + cluster_effects_nm + county_effects_nm); 
  // }
  
  if (!suppress_shocks) {
    vector[num_shocks] all_u_sd = append_row(1, ub_ur_sd);
    
    L_all_u_vcov = diag_pre_multiply(all_u_sd, L_all_u_corr);
    total_error_sd = sqrt(sum(L_all_u_vcov * L_all_u_vcov'));
  } else {
    total_error_sd = 1;
  }
  
  for (mix_index in 1:num_v_mix) {
    if (mix_index == 1) {
      structural_cluster_takeup_prob[1] = Phi(- structural_cluster_obs_v / total_error_sd)';
      structural_cluster_takeup_prob_nm[1] = Phi(- structural_cluster_obs_v_nm / total_error_sd)';
      
      // structural_obs_takeup_prob[1] = structural_cluster_takeup_prob[1, obs_cluster_id]; // BUGBUG no name matching yet 
    } else {
      structural_cluster_takeup_prob[mix_index] = Phi(- (structural_cluster_obs_v - v_mix_mean[mix_index - 1]) / v_mix_sd[mix_index - 1])';
      structural_cluster_takeup_prob_nm[mix_index] = Phi(- (structural_cluster_obs_v_nm - v_mix_mean[mix_index - 1]) / v_mix_sd[mix_index - 1])';
      
      // structural_obs_takeup_prob[mix_index] = structural_cluster_takeup_prob[mix_index, obs_cluster_id]; // BUGBUG no name matching yet 
    }
  }
}

model {
  beta_control[1] ~ normal(0, beta_control_sd);
  
  if (use_cost_model == COST_MODEL_TYPE_DISCRETE) {
    beta_control[2:] ~ normal(0, beta_far_effect_sd);
  }
  
  beta_ink_effect ~ normal(0, beta_ink_effect_sd);
  beta_calendar_effect ~ normal(0, beta_calendar_effect_sd);
  beta_bracelet_effect ~ normal(0, beta_bracelet_effect_sd);
  
  beta_salience ~ normal(0, 1);
  dist_beta_salience ~ normal(0, 1);
  dist_quadratic_beta_salience ~ normal(0, 1);
  
  v_mu ~ normal(0, 1);
  
  // beta_salience_nm_effect ~ normal(0, 1);
  
  // if (num_treatments_name_matched > 0) {
  //   beta_nm_effect ~ normal(0, 1);
  // }

  if (!suppress_reputation) { 
    mu_rep_raw ~ normal(0, mu_rep_sd);
    
    if (use_mu_cluster_effects) {
      to_vector(mu_cluster_effects_raw) ~ normal(0, 1);
      mu_cluster_effects_sd ~ normal(0, 0.05);
    }
    
    if (use_mu_county_effects) {
      to_vector(mu_county_effects_raw) ~ normal(0, 1);
      mu_county_effects_sd ~ normal(0, 0.1);
    }
  }

  if (num_treatments_param_kappa > 0) {
    dist_cost_k_raw ~ std_normal();
  } 
  
  u_splines_v_sigma ~ normal(0, u_splines_v_sigma_sd);
  
  if (num_treatments_param > 0) {
    dist_beta_v ~ normal(0, 1);
    
    if (use_param_dist_cluster_effects) {
      to_vector(dist_beta_cluster_raw) ~ normal(0, 1);
      dist_beta_cluster_sd ~ normal(0, 0.25);
    }
    
    if (use_param_dist_county_effects) {
      to_vector(dist_beta_county_raw) ~ normal(0, 1);
      dist_beta_county_sd ~ normal(0, 0.25);
    }
  }
  
  if (num_treatments_semiparam > 0) {
    to_vector(u_splines_v_raw) ~ normal(0, 1);
  }
  
  if (num_treatments_param_quadratic > 0) {
    dist_quadratic_beta_v ~ normal(0, 1);
  }
  
  if (use_cluster_effects) {
    to_vector(structural_beta_cluster_raw) ~ std_normal();
    structural_beta_cluster_sd ~ normal(0, structural_beta_cluster_sd_sd);
    
    // if (use_name_matched_obs) {
    //   to_vector(beta_nm_effect_cluster_raw) ~ normal(0, 1);
    //   beta_nm_effect_cluster_sd ~ normal(0, 0.25);
    // }
  }
  
  if (use_county_effects) {
    to_vector(structural_beta_county_raw) ~ std_normal();
    structural_beta_county_sd ~ normal(0, structural_beta_county_sd_sd);
    
    // if (use_name_matched_obs) {
    //   to_vector(beta_nm_effect_county_raw) ~ normal(0, 1);
    //   beta_nm_effect_county_sd ~ normal(0, 0.25);
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
  
  for (dist_group_index in 1:num_discrete_dist) { 
    for (dist_group_mix_index in 1:num_dist_group_mix) {
      group_dist_mean[dist_group_index, dist_group_mix_index] ~ normal(0, 1) T[0, ];
    }
    
    group_dist_sd[dist_group_index] ~ normal(0, 1);
  }
  
  // Distance Likelihood
  
  for (cluster_index in 1:num_clusters) {
    int dist_group_pos = 1;
    int curr_assigned_dist_group = cluster_treatment_map[cluster_assigned_dist_group_treatment[cluster_index], 2];
    
    vector[num_dist_group_mix] group_dist_mix_lp;
    
    for (group_dist_mix_index in 1:num_dist_group_mix) {
      group_dist_mix_lp[group_dist_mix_index] = 
        normal_lpdf(cluster_standard_dist[cluster_index] | group_dist_mean[curr_assigned_dist_group, group_dist_mix_index], group_dist_sd[curr_assigned_dist_group, group_dist_mix_index])  
        + normal_lccdf(0 | group_dist_mean[curr_assigned_dist_group, group_dist_mix_index], group_dist_sd[curr_assigned_dist_group, group_dist_mix_index]) +
        + log(group_dist_mix[curr_assigned_dist_group, group_dist_mix_index]); 
    }
    
    target += log_sum_exp(group_dist_mix_lp);
    
    for (dist_group_index in 1:num_discrete_dist) {
      if (dist_group_index != curr_assigned_dist_group) {
        
        vector[num_dist_group_mix] missing_group_dist_mix_lp;
        
        for (group_dist_mix_index in 1:num_dist_group_mix) {
          missing_group_dist_mix_lp[group_dist_mix_index] =
            normal_lpdf(missing_cluster_standard_dist[cluster_index, dist_group_pos] | group_dist_mean[dist_group_index, group_dist_mix_index], group_dist_sd[dist_group_index, group_dist_mix_index]);
            
          if (missing_cluster_standard_dist[cluster_index, dist_group_pos] < 0) {
            missing_group_dist_mix_lp[group_dist_mix_index] += negative_infinity();
          } else {
            missing_group_dist_mix_lp[group_dist_mix_index] += normal_lccdf(0 | group_dist_mean[dist_group_index, group_dist_mix_index], group_dist_sd[dist_group_index, group_dist_mix_index])
              + log(group_dist_mix[dist_group_index, group_dist_mix_index]); 
          }
        }
        
        target += log_sum_exp(missing_group_dist_mix_lp);
        
        dist_group_pos += 1;
      }
    }
  }
  
  if (!predict_prior) {
    // Take-up Likelihood 
    
    if (num_v_mix == 1) {
      if (use_binomial) {
        if (use_name_matched_obs) {
          reject("Cannot use binomial model if using name-matched observations.");
        }
        
        cluster_takeup_count[included_clusters] ~ binomial(cluster_size[included_clusters], structural_cluster_takeup_prob[1, included_clusters]);
      } else {
        takeup[included_monitored_obs] ~ bernoulli(structural_cluster_takeup_prob[1, obs_cluster_id[included_monitored_obs]]);
        
        if (use_name_matched_obs) {
          takeup[included_name_matched_obs] ~ bernoulli(structural_cluster_takeup_prob_nm[1, obs_cluster_id[included_name_matched_obs]]);
        }
      }
    } else {
      if (use_binomial) {
        if (use_name_matched_obs) {
          reject("Cannot use binomial model if using name-matched observations.");
        }
        
        cluster_takeup_count[included_clusters] ~ mixed_binomial(lambda_v_mix, cluster_size[included_clusters], structural_cluster_takeup_prob[, included_clusters]);
      } else {
        takeup[included_monitored_obs] ~ mixed_bernoulli(lambda_v_mix, structural_cluster_takeup_prob[, obs_cluster_id[included_monitored_obs]]);
        
        if (use_name_matched_obs) {
          takeup[included_name_matched_obs] ~ mixed_bernoulli(lambda_v_mix, structural_cluster_takeup_prob_nm[, obs_cluster_id[included_name_matched_obs]]);
        }
      }
    }
  }
}

generated quantities { 
  matrix[num_shocks, num_shocks] all_u_corr;
 
  vector[num_clusters] cluster_cf_benefit_cost[num_dist_group_treatments]; 
  
  vector[generate_rep ? num_clusters : 0] cluster_rep_benefit_cost; 
  matrix[generate_sim && use_cost_model != COST_MODEL_TYPE_DISCRETE ? num_grid_obs : 0, num_treatments] sim_benefit_cost; 
  
  // Cross Validation
  vector[use_binomial || cluster_log_lik ? num_included_clusters : num_included_obs] log_lik = 
    rep_vector(negative_infinity(), use_binomial || cluster_log_lik ? num_included_clusters : num_included_obs); 
  vector[use_binomial || cluster_log_lik ? num_excluded_clusters : num_excluded_obs] log_lik_heldout = 
    rep_vector(negative_infinity(), use_binomial || cluster_log_lik ? num_excluded_clusters : num_excluded_obs); 
    
  if (num_shocks > 0) {
    all_u_corr = L_all_u_corr * L_all_u_corr';
  }
  
  {
    int treatment_cluster_pos = 1;
    int cluster_treatment_cf_pos = 1;
    
    for (treatment_index in 1:num_dist_group_treatments) {
      vector[num_clusters] treatment_dist_cost;
     
      if (use_cost_model == COST_MODEL_TYPE_DISCRETE) {
        treatment_dist_cost = rep_vector(0, num_clusters);
      } else {                                                            
        int curr_assigned_treatment = cluster_treatment_map[treatment_index, 1];
        int curr_assigned_dist_group = cluster_treatment_map[treatment_index, 2];
        
        treatment_dist_cost = param_dist_cost(cf_cluster_standard_dist[, curr_assigned_dist_group], // cluster_standard_dist, 
                                              cluster_linear_dist_cost[, treatment_index],
                                              rep_vector(quadratic_dist_cost[treatment_index], num_clusters),
                                              u_splines_v[rep_array(curr_assigned_treatment, num_clusters)],
                                              Z_splines_v[rep_array(curr_assigned_treatment, num_clusters)]);
      }
                                                                   
      cluster_cf_benefit_cost[treatment_index] =
        // Not using "restricted" design matrix because restrictions are only on the top-level parameters not village and county level effects
        (structural_beta_cluster + structural_beta_county[cluster_county_id]) * treatment_map_design_matrix[treatment_index]'
        + rep_vector(structural_treatment_effect[treatment_index], num_clusters) 
        - treatment_dist_cost;
    }
  }
  
  // Cross Validation 
  
  if (use_binomial) {
    for (cluster_index_index in 1:num_included_clusters) {
      int cluster_index = included_clusters[cluster_index_index];
      
      if (num_v_mix == 1) {
        log_lik[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], structural_cluster_takeup_prob[1, cluster_index]);
      } else {
        log_lik[cluster_index_index] = mixed_binomial_lpmf({ cluster_takeup_count[cluster_index] } | lambda_v_mix, { cluster_size[cluster_index] }, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
    }
    
    for (cluster_index_index in 1:num_excluded_clusters) {
      int cluster_index = excluded_clusters[cluster_index_index];
      
      if (num_v_mix == 1) {
        log_lik_heldout[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], structural_cluster_takeup_prob[1, cluster_index]);
      } else {
        log_lik_heldout[cluster_index_index] = mixed_binomial_lpmf({ cluster_takeup_count[cluster_index] } | lambda_v_mix, { cluster_size[cluster_index] }, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
    }
  } else {
    vector[num_clusters] temp_log_lik = rep_vector(0, num_clusters);
    
    for (obs_index_index in 1:num_included_monitored_obs) {
      int obs_index = included_monitored_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      real curr_log_lik;
      
      if (num_v_mix == 1) {
        curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[1, cluster_index:cluster_index]);
      } else {
        curr_log_lik = mixed_bernoulli_lpmf({ takeup[obs_index] } | lambda_v_mix, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
      
      if (cluster_log_lik) {
        temp_log_lik[cluster_index] += curr_log_lik;
      } else {
        log_lik[obs_index_index] = curr_log_lik;
      }
    }
    
    for (obs_index_index in 1:num_included_name_matched_obs) {
      int obs_index = included_name_matched_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      real curr_log_lik;
      
      if (num_v_mix == 1) {
        curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob_nm[1, cluster_index:cluster_index]);
      } else {
        curr_log_lik = mixed_bernoulli_lpmf({ takeup[obs_index] } | lambda_v_mix, structural_cluster_takeup_prob_nm[, cluster_index:cluster_index]);
      }
      
      if (cluster_log_lik) {
        temp_log_lik[cluster_index] += curr_log_lik;
      } else {
        log_lik[obs_index_index] = curr_log_lik;
      }
    }
    
    for (obs_index_index in 1:num_excluded_monitored_obs) {
      int obs_index = excluded_monitored_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      real curr_log_lik;
      
      if (num_v_mix == 1) {
        curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[1, cluster_index:cluster_index]);
      } else {
        curr_log_lik = mixed_bernoulli_lpmf({ takeup[obs_index] } | lambda_v_mix, structural_cluster_takeup_prob[, cluster_index:cluster_index]);
      }
      
      if (cluster_log_lik) {
        temp_log_lik[cluster_index] += curr_log_lik;
      } else {
        log_lik_heldout[obs_index_index] = curr_log_lik;
      }
    }
    
    for (obs_index_index in 1:num_excluded_name_matched_obs) {
      int obs_index = excluded_name_matched_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      real curr_log_lik;
      
      if (num_v_mix == 1) {
        curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob_nm[1, cluster_index:cluster_index]);
      } else {
        curr_log_lik = mixed_bernoulli_lpmf({ takeup[obs_index] } | lambda_v_mix, structural_cluster_takeup_prob_nm[, cluster_index:cluster_index]);
      }
      
      if (cluster_log_lik) {
        temp_log_lik[cluster_index] += curr_log_lik;
      } else {
        log_lik_heldout[obs_index_index] = curr_log_lik;
      }
    }
    
    if (cluster_log_lik) {
      log_lik = temp_log_lik[included_clusters];
      log_lik_heldout = temp_log_lik[excluded_clusters];
    }
  }
 
  // Simulation 
  
  if (generate_sim && use_cost_model != COST_MODEL_TYPE_DISCRETE) {
    vector[num_treatments] sim_beta_cluster = use_cluster_effects ? to_vector(normal_rng(rep_array(0, num_treatments), structural_beta_cluster_sd')) : rep_vector(0, num_treatments);
    vector[num_treatments] sim_beta_county = use_county_effects ? to_vector(normal_rng(rep_array(0, num_treatments), structural_beta_county_sd')) : rep_vector(0, num_treatments);
      
    if (use_cost_model == COST_MODEL_TYPE_PARAM_KAPPA) { 
      reject("Not supported yet.");
    } else {
      if (in_array(use_cost_model, { COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE, COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE, COST_MODEL_TYPE_SEMIPARAM_SALIENCE })) {
        reject("Not supported yet.");
      } else {
        if (use_single_cost_model) {
          if (use_param_dist_cluster_effects) {
            reject("Not supported yet.");
          }
          
          if (use_param_dist_county_effects) {
            reject("Not supported yet.");
          }
        } else {
          if (use_param_dist_cluster_effects) {
            reject("Not supported yet.");
          }
          
          if (use_param_dist_county_effects) {
            reject("Not supported yet.");
          }
        }
      } 
    }
    
    for (treatment_index in 1:num_treatments) {
      sim_benefit_cost[, treatment_index] = 
        structural_treatment_effect[treatment_index] +
        treatment_map_design_matrix[treatment_index, :num_treatments] * (sim_beta_cluster + sim_beta_county) +
        param_dist_cost(grid_dist,
                        [ linear_dist_cost[treatment_index] ]',
                        [ quadratic_dist_cost[treatment_index] ]',
                        rep_matrix(u_splines_v[treatment_index], num_grid_obs),
                        Z_grid_v);
    }
  }
  
  if (generate_rep) {
    for (cluster_index in 1:num_clusters) {
      vector[num_dist_group_treatments] rep_beta_cluster = rep_vector(0, num_dist_group_treatments);
        
      vector[num_dist_group_treatments] rep_beta_county = rep_vector(0, num_dist_group_treatments);
      
      int curr_assigned_treatment = cluster_treatment_map[cluster_assigned_dist_group_treatment[cluster_index], 1];
      
      if (use_cost_model != COST_MODEL_TYPE_DISCRETE) {
        real rep_dist_cost;
        
        rep_beta_cluster[1:num_treatments] = use_cluster_effects ? to_vector(normal_rng(rep_array(0, num_treatments), structural_beta_cluster_sd')) : rep_vector(0, num_treatments);
        rep_beta_cluster[(num_treatments + 2):num_dist_group_treatments] = rep_beta_cluster[2:num_treatments];
        
        rep_beta_county[1:num_treatments] = use_county_effects ? to_vector(normal_rng(rep_array(0, num_treatments), structural_beta_county_sd')) : rep_vector(0, num_treatments);
        rep_beta_county[(num_treatments + 2):num_dist_group_treatments] = rep_beta_county[2:num_treatments];
        
        if (use_cost_model == COST_MODEL_TYPE_PARAM_KAPPA) { 
          rep_dist_cost = param_kappa_dist_cost([ cluster_standard_dist[cluster_index] ]', [ dist_cost_k[curr_assigned_treatment] ]')[1];
        } else {
         real rep_linear_dist_cost = linear_dist_cost[cluster_assigned_dist_group_treatment[cluster_index]];
          
          if (in_array(use_cost_model, { COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE, COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE, COST_MODEL_TYPE_SEMIPARAM_SALIENCE })) {
            if (use_param_dist_cluster_effects) {
              rep_linear_dist_cost += normal_rng(0, dist_beta_cluster_sd[1]); 
            }
            
            if (use_param_dist_county_effects) {
              rep_linear_dist_cost += normal_rng(0, dist_beta_county_sd[1]); 
            }
          } else {
            if (use_single_cost_model) {
              if (use_param_dist_cluster_effects) {
                rep_linear_dist_cost += normal_rng(0, dist_beta_cluster_sd[1]); 
              }
              
              if (use_param_dist_county_effects) {
                rep_linear_dist_cost += normal_rng(0, dist_beta_county_sd[1]); 
              }
            } else {
              vector[num_treatments] rep_dist_beta = rep_vector(0, num_treatments);
              
              if (use_param_dist_cluster_effects) {
                rep_dist_beta += to_vector(normal_rng(rep_array(0, num_treatments), dist_beta_cluster_sd'));
              }
              
              if (use_param_dist_county_effects) {
                rep_dist_beta += to_vector(normal_rng(rep_array(0, num_treatments), dist_beta_county_sd'));
              }
              
              rep_linear_dist_cost += treatment_map_design_matrix[cluster_assigned_dist_group_treatment[cluster_index]] * append_row(rep_dist_beta, rep_dist_beta);
            }
          } 
            
          rep_dist_cost = param_dist_cost([ cluster_standard_dist[cluster_index] ]', 
                                          [ rep_linear_dist_cost ]',
                                          [ quadratic_dist_cost[cluster_assigned_dist_group_treatment[cluster_index]] ]',
                                          to_matrix(u_splines_v[curr_assigned_treatment]),
                                          to_matrix(Z_splines_v[curr_assigned_treatment]))[1];
        }
        
        cluster_rep_benefit_cost[cluster_index] = structural_treatment_effect[cluster_assigned_dist_group_treatment[cluster_index]] 
          + treatment_map_design_matrix[cluster_assigned_dist_group_treatment[cluster_index]] * (rep_beta_cluster + rep_beta_county) - rep_dist_cost;
      } else {
        rep_beta_cluster = use_cluster_effects ? to_vector(normal_rng(rep_array(0, num_dist_group_treatments), structural_beta_cluster_sd')) : rep_vector(0, num_dist_group_treatments);
        rep_beta_county = use_county_effects ? to_vector(normal_rng(rep_array(0, num_dist_group_treatments), structural_beta_county_sd')) : rep_vector(0, num_dist_group_treatments);
        
        cluster_rep_benefit_cost[cluster_index] = structural_treatment_effect[cluster_assigned_dist_group_treatment[cluster_index]]
          + treatment_map_design_matrix[cluster_assigned_dist_group_treatment[cluster_index]] * (rep_beta_cluster + rep_beta_county);
      }
    }
  }
}
