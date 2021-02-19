functions {
#include takeup_functions.stan
}

data {
#include takeup_data_sec.stan
#include wtp_data.stan

  int<lower = 0, upper = 1> use_wtp_model;
}

transformed data {
#include wtp_transformed_data.stan
#include takeup_transformed_data_declare.stan
#include takeup_transformed_data_define.stan
}

parameters {
#include wtp_parameters.stan
  
  // Levels: control ink calendar bracelet
  real beta_control;
  real beta_ink_effect;
  real<lower = (use_private_incentive_restrictions ? 0 : negative_infinity())> beta_calendar_effect;
  real<lower = (use_private_incentive_restrictions ? 0 : negative_infinity())> beta_bracelet_effect;
  
  matrix[use_cluster_effects ? num_clusters : 0, num_treatments] structural_beta_cluster_raw;
  row_vector<lower = 0>[use_cluster_effects ? num_treatments : 0] structural_beta_cluster_sd;
  
  matrix[use_county_effects ? num_counties : 0, num_treatments] structural_beta_county_raw;
  row_vector<lower = 0>[use_county_effects ? num_treatments : 0] structural_beta_county_sd;
  
  // Salience
  
  real<lower = 0> beta_salience;
  real<lower = 0> dist_beta_salience;
  real<lower = 0> dist_quadratic_beta_salience;
  
  // Reputational Returns
  
  real v_mu;
  
  row_vector<lower = 0>[suppress_reputation && !use_dist_salience ? 0 : num_treatments] mu_rep_raw;
  real<lower = 0> u_sd;
  cholesky_factor_corr[suppress_shocks ? 0 : 2] L_all_u_corr;
  
  matrix[!use_mu_cluster_effects || (suppress_reputation && !use_dist_salience) ? 0 : num_clusters, num_treatments] mu_cluster_effects_raw;
  row_vector<lower = 0>[!use_mu_cluster_effects || (suppress_reputation && !use_dist_salience) ? 0 : num_treatments] mu_cluster_effects_sd;
  
  matrix[!use_mu_county_effects || (suppress_reputation && !use_dist_salience) ? 0 : num_counties, num_treatments] mu_county_effects_raw;
  row_vector<lower = 0>[!use_mu_county_effects || (suppress_reputation && !use_dist_salience) ? 0 : num_treatments] mu_county_effects_sd;
  
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
  
  ordered[num_dist_group_mix] group_dist_mean[num_discrete_dist];
  vector<lower = 0>[num_dist_group_mix] group_dist_sd[num_discrete_dist];
  
  matrix<lower = 0>[num_clusters, num_discrete_dist - 1] missing_cluster_standard_dist; 
  
  // WTP valuation parameters
  real<lower = 0> wtp_value_utility;
}

transformed parameters {
#include wtp_transformed_parameters.stan
  
  vector[num_dist_group_treatments] beta;
  
  vector[num_dist_group_treatments] structural_treatment_effect;
  vector[num_clusters] structural_cluster_benefit_cost;
  matrix[num_clusters, num_dist_group_treatments] structural_beta_cluster = rep_matrix(0, num_clusters, num_dist_group_treatments);
  matrix[num_counties, num_dist_group_treatments] structural_beta_county = rep_matrix(0, num_counties, num_dist_group_treatments);
 
  row_vector<lower = 0>[suppress_reputation && !use_dist_salience ? 0 : num_treatments] mu_rep = mu_rep_raw;
  matrix<lower = 0>[!suppress_reputation || use_dist_salience ? num_clusters : 0, num_treatments] cluster_mu_rep;
  
  vector[num_clusters] structural_cluster_obs_v = rep_vector(0, num_clusters);
  row_vector<lower = 0, upper = 1>[fit_model_to_data ? num_clusters : 0] structural_cluster_takeup_prob;
  
  vector<lower = 0>[num_treatments_param_kappa] dist_cost_k;
  vector[num_dist_group_treatments] linear_dist_cost = rep_vector(0, num_dist_group_treatments);
  vector[num_dist_group_treatments] quadratic_dist_cost = rep_vector(0, num_dist_group_treatments);
  matrix[num_clusters, num_dist_group_treatments] cluster_linear_dist_cost = rep_matrix(0, num_clusters, num_dist_group_treatments);
  matrix[num_clusters, num_dist_group_treatments] cluster_quadratic_dist_cost = rep_matrix(0, num_clusters, num_dist_group_treatments);
  matrix[num_treatments, num_knots_v] u_splines_v = rep_matrix(0, num_treatments, num_knots_v);
  vector[num_clusters] cluster_dist_cost = rep_vector(0, num_clusters);
  
  matrix<lower = 0>[num_clusters, num_discrete_dist] cf_cluster_standard_dist; 
  
  real total_error_sd;
  
  if (!suppress_shocks) {
    vector[2] all_u_sd = [ 1.0, u_sd ]';
    matrix[2, 2] L_all_u_vcov = diag_pre_multiply(all_u_sd, L_all_u_corr);
    matrix[2, 2] all_u_vcov = L_all_u_vcov * L_all_u_vcov'; 
    total_error_sd = sqrt(sum(all_u_vcov) + (use_wtp_model ? square(wtp_sigma * wtp_value_utility) : 0));
  } else {
    total_error_sd = sqrt(1 + (use_wtp_model ? square(wtp_sigma * wtp_value_utility) : 0));
  }
  
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
    if (dist_index > 1) {
      beta[(num_treatments * (dist_index - 1) + 1):(num_treatments * dist_index)] = rep_vector(0, num_treatments); 
    } else if (use_wtp_model) { 
      beta[(num_treatments * (dist_index - 1) + 1):(num_treatments * dist_index)] = 
        [ beta_control, beta_ink_effect, beta_bracelet_effect + wtp_value_utility * hyper_wtp_mu, beta_bracelet_effect ]';
    } else {
      beta[(num_treatments * (dist_index - 1) + 1):(num_treatments * dist_index)] = 
        [ beta_control, beta_ink_effect, beta_calendar_effect, beta_bracelet_effect ]';
    }
  }
 
  structural_treatment_effect = restricted_treatment_map_design_matrix * beta;
  
  // Levels: control ink calendar bracelet
 
  if (!suppress_reputation || use_dist_salience) { 
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
      
      if (use_param_dist_cluster_effects) {
        // TODO
      }
      
      if (use_param_dist_county_effects) {
        // TODO
      }
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
        
        if (use_param_dist_cluster_effects) {
          // TODO
        }
        
        if (use_param_dist_county_effects) {
          // TODO
        }
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
        
        if (use_param_dist_cluster_effects) {
          // TODO
        }
        
        if (use_param_dist_county_effects) {
          // TODO
        }
      }
    }
  } 
  
  cluster_linear_dist_cost += rep_matrix(linear_dist_cost', num_clusters);
  cluster_quadratic_dist_cost += rep_matrix(quadratic_dist_cost', num_clusters);
  
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
                                      to_vector(cluster_quadratic_dist_cost)[long_cluster_by_treatment_index],
                                      u_splines_v[cluster_treatment_map[cluster_assigned_dist_group_treatment, 1]],
                                      Z_splines_v[cluster_treatment_map[cluster_assigned_dist_group_treatment, 1]]);
  
  structural_cluster_benefit_cost = structural_treatment_effect[cluster_assigned_dist_group_treatment] - cluster_dist_cost;
  
  if (use_cluster_effects) {
    vector[num_clusters] cluster_effects;
    
    structural_beta_cluster[, 1:num_treatments] = structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_clusters);
    
    if (use_wtp_model) { // Bracelet and Calendar are the same
      structural_beta_cluster[, 3] = structural_beta_cluster[, 4];
    } 
    
    cluster_effects = rows_dot_product(cluster_treatment_design_matrix, structural_beta_cluster); 
    structural_cluster_benefit_cost += cluster_effects;
  }
  
  if (use_county_effects) {
    vector[num_clusters] county_effects;
    
    structural_beta_county[, 1:num_treatments] = structural_beta_county_raw .* rep_matrix(structural_beta_county_sd, num_counties);
    
    if (use_wtp_model) { // Calendar = Bracelet + strata_effect
      structural_beta_county[, 3] = structural_beta_county[, 4] + wtp_value_utility * strata_effect_wtp_mu; 
    } 
    
    county_effects = rows_dot_product(cluster_treatment_design_matrix, structural_beta_county[cluster_county_id]); 
    structural_cluster_benefit_cost += county_effects;
  }

  if (fit_model_to_data) {
    if (suppress_reputation) {
      structural_cluster_obs_v = - structural_cluster_benefit_cost;
    } else {
      if (multithreaded) {
        structural_cluster_obs_v = map_find_fixedpoint_solution(
          structural_cluster_benefit_cost, 
          to_vector(cluster_mu_rep)[long_cluster_by_incentive_treatment_index],
          total_error_sd,
          u_sd,
          
          use_u_in_delta,
          alg_sol_rel_tol, // 1e-10,
          alg_sol_f_tol, // 1e-5,
          alg_sol_max_steps
        ); 
      } else {
        for (cluster_index in 1:num_clusters) {
          structural_cluster_obs_v[cluster_index] = find_fixedpoint_solution(
            structural_cluster_benefit_cost[cluster_index],
            cluster_mu_rep[cluster_index, cluster_treatment_map[cluster_assigned_dist_group_treatment, 1]],
            total_error_sd,
            u_sd,
  
            use_u_in_delta,
            alg_sol_rel_tol, // 1e-10,
            alg_sol_f_tol, // 1e-5,
            alg_sol_max_steps
          );
        }
      }
    }
    
    structural_cluster_takeup_prob = Phi_approx(- structural_cluster_obs_v / total_error_sd)';
  }
}

model {
#include wtp_model_section.stan
  
  wtp_value_utility ~ normal(0, 0.1);

  beta_control ~ normal(0, beta_control_sd);
  
  beta_ink_effect ~ normal(0, beta_ink_effect_sd);
  beta_calendar_effect ~ normal(0, beta_calendar_effect_sd);
  beta_bracelet_effect ~ normal(0, beta_bracelet_effect_sd);
  
  beta_salience ~ normal(0, 1);
  dist_beta_salience ~ normal(0, 1);
  dist_quadratic_beta_salience ~ normal(0, 1);
  
  v_mu ~ normal(0, 1);
  
  if (!suppress_reputation || use_dist_salience) { 
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

  u_splines_v_sigma ~ normal(0, u_splines_v_sigma_sd);
  
  if (num_treatments_param > 0) {
    dist_beta_v ~ normal(0, dist_beta_v_sd);
    
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
  }
  
  if (use_county_effects) {
    to_vector(structural_beta_county_raw) ~ std_normal();
    structural_beta_county_sd ~ normal(0, structural_beta_county_sd_sd);
  }
  
  u_sd ~ normal(0, 1);
  L_all_u_corr ~ lkj_corr_cholesky(2);
  
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
  
  if (fit_model_to_data) {
    // Take-up Likelihood 
    if (use_binomial) {
      cluster_takeup_count[included_clusters] ~ binomial(cluster_size[included_clusters], structural_cluster_takeup_prob[included_clusters]);
    } else {
      takeup[included_monitored_obs] ~ bernoulli(structural_cluster_takeup_prob[obs_cluster_id[included_monitored_obs]]);
    }
  }
}

generated quantities {
  matrix[num_clusters, num_dist_group_treatments] structural_cluster_benefit = 
        rep_matrix(structural_treatment_effect', num_clusters) + 
        (structural_beta_cluster + structural_beta_county[cluster_county_id]) * treatment_map_design_matrix';
  
  matrix[2, 2] all_u_corr;
 
  vector[num_clusters] cluster_cf_benefit_cost[num_dist_group_treatments]; 
  
  matrix[num_dist_group_treatments, num_treatments] cluster_cf_cutoff[num_clusters]; 
  
  vector[generate_rep ? num_clusters : 0] cluster_rep_benefit_cost; 
  vector[generate_rep ? num_clusters : 0] cluster_rep_cutoff; 
  matrix[generate_sim ? num_grid_obs : 0, num_treatments] sim_benefit_cost; 
  
  // Cross Validation
  vector[cross_validate ? (use_binomial || cluster_log_lik ? num_included_clusters : num_included_obs) : 0] log_lik;
  vector[cross_validate ? (use_binomial || cluster_log_lik ? num_excluded_clusters : num_excluded_obs) : 0] log_lik_heldout;
  
#include wtp_generated_quantities.stan
    
  if (!suppress_shocks) {
    all_u_corr = L_all_u_corr * L_all_u_corr';
  }
  
  
  {
    int treatment_cluster_pos = 1;
    int cluster_treatment_cf_pos = 1;
    
    for (treatment_index in 1:num_dist_group_treatments) {
      vector[num_clusters] treatment_dist_cost;
     
      int curr_assigned_treatment = cluster_treatment_map[treatment_index, 1];
      int curr_assigned_dist_group = cluster_treatment_map[treatment_index, 2];
      
      treatment_dist_cost = param_dist_cost(cf_cluster_standard_dist[, curr_assigned_dist_group], // cluster_standard_dist, 
                                            cluster_linear_dist_cost[, treatment_index],
                                            cluster_quadratic_dist_cost[, treatment_index],
                                            u_splines_v[rep_array(curr_assigned_treatment, num_clusters)],
                                            Z_splines_v[rep_array(curr_assigned_treatment, num_clusters)]);
                                                                 
      cluster_cf_benefit_cost[treatment_index] = structural_cluster_benefit[, treatment_index] - treatment_dist_cost;
      
      for (cluster_index in 1:num_clusters) {
        for (mu_treatment_index in 1:num_treatments) {
          cluster_cf_cutoff[cluster_index, treatment_index, mu_treatment_index] = find_fixedpoint_solution(
            cluster_cf_benefit_cost[treatment_index, cluster_index],
            cluster_mu_rep[cluster_index, mu_treatment_index],
            total_error_sd,
            u_sd,
            
            use_u_in_delta, 
            alg_sol_rel_tol, 
            alg_sol_f_tol, 
            alg_sol_max_steps
          ); 
        }
      }
    }
  }
  
  // Cross Validation 
  if (cross_validate) {
    log_lik = rep_vector(negative_infinity(), use_binomial || cluster_log_lik ? num_included_clusters : num_included_obs); 
    log_lik_heldout = rep_vector(negative_infinity(), use_binomial || cluster_log_lik ? num_excluded_clusters : num_excluded_obs); 
    
    if (use_binomial) {
      for (cluster_index_index in 1:num_included_clusters) {
        int cluster_index = included_clusters[cluster_index_index];
        
        log_lik[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], structural_cluster_takeup_prob[cluster_index]);
      }
      
      for (cluster_index_index in 1:num_excluded_clusters) {
        int cluster_index = excluded_clusters[cluster_index_index];
        
        log_lik_heldout[cluster_index_index] = binomial_lpmf(cluster_takeup_count[cluster_index] | cluster_size[cluster_index], structural_cluster_takeup_prob[cluster_index]);
      }
    } else {
      vector[num_clusters] temp_log_lik = rep_vector(0, num_clusters);
      
      for (obs_index_index in 1:num_included_monitored_obs) {
        int obs_index = included_monitored_obs[obs_index_index];
        int cluster_index = obs_cluster_id[obs_index];
        
        real curr_log_lik;
        
        curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[cluster_index:cluster_index]);
        
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
        
        curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[cluster_index:cluster_index]);
        
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
  }   
  
  // Simulation 
  if (generate_sim) {
    vector[num_treatments] sim_beta_cluster = use_cluster_effects ? to_vector(normal_rng(rep_array(0, num_treatments), structural_beta_cluster_sd')) : rep_vector(0, num_treatments);
    vector[num_treatments] sim_beta_county = use_county_effects ? to_vector(normal_rng(rep_array(0, num_treatments), structural_beta_county_sd')) : rep_vector(0, num_treatments);
      
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

  // Replicated Data  
  if (generate_rep) {
    for (cluster_index in 1:num_clusters) {
      vector[num_dist_group_treatments] rep_beta_cluster = rep_vector(0, num_dist_group_treatments);
        
      vector[num_dist_group_treatments] rep_beta_county = rep_vector(0, num_dist_group_treatments);
      
      int curr_assigned_treatment = cluster_treatment_map[cluster_assigned_dist_group_treatment[cluster_index], 1];
      
      real rep_dist_cost;
      
      real rep_linear_dist_cost = linear_dist_cost[cluster_assigned_dist_group_treatment[cluster_index]];
      real rep_quadratic_dist_cost = quadratic_dist_cost[cluster_assigned_dist_group_treatment[cluster_index]];
      
      rep_beta_cluster[1:num_treatments] = use_cluster_effects ? to_vector(normal_rng(rep_array(0, num_treatments), structural_beta_cluster_sd')) : rep_vector(0, num_treatments);
      rep_beta_cluster[(num_treatments + 2):num_dist_group_treatments] = rep_beta_cluster[2:num_treatments];
      
      rep_beta_county[1:num_treatments] = use_county_effects ? to_vector(normal_rng(rep_array(0, num_treatments), structural_beta_county_sd')) : rep_vector(0, num_treatments);
      rep_beta_county[(num_treatments + 2):num_dist_group_treatments] = rep_beta_county[2:num_treatments];
        
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
          // TODO rep_quadratic_dist_cost
        }
      } 
        
      rep_dist_cost = param_dist_cost([ cluster_standard_dist[cluster_index] ]', 
                                      [ rep_linear_dist_cost ]',
                                      [ rep_quadratic_dist_cost ]', 
                                      to_matrix(u_splines_v[curr_assigned_treatment]),
                                      to_matrix(Z_splines_v[curr_assigned_treatment]))[1];
      
      cluster_rep_benefit_cost[cluster_index] = structural_treatment_effect[cluster_assigned_dist_group_treatment[cluster_index]] 
        + treatment_map_design_matrix[cluster_assigned_dist_group_treatment[cluster_index]] * (rep_beta_cluster + rep_beta_county) - rep_dist_cost;
        
      cluster_rep_cutoff[cluster_index] = find_fixedpoint_solution(
          cluster_rep_benefit_cost[cluster_index],
          cluster_mu_rep[cluster_index, curr_assigned_treatment],
          total_error_sd,
          u_sd,
          
          use_u_in_delta,
          alg_sol_rel_tol,
          alg_sol_f_tol,
          alg_sol_max_steps
        );
    }
  }
}
