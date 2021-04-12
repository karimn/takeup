functions {
#include /../multilvlr/util.stan
#include beliefs_functions.stan
#include takeup_functions.stan
}

data {
#include base_data_sec.stan
#include takeup_data_sec.stan
#include wtp_data.stan
#include beliefs_data_sec.stan
#include dist_data_sec.stan

  int MIN_COST_MODEL_TYPE_VALUE;
  int MAX_COST_MODEL_TYPE_VALUE;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_KAPPA;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_LINEAR;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_QUADRATIC;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_SEMIPARAM;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_SEMIPARAM_SALIENCE;
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_DISCRETE;
  
  int<lower = 1> num_v_mix;

  int<lower = 0, upper = 1> use_wtp_model;
  int<lower = 0, upper = 1> use_homoskedastic_shocks;
  int<lower = 0, upper = 1> use_restricted_mu;
  
  int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> use_cost_model;
  int<lower = 0, upper = 1> suppress_reputation;
  int<lower = 0, upper = 1> use_u_in_delta;
  
  int<lower = 0, upper = 1> use_mu_cluster_effects;
  int<lower = 0, upper = 1> use_mu_county_effects;
  int<lower = 0, upper = 1> use_param_dist_cluster_effects; // These are used for parameteric (linear, quadratic) distance cost models only
  int<lower = 0, upper = 1> use_param_dist_county_effects;
  
  // Rate of Change
  
  int<lower = 1, upper = num_treatments> roc_compare_treatment_id_left;
  int<lower = 1, upper = num_treatments> roc_compare_treatment_id_right;
  int<lower = 0> num_roc_distances;
  vector<lower = 0>[num_roc_distances] roc_distances; // These should be standardized already
 
  // Sim Delta
  
  int<lower = 0> num_sim_delta_w;
  vector[num_sim_delta_w] sim_delta_w;
  
  // Hyperparam
  
  real<lower = 0> mu_rep_sd;
  // real<lower = 0> mu_beliefs_effects_sd;
  // real<lower = 0> mu_beliefs_effects_lambda;
}

transformed data {
#include base_transformed_data_declare.stan 
#include wtp_transformed_data.stan
#include takeup_transformed_data_declare.stan
#include beliefs_transformed_data_declare.stan

#include takeup_transformed_data_define.stan

  real dummy_xr[1] = { 1.0 }; 
  int dummy_xi[1] = { 1 }; 
  
  int<lower = 0, upper = num_treatments> num_treatment_shocks = num_treatments - (use_wtp_model ? 1 : 0);
  
  int<lower = 0, upper = 1> use_dist_salience = in_array(use_cost_model, 
                                                         { COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE, COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE, COST_MODEL_TYPE_SEMIPARAM_SALIENCE });
  int<lower = 0, upper = 1> use_semiparam = in_array(use_cost_model, { COST_MODEL_TYPE_SEMIPARAM, COST_MODEL_TYPE_SEMIPARAM_SALIENCE }); 
  
  int num_treatments_param_kappa = use_cost_model == COST_MODEL_TYPE_PARAM_KAPPA ? num_treatments : 0;
  
  int num_treatments_param = 0;
  int num_treatments_param_quadratic = 0; 
  int num_treatments_semiparam = 0;
  
  if (use_single_cost_model || in_array(use_cost_model, { COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE, 
                                                          COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE, 
                                                          COST_MODEL_TYPE_SEMIPARAM_SALIENCE })) {
    num_treatments_param = 1;
  } else if (in_array(use_cost_model, { COST_MODEL_TYPE_SEMIPARAM, COST_MODEL_TYPE_PARAM_LINEAR, COST_MODEL_TYPE_PARAM_QUADRATIC })) {
    num_treatments_param = num_treatments;
  } 
  
  if ((use_single_cost_model && use_cost_model == COST_MODEL_TYPE_PARAM_QUADRATIC) || use_cost_model == COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE) {
    num_treatments_param_quadratic = 1;
  } else if (in_array(use_cost_model, { COST_MODEL_TYPE_PARAM_QUADRATIC })) {
    num_treatments_param_quadratic = num_treatments;
  } 
  
  if (use_semiparam) {
    if (use_single_cost_model || use_cost_model == COST_MODEL_TYPE_SEMIPARAM_SALIENCE) {
      num_treatments_semiparam = 1;
    } else {
      num_treatments_semiparam = num_treatments;
    }
  }
}

parameters {
#include dist_parameters_sec.stan
#include wtp_parameters.stan
#include beliefs_parameters_sec.stan
  
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
  
  real<lower = 0> base_mu_rep;
  // real<lower = 0> mu_beliefs_effect;
  vector<lower = 0>[use_homoskedastic_shocks ? 1 : num_treatment_shocks] raw_u_sd;
  
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
  
  // WTP valuation parameters
  real<lower = 0> wtp_value_utility;
}

transformed parameters {
#include wtp_transformed_parameters.stan
#include beliefs_transformed_parameters_declare.stan
  
  vector[num_dist_group_treatments] beta;
  
  vector[num_dist_group_treatments] structural_treatment_effect;
  vector[num_clusters] structural_cluster_benefit_cost;
  matrix[num_clusters, num_dist_group_treatments] structural_beta_cluster = rep_matrix(0, num_clusters, num_dist_group_treatments);
  matrix[num_counties, num_dist_group_treatments] structural_beta_county = rep_matrix(0, num_counties, num_dist_group_treatments);
 
  // row_vector[suppress_reputation && !use_dist_salience ? 0 : num_dist_group_treatments - 1] mu_rep_effect;
  // row_vector<lower = 0>[suppress_reputation && !use_dist_salience ? 0 : num_dist_group_treatments] mu_rep;
  // matrix<lower = 0>[!suppress_reputation || use_dist_salience ? num_clusters : 0, num_dist_group_treatments] cluster_mu_rep;
  vector<lower = 0>[!suppress_reputation || use_dist_salience ? num_clusters : 0] obs_cluster_mu_rep;
  
  vector[num_clusters] structural_cluster_obs_v = rep_vector(0, num_clusters);
  row_vector<lower = 0, upper = 1>[fit_model_to_data ? num_clusters : 0] structural_cluster_takeup_prob;
  
  vector<lower = 0>[num_treatments_param_kappa] dist_cost_k;
  vector[num_dist_group_treatments] linear_dist_cost = rep_vector(0, num_dist_group_treatments);
  vector[num_dist_group_treatments] quadratic_dist_cost = rep_vector(0, num_dist_group_treatments);
  matrix[num_clusters, num_dist_group_treatments] cluster_linear_dist_cost = rep_matrix(0, num_clusters, num_dist_group_treatments);
  matrix[num_clusters, num_dist_group_treatments] cluster_quadratic_dist_cost = rep_matrix(0, num_clusters, num_dist_group_treatments);
  matrix[num_treatments, num_knots_v] u_splines_v = rep_matrix(0, num_treatments, num_knots_v);
  vector[num_clusters] cluster_dist_cost = rep_vector(0, num_clusters);
  
  vector<lower = 0>[num_treatments] u_sd;
  vector<lower = 0>[num_treatments] total_error_sd;
  
#include beliefs_transformed_parameters_define.stan
  
  if (use_homoskedastic_shocks) {
    u_sd = rep_vector(use_wtp_model ? sqrt(square(raw_u_sd[1]) + square(wtp_sigma * wtp_value_utility)) : raw_u_sd[1], num_treatments);  
  } else {
    u_sd[{ 1, 2, BRACELET_TREATMENT_INDEX }] = raw_u_sd[{ 1, 2, use_wtp_model ? num_treatment_shocks : BRACELET_TREATMENT_INDEX }];
    u_sd[CALENDAR_TREATMENT_INDEX] = use_wtp_model ? raw_u_sd[num_treatment_shocks] : raw_u_sd[CALENDAR_TREATMENT_INDEX];  
    // u_sd[CALENDAR_TREATMENT_INDEX] = use_wtp_model ? sqrt(square(raw_u_sd[num_treatment_shocks]) + square(wtp_sigma * wtp_value_utility)) : raw_u_sd[CALENDAR_TREATMENT_INDEX];  
  }
    
  total_error_sd = sqrt(1 + square(u_sd));
  
  for (dist_index in 1:num_discrete_dist) {
    if (dist_index > 1) {
      beta[(num_treatments + 1):] = rep_vector(0, num_treatments); 
    } else if (use_wtp_model) { 
      beta[1:2] = [ beta_control, beta_ink_effect ]';
      beta[CALENDAR_TREATMENT_INDEX] = beta_bracelet_effect + wtp_value_utility * hyper_wtp_mu;
      beta[BRACELET_TREATMENT_INDEX] = beta_bracelet_effect;
    } else {
      beta[1:num_treatments] = [ beta_control, beta_ink_effect, beta_calendar_effect, beta_bracelet_effect ]';
    }
  }
 
  structural_treatment_effect = restricted_treatment_map_design_matrix * beta;
  
  // Levels: control ink calendar bracelet
 
  if (!suppress_reputation || use_dist_salience) { 
    obs_cluster_mu_rep = calculate_mu_rep(
      cluster_incentive_treatment_id, cluster_standard_dist, 
      // base_mu_rep, mu_beliefs_effect, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord);
      base_mu_rep, 1, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord);
    
    if (use_mu_cluster_effects) {
      // I shouldn't be using this anymore since we're using beliefs 
      
      // matrix[num_clusters, num_treatments] mu_cluster_effects =  mu_cluster_effects_raw .* rep_matrix(mu_cluster_effects_sd, num_clusters);
      // 
      // cluster_mu_rep = cluster_mu_rep .* exp(mu_cluster_effects);
    }
    
    if (use_mu_county_effects) {
      // matrix[num_counties, num_treatments] mu_county_effects =  mu_county_effects_raw .* rep_matrix(mu_county_effects_sd, num_counties);
      // 
      // cluster_mu_rep = cluster_mu_rep .* exp(mu_county_effects[cluster_county_id]);
    } 
  }
  
  if (in_array(use_cost_model, { COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE, COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE, COST_MODEL_TYPE_SEMIPARAM_SALIENCE })) {
    // linear_dist_cost = rep_vector(dist_beta_v[1], num_dist_group_treatments) + dist_beta_salience * mu_rep';  // append_col(mu_rep, mu_rep)';
    reject("Salience not yet supported."); // TODO Fix this
    
    if (use_param_dist_cluster_effects) {
      cluster_linear_dist_cost += rep_matrix(dist_beta_cluster_raw[, 1] * dist_beta_cluster_sd[1], num_dist_group_treatments);
    }
    
    if (use_param_dist_county_effects) {
      cluster_linear_dist_cost += rep_matrix((dist_beta_county_raw[, 1] * dist_beta_county_sd[1])[cluster_county_id], num_dist_group_treatments);
    }
    
    if (use_cost_model == COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE) {
      // quadratic_dist_cost = rep_vector(dist_quadratic_beta_v[1], num_dist_group_treatments) + dist_quadratic_beta_salience * mu_rep'; // append_col(mu_rep, mu_rep)';
      reject("Salience not yet supported."); // TODO Fix this
      
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
        // u_splines_v[treatment_index] = u_splines_v_raw[1] * u_splines_v_sigma * mu_rep[treatment_index]; // BUGBUG This shouldn't work anymore
        reject("Salience not yet supported."); // TODO Fix this
      } else {
        u_splines_v[treatment_index] = u_splines_v_raw[treatment_index] * u_splines_v_sigma;
      }
    }        
  } 
    
  cluster_dist_cost = param_dist_cost_with_splines(
                                      cluster_standard_dist, 
                                      to_vector(cluster_linear_dist_cost)[long_cluster_by_treatment_index],
                                      to_vector(cluster_quadratic_dist_cost)[long_cluster_by_treatment_index],
                                      u_splines_v[cluster_incentive_treatment_id],
                                      Z_splines_v[cluster_incentive_treatment_id]);
  
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
          // to_vector(cluster_mu_rep)[long_cluster_by_treatment_index],
          obs_cluster_mu_rep,
          total_error_sd[cluster_incentive_treatment_id],
          u_sd[cluster_incentive_treatment_id],
          
          use_u_in_delta,
          alg_sol_rel_tol, // 1e-10,
          alg_sol_f_tol, // 1e-5,
          alg_sol_max_steps
        ); 
      } else {
        for (cluster_index in 1:num_clusters) {
          structural_cluster_obs_v[cluster_index] = find_fixedpoint_solution(
            structural_cluster_benefit_cost[cluster_index],
            // cluster_mu_rep[cluster_index, cluster_assigned_dist_group_treatment[cluster_index]],
            obs_cluster_mu_rep[cluster_index],
            total_error_sd[cluster_incentive_treatment_id[cluster_index]],
            u_sd[cluster_incentive_treatment_id[cluster_index]],

            use_u_in_delta,
            alg_sol_rel_tol, // 1e-10,
            alg_sol_f_tol, // 1e-5,
            alg_sol_max_steps
          );
        }
      }
    }
    
    structural_cluster_takeup_prob = Phi_approx(- structural_cluster_obs_v ./ total_error_sd[cluster_incentive_treatment_id])';
  }
}

model {
#include wtp_model_section.stan
#include beliefs_model_sec.stan
#include dist_model_sec.stan
  
  wtp_value_utility ~ normal(0, 0.1);

  beta_control ~ normal(0, beta_control_sd);
  
  beta_ink_effect ~ normal(0, beta_ink_effect_sd);
  beta_calendar_effect ~ normal(0, beta_calendar_effect_sd);
  beta_bracelet_effect ~ normal(0, beta_bracelet_effect_sd);
  
  beta_salience ~ normal(0, 1);
  dist_beta_salience ~ normal(0, 1);
  dist_quadratic_beta_salience ~ normal(0, 1);
  
  base_mu_rep ~ normal(0, mu_rep_sd);
  // mu_beliefs_effect ~ normal(0, mu_beliefs_effects_sd);
  // mu_beliefs_effect ~ exponential(mu_beliefs_effects_lambda);
  
  if (!suppress_reputation || use_dist_salience) { 
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
  
  raw_u_sd ~ normal(0, 1);
  
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
#include beliefs_generated_quantities_declare.stan
#include dist_generated_quantities_declare.stan

  matrix[num_clusters, num_dist_group_treatments] structural_cluster_benefit = 
        rep_matrix(structural_treatment_effect', num_clusters) + 
        (structural_beta_cluster + structural_beta_county[cluster_county_id]) * treatment_map_design_matrix';
  
  // matrix<lower = 0>[num_clusters, num_discrete_dist] cf_cluster_standard_dist; 
 
  vector[num_clusters] cluster_cf_benefit_cost[num_dist_group_treatments]; 
  
  vector[num_clusters] cluster_cf_cutoff[num_dist_group_treatments, num_treatments]; 
  
  vector[generate_rep ? num_clusters : 0] cluster_rep_benefit_cost; 
  vector[generate_rep ? num_clusters : 0] cluster_rep_cutoff; 
  matrix[generate_sim ? num_grid_obs : 0, num_treatments] sim_benefit_cost; 
  
  matrix[num_clusters, num_roc_distances] cluster_roc_left;
  matrix[num_clusters, num_roc_distances] cluster_roc_right;
  matrix[num_clusters, num_roc_distances] cluster_roc_diff;
  matrix[num_clusters, num_roc_distances] cluster_w_cutoff_left;
  matrix[num_clusters, num_roc_distances] cluster_w_cutoff_right;
  matrix[num_clusters, num_roc_distances] cluster_takeup_prop_left;
  matrix[num_clusters, num_roc_distances] cluster_takeup_prop_right;
  matrix[num_clusters, num_roc_distances] cluster_social_multiplier_left;
  matrix[num_clusters, num_roc_distances] cluster_social_multiplier_right;
  matrix[num_clusters, num_roc_distances] cluster_rep_return_left;
  matrix[num_clusters, num_roc_distances] cluster_rep_return_right;
  
  vector[num_sim_delta_w] sim_delta;
  
  // Cross Validation
  vector[cross_validate ? (use_binomial || cluster_log_lik ? num_included_clusters : num_included_obs) : 0] log_lik;
  vector[cross_validate ? (use_binomial || cluster_log_lik ? num_excluded_clusters : num_excluded_obs) : 0] log_lik_heldout;
  
#include dist_generated_quantities_define.stan
#include wtp_generated_quantities.stan
#include beliefs_generated_quantities_define.stan

  {
    int treatment_cluster_pos = 1;
    int cluster_treatment_cf_pos = 1;
    
    for (treatment_index in 1:num_dist_group_treatments) {
      vector[num_clusters] treatment_dist_cost;
     
      int curr_assigned_treatment = cluster_treatment_map[treatment_index, 1];
      int curr_assigned_dist_group = cluster_treatment_map[treatment_index, 2];
      
      treatment_dist_cost = param_dist_cost_with_splines(
                                            all_cluster_standard_dist[, curr_assigned_dist_group],
                                            cluster_linear_dist_cost[, treatment_index],
                                            cluster_quadratic_dist_cost[, treatment_index],
                                            u_splines_v[rep_array(curr_assigned_treatment, num_clusters)],
                                            Z_splines_v[rep_array(curr_assigned_treatment, num_clusters)]);
                                                                 
      cluster_cf_benefit_cost[treatment_index] = structural_cluster_benefit[, treatment_index] - treatment_dist_cost;
      
      for (mu_treatment_index in 1:num_treatments) {
        vector[num_clusters] curr_cluster_mu_rep = calculate_mu_rep(
          { mu_treatment_index }, all_cluster_standard_dist[, curr_assigned_dist_group], 
          // base_mu_rep, mu_beliefs_effect, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord);
          base_mu_rep, 1, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord);
        
        if (multithreaded) {
          cluster_cf_cutoff[treatment_index, mu_treatment_index] = map_find_fixedpoint_solution(
            cluster_cf_benefit_cost[treatment_index],
            curr_cluster_mu_rep,
            rep_vector(total_error_sd[curr_assigned_treatment], num_clusters),
            rep_vector(u_sd[curr_assigned_treatment], num_clusters),
            
            use_u_in_delta, 
            alg_sol_rel_tol, 
            alg_sol_f_tol, 
            alg_sol_max_steps
          );
        } else {
          for (cluster_index in 1:num_clusters) {
            cluster_cf_cutoff[treatment_index, mu_treatment_index, cluster_index] = find_fixedpoint_solution(
              cluster_cf_benefit_cost[treatment_index, cluster_index],
              curr_cluster_mu_rep[cluster_index],
              total_error_sd[curr_assigned_treatment],
              u_sd[curr_assigned_treatment],
              
              use_u_in_delta, 
              alg_sol_rel_tol, 
              alg_sol_f_tol, 
              alg_sol_max_steps
            ); 
          }
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
        param_dist_cost_with_splines(
                        grid_dist,
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
      
      int curr_assigned_treatment = cluster_incentive_treatment_id[cluster_index];
      
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
        
      rep_dist_cost = param_dist_cost_with_splines(
                                      [ cluster_standard_dist[cluster_index] ]', 
                                      [ rep_linear_dist_cost ]',
                                      [ rep_quadratic_dist_cost ]', 
                                      to_matrix(u_splines_v[curr_assigned_treatment]),
                                      to_matrix(Z_splines_v[curr_assigned_treatment]))[1];
      
      cluster_rep_benefit_cost[cluster_index] = structural_treatment_effect[cluster_assigned_dist_group_treatment[cluster_index]] 
        + treatment_map_design_matrix[cluster_assigned_dist_group_treatment[cluster_index]] * (rep_beta_cluster + rep_beta_county) - rep_dist_cost;
        
      cluster_rep_cutoff[cluster_index] = find_fixedpoint_solution(
          cluster_rep_benefit_cost[cluster_index],
          obs_cluster_mu_rep[cluster_index],
          total_error_sd[curr_assigned_treatment],
          u_sd[curr_assigned_treatment],
          
          use_u_in_delta,
          alg_sol_rel_tol,
          alg_sol_f_tol,
          alg_sol_max_steps
        );
    }
  }
 
  if (multithreaded) { 
    for (roc_dist_index in 1:num_roc_distances) {
      vector[num_clusters] roc_cluster_dist = rep_vector(roc_distances[roc_dist_index], num_clusters);
      
      vector[num_clusters] curr_net_benefit_right = 
        structural_cluster_benefit[, roc_compare_treatment_id_right] - param_dist_cost(roc_cluster_dist,
                                                                                       cluster_linear_dist_cost[, roc_compare_treatment_id_right],
                                                                                       cluster_quadratic_dist_cost[, roc_compare_treatment_id_right]);
                                                                                       
      vector[num_clusters] curr_net_benefit_left = 
        structural_cluster_benefit[, roc_compare_treatment_id_left] - param_dist_cost(roc_cluster_dist,
                                                                                       cluster_linear_dist_cost[, roc_compare_treatment_id_left],
                                                                                       cluster_quadratic_dist_cost[, roc_compare_treatment_id_left]);
                                                                                       
      matrix[num_clusters, 2] curr_cluster_mu_rep_left = calculate_mu_rep_deriv(
        roc_compare_treatment_id_left, roc_cluster_dist,
        base_mu_rep, 1, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord);
        // base_mu_rep, mu_beliefs_effect, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord);
        
      matrix[num_clusters, 2] curr_cluster_mu_rep_right = calculate_mu_rep_deriv(
        roc_compare_treatment_id_right, roc_cluster_dist,
        base_mu_rep, 1, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord);
        // base_mu_rep, mu_beliefs_effect, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord);
       
      matrix[num_clusters, 5] roc_results = map_calculate_roc_diff(
        roc_compare_treatment_id_left, roc_compare_treatment_id_right, 
        curr_net_benefit_left,
        curr_net_benefit_right,
        
        rep_vector(total_error_sd[roc_compare_treatment_id_right], num_clusters),
        rep_vector(u_sd[roc_compare_treatment_id_right], num_clusters),
        
        cluster_linear_dist_cost[, roc_compare_treatment_id_right],
        
        curr_cluster_mu_rep_left[, 1], curr_cluster_mu_rep_right[, 1],
        curr_cluster_mu_rep_left[, 2], curr_cluster_mu_rep_right[, 2],
        
        use_u_in_delta,
        alg_sol_rel_tol,
        alg_sol_f_tol,
        alg_sol_max_steps
      );
        
      cluster_w_cutoff_left[, roc_dist_index] = roc_results[, 1];
      cluster_w_cutoff_right[, roc_dist_index] = roc_results[, 2];
      cluster_takeup_prop_left[, roc_dist_index] = 1 - Phi_approx(roc_results[, 1] / total_error_sd[1]);
      cluster_takeup_prop_right[, roc_dist_index] = 1 - Phi_approx(roc_results[, 2] / total_error_sd[1]);
      cluster_social_multiplier_left[, roc_dist_index] = - roc_results[, 3];
      cluster_social_multiplier_right[, roc_dist_index] = - roc_results[, 5];
      cluster_rep_return_left[, roc_dist_index] = roc_results[, 3] .* curr_cluster_mu_rep_left[, 1]; 
      cluster_rep_return_right[, roc_dist_index] = roc_results[, 5] .* curr_cluster_mu_rep_right[, 1];
      cluster_roc_left[, roc_dist_index] = roc_results[, 7]; 
      cluster_roc_right[, roc_dist_index] = roc_results[, 8]; 
      cluster_roc_diff[, roc_dist_index] = roc_results[, 9]; 
    }
  }
  
  for (sim_delta_index in 1:num_sim_delta_w) {
    sim_delta[sim_delta_index] = expected_delta(sim_delta_w[sim_delta_index], total_error_sd[1], u_sd[1], dummy_xr, dummy_xi);
  }
}

