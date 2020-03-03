#include /takeup_header.stan

parameters {
  // Levels: control ink calendar bracelet
  vector[num_discrete_dist] beta_control;
  vector[num_discrete_dist] beta_ink_effect;
  vector[num_discrete_dist] beta_calendar_effect;
  vector[num_discrete_dist] beta_bracelet_effect;
  
  matrix[use_cluster_effects ? num_clusters : 0, num_dist_group_treatments] structural_beta_cluster_raw;
  row_vector<lower = 0>[use_cluster_effects ? num_dist_group_treatments: 0] structural_beta_cluster_sd;
  
  matrix[use_county_effects ? num_counties : 0, num_dist_group_treatments] structural_beta_county_raw;
  row_vector<lower = 0>[use_county_effects ? num_dist_group_treatments : 0] structural_beta_county_sd;
  
  // Name Matched
  
  // vector[num_treatments_name_matched] beta_nm_effect;
  // 
  // matrix[use_cluster_effects && use_name_matched_obs ? num_clusters : 0, (use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_treatments * num_discrete_dist : num_treatments)] beta_nm_effect_cluster_raw;
  // row_vector<lower = 0>[use_cluster_effects && use_name_matched_obs ? (use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_treatments * num_discrete_dist : num_treatments) : 0] beta_nm_effect_cluster_sd;
  // 
  // matrix[use_county_effects && use_name_matched_obs ? num_clusters : 0, (use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_treatments * num_discrete_dist : num_treatments)] beta_nm_effect_county_raw;
  // row_vector<lower = 0>[use_county_effects && use_name_matched_obs ? (use_cost_model == COST_MODEL_TYPE_DISCRETE ? num_treatments * num_discrete_dist : num_treatments) : 0] beta_nm_effect_county_sd;
}

transformed parameters {
  vector[num_dist_group_treatments] beta; 
  vector[num_dist_group_treatments] structural_treatment_effect;
  vector[num_clusters] structural_cluster_benefit_cost;
  matrix[num_clusters, num_dist_group_treatments] structural_beta_cluster = rep_matrix(0, num_clusters, num_dist_group_treatments);
  matrix[num_counties, num_dist_group_treatments] structural_beta_county = rep_matrix(0, num_counties, num_dist_group_treatments);
  
  vector[num_clusters] structural_cluster_obs_v = rep_vector(0, num_clusters);
  vector<lower = 0, upper = 1>[num_clusters] structural_cluster_takeup_prob;
  
  for (dist_index in 1:num_discrete_dist) {
    beta[(num_treatments * (dist_index - 1) + 1):(num_treatments * dist_index)] = [ beta_control[dist_index], beta_ink_effect[dist_index], beta_calendar_effect[dist_index], beta_bracelet_effect[dist_index] ]';
  }
 
  structural_treatment_effect = treatment_map_design_matrix * beta;
  structural_cluster_benefit_cost = structural_treatment_effect[cluster_assigned_dist_group_treatment];
  
  if (use_cluster_effects) {
    vector[num_clusters] cluster_effects;
    
    structural_beta_cluster = structural_beta_cluster_raw .* rep_matrix(structural_beta_cluster_sd, num_clusters);
    
    cluster_effects = rows_dot_product(cluster_treatment_design_matrix, structural_beta_cluster); 
    structural_cluster_benefit_cost += cluster_effects;
  }
  
  if (use_county_effects) {
    vector[num_clusters] county_effects;
    
    structural_beta_county = structural_beta_county_raw .* rep_matrix(structural_beta_county_sd, num_counties);
    
    county_effects = rows_dot_product(cluster_treatment_design_matrix, structural_beta_county[cluster_county_id]); 
    structural_cluster_benefit_cost += county_effects;
  }

  structural_cluster_obs_v = - structural_cluster_benefit_cost;
  
  structural_cluster_takeup_prob = Phi(- structural_cluster_obs_v);
}

model {
  beta_control[1] ~ normal(0, beta_control_sd);
  beta_control[2:] ~ normal(0, beta_far_effect_sd);
  
  beta_ink_effect ~ normal(0, beta_ink_effect_sd);
  beta_calendar_effect ~ normal(0, beta_calendar_effect_sd);
  beta_bracelet_effect ~ normal(0, beta_bracelet_effect_sd);
  
  if (use_cluster_effects) {
    to_vector(structural_beta_cluster_raw) ~ std_normal();
    structural_beta_cluster_sd ~ normal(0, structural_beta_cluster_sd_sd);
  }
  
  if (use_county_effects) {
    to_vector(structural_beta_county_raw) ~ std_normal();
    structural_beta_county_sd ~ normal(0, structural_beta_county_sd_sd);
  }
  
  if (!predict_prior) {
    // Take-up Likelihood 
    
    if (use_binomial) {
      cluster_takeup_count[included_clusters] ~ binomial(cluster_size[included_clusters], structural_cluster_takeup_prob[included_clusters]);
    } else {
      takeup[included_monitored_obs] ~ bernoulli(structural_cluster_takeup_prob[obs_cluster_id[included_monitored_obs]]);
    }
  }
}

generated quantities { 
  vector[num_clusters] cluster_cf_benefit_cost[num_dist_group_treatments]; 
  
  vector[generate_rep ? num_clusters : 0] cluster_rep_benefit_cost; 
  
  // Cross Validation
  vector[use_binomial || cluster_log_lik ? num_included_clusters : num_included_obs] log_lik = rep_vector(negative_infinity(), use_binomial || cluster_log_lik ? num_included_clusters : num_included_obs); 
  vector[use_binomial || cluster_log_lik ? num_excluded_clusters : num_excluded_obs] log_lik_heldout = rep_vector(negative_infinity(), use_binomial || cluster_log_lik ? num_excluded_clusters : num_excluded_obs); 
    
  {
    int treatment_cluster_pos = 1;
    int cluster_treatment_cf_pos = 1;
    
    for (treatment_index in 1:num_dist_group_treatments) {
      vector[num_clusters] treatment_dist_cost;
     
      treatment_dist_cost = rep_vector(0, num_clusters);
                                                                   
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
      
      real curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[cluster_index:cluster_index]);
      
      if (cluster_log_lik) {
        temp_log_lik[cluster_index] += curr_log_lik;
      } else {
        log_lik[obs_index_index] = curr_log_lik;
      }
    }
    
    for (obs_index_index in 1:num_excluded_monitored_obs) {
      int obs_index = excluded_monitored_obs[obs_index_index];
      int cluster_index = obs_cluster_id[obs_index];
      
      real  curr_log_lik = bernoulli_lpmf(takeup[obs_index] | structural_cluster_takeup_prob[cluster_index:cluster_index]);
      
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
  
  if (generate_rep) {
    for (cluster_index in 1:num_clusters) {
      vector[num_dist_group_treatments] rep_beta_cluster = rep_vector(0, num_dist_group_treatments);
      vector[num_dist_group_treatments] rep_beta_county = rep_vector(0, num_dist_group_treatments);
      
      int curr_assigned_treatment = cluster_treatment_map[cluster_assigned_dist_group_treatment[cluster_index], 1];
      
      rep_beta_cluster = use_cluster_effects ? to_vector(normal_rng(rep_array(0, num_dist_group_treatments), structural_beta_cluster_sd')) : rep_vector(0, num_dist_group_treatments);
      rep_beta_county = use_county_effects ? to_vector(normal_rng(rep_array(0, num_dist_group_treatments), structural_beta_county_sd')) : rep_vector(0, num_dist_group_treatments);
        
      cluster_rep_benefit_cost[cluster_index] = structural_treatment_effect[cluster_assigned_dist_group_treatment[cluster_index]]
        + treatment_map_design_matrix[cluster_assigned_dist_group_treatment[cluster_index]] * (rep_beta_cluster + rep_beta_county);
    }
  }
}
