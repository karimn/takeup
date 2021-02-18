for (grid_index in 1:num_grid_obs) {
  sim_grid_mu[grid_index] = rep_vector(0, 2);
}

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

if (num_excluded_obs > 0) {
  excluded_obs = which(obs_cluster_id, included_clusters, 0);
  excluded_monitored_obs = which(monitored_obs, excluded_obs, 1); 
  excluded_name_matched_obs = which(name_matched_obs, excluded_obs, 1);
}

{
  int max_close_index = num_dist_group_treatments / num_discrete_dist;
  
  treatment_map_design_matrix[, 1] = rep_vector(1, num_dist_group_treatments);
  
  for(treatment_index in 1:num_dist_group_treatments) {
    treatment_map_design_matrix[treatment_index, treatment_index] = 1;
    
    if (treatment_index > max_close_index) {
      treatment_map_design_matrix[treatment_index, max_close_index + 1] = 1;
      treatment_map_design_matrix[treatment_index, treatment_index - max_close_index] = 1;
    }
  }
  
  // print(treatment_map_design_matrix);
}

restricted_treatment_map_design_matrix = treatment_map_design_matrix;

if (use_private_incentive_restrictions) {
  // restricted_treatment_map_design_matrix[3,4] = 1; // Calendar effect is > then bracelet
  // restricted_treatment_map_design_matrix[7,8] = 1; // Calendar effect is > then bracelet
}

cluster_treatment_design_matrix = 
  treatment_map_design_matrix[cluster_assigned_dist_group_treatment]; 