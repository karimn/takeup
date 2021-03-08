for (grid_index in 1:num_grid_obs) {
  sim_grid_mu[grid_index] = rep_vector(0, 2);
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