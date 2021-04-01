for (treatment_index in 1:num_dist_group_treatments) {
  int curr_incentive_treatment_id = cluster_treatment_map[treatment_index, 1];
  int curr_assigned_dist_id = cluster_treatment_map[treatment_index, 2];
  
  obs_prob_1ord[, treatment_index] = inv_logit(
    calculate_beliefs_latent_predictor(
      beliefs_treatment_map_design_matrix[curr_incentive_treatment_id:curr_incentive_treatment_id], 
      centered_obs_beta_1ord, 
      centered_obs_dist_beta_1ord, 
      missing_cluster_standard_dist[beliefs_cluster_index, curr_assigned_dist_id]
    )
  );
      
    // centered_obs_beta_1ord * beliefs_treatment_map_design_matrix[curr_incentive_treatment_id]' + 
    // (centered_obs_dist_beta_1ord * beliefs_treatment_map_design_matrix[curr_incentive_treatment_id]') .* missing_cluster_standard_dist[beliefs_cluster_index, curr_assigned_dist_id]);
    
  // obs_prob_2ord[, treatment_index] = inv_logit(
  //   centered_obs_beta_2ord * beliefs_treatment_map_design_matrix[curr_incentive_treatment_id]' + 
  //   (centered_obs_dist_beta_2ord * beliefs_treatment_map_design_matrix[curr_incentive_treatment_id]') .* missing_cluster_standard_dist[beliefs_cluster_index, curr_assigned_dist_id]);
    
  obs_prob_2ord[, treatment_index] = inv_logit(
    calculate_beliefs_latent_predictor(
      beliefs_treatment_map_design_matrix[curr_incentive_treatment_id:curr_incentive_treatment_id], 
      centered_obs_beta_2ord, 
      centered_obs_dist_beta_2ord, 
      missing_cluster_standard_dist[beliefs_cluster_index, curr_assigned_dist_id]
    )
  );
  
  prob_1ord[treatment_index] = mean(obs_prob_1ord[, treatment_index]);
  prob_2ord[treatment_index] = mean(obs_prob_2ord[, treatment_index]);
}

ate_1ord = prob_1ord[beliefs_ate_pairs[, 1]] - prob_1ord[beliefs_ate_pairs[, 2]];
ate_2ord = prob_2ord[beliefs_ate_pairs[, 1]] - prob_2ord[beliefs_ate_pairs[, 2]];
  