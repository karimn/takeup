matrix<lower = 0, upper = 1>[num_beliefs_obs, num_dist_group_treatments] obs_prob_1ord; // = inv_logit(centered_obs_beta_1ord * beliefs_treatment_map_design_matrix');
vector<lower = -1, upper = 1>[num_dist_group_treatments] prob_1ord;
vector<lower = -1, upper = 1>[num_beliefs_ate_pairs] ate_1ord;


matrix[num_beliefs_obs, num_dist_group_treatments] obs_lin_pred_1ord; // = inv_logit(centered_obs_beta_1ord * beliefs_treatment_map_design_matrix');

matrix<lower = 0, upper = 1>[num_beliefs_obs, num_dist_group_treatments] obs_prob_2ord; // = inv_logit(centered_obs_beta_2ord * beliefs_treatment_map_design_matrix');
vector<lower = -1, upper = 1>[num_dist_group_treatments] prob_2ord;
vector<lower = -1, upper = 1>[num_beliefs_ate_pairs] ate_2ord;
  
for (treatment_index in 1:num_dist_group_treatments) {
  int curr_incentive_treatment_id = cluster_treatment_map[treatment_index, 1];
  int curr_assigned_dist_id = cluster_treatment_map[treatment_index, 2];
  obs_lin_pred_1ord[, treatment_index] = calculate_beliefs_latent_predictor(
      beliefs_treatment_map_design_matrix[curr_incentive_treatment_id:curr_incentive_treatment_id], 
      centered_obs_beta_1ord, 
      centered_obs_dist_beta_1ord, 
      all_cluster_standard_dist[beliefs_cluster_index, curr_assigned_dist_id]
    );

  obs_prob_1ord[, treatment_index] = inv_logit(
    obs_lin_pred_1ord[, treatment_index]
  );
      
  obs_prob_2ord[, treatment_index] = inv_logit(
    calculate_beliefs_latent_predictor(
      beliefs_treatment_map_design_matrix[curr_incentive_treatment_id:curr_incentive_treatment_id], 
      centered_obs_beta_2ord, 
      centered_obs_dist_beta_2ord, 
      all_cluster_standard_dist[beliefs_cluster_index, curr_assigned_dist_id]
    )
  );
  
  prob_1ord[treatment_index] = mean(obs_prob_1ord[, treatment_index]);
  prob_2ord[treatment_index] = mean(obs_prob_2ord[, treatment_index]);
}

ate_1ord = prob_1ord[beliefs_ate_pairs[, 1]] - prob_1ord[beliefs_ate_pairs[, 2]];
ate_2ord = prob_2ord[beliefs_ate_pairs[, 1]] - prob_2ord[beliefs_ate_pairs[, 2]];