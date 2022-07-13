matrix[num_clusters, num_roc_distances] cluster_roc_left = rep_matrix(0, num_clusters, num_roc_distances);
matrix[num_clusters, num_roc_distances] cluster_roc_right = rep_matrix(0, num_clusters, num_roc_distances);
matrix[num_clusters, num_roc_distances] cluster_roc_diff = rep_matrix(0, num_clusters, num_roc_distances);
matrix[num_clusters, num_roc_distances] cluster_w_cutoff_left = rep_matrix(0, num_clusters, num_roc_distances);
matrix[num_clusters, num_roc_distances] cluster_w_cutoff_right = rep_matrix(0, num_clusters, num_roc_distances);
matrix[num_clusters, num_roc_distances] cluster_takeup_prop_left = rep_matrix(0, num_clusters, num_roc_distances);
matrix[num_clusters, num_roc_distances] cluster_takeup_prop_right = rep_matrix(0, num_clusters, num_roc_distances);
matrix[num_clusters, num_roc_distances] cluster_social_multiplier_left = rep_matrix(0, num_clusters, num_roc_distances);
matrix[num_clusters, num_roc_distances] cluster_social_multiplier_right = rep_matrix(0, num_clusters, num_roc_distances);
matrix[num_clusters, num_roc_distances] cluster_rep_return_left = rep_matrix(0, num_clusters, num_roc_distances);
matrix[num_clusters, num_roc_distances] cluster_rep_return_right = rep_matrix(0, num_clusters, num_roc_distances);
 
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
      
    matrix[num_clusters, 2] curr_cluster_mu_rep_right = calculate_mu_rep_deriv(
      roc_compare_treatment_id_right, roc_cluster_dist,
      base_mu_rep, 1, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord);
     
    matrix[num_clusters, 9] roc_results = map_calculate_roc_diff(
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