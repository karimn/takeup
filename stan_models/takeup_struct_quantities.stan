array[num_roc_distances] matrix[num_clusters, num_treatments] cluster_roc; 
array[num_roc_distances] matrix[num_clusters, num_treatments - 1] cluster_roc_diff; 
array[num_roc_distances] matrix[num_clusters, num_treatments] cluster_w_cutoff;
array[num_roc_distances] matrix[num_clusters, num_treatments] cluster_takeup_prop;
array[num_roc_distances] matrix[num_clusters, num_treatments] cluster_social_multiplier;
array[num_roc_distances] matrix[num_clusters, num_treatments] cluster_rep_return;
array[num_roc_distances] matrix[num_clusters, num_treatments] cluster_rep_return_dist;
 
if (multithreaded) { 
  for (roc_dist_index in 1:num_roc_distances) {
    vector[num_clusters] roc_cluster_dist = rep_vector(roc_distances[roc_dist_index], num_clusters);
    
    matrix[num_clusters, num_treatments] curr_net_benefit = 
      structural_cluster_benefit[, :num_treatments] - param_dist_cost(roc_distances[roc_dist_index], cluster_linear_dist_cost[, :num_treatments], cluster_quadratic_dist_cost[, :num_treatments]);
    
    int diff_index = 1;
    // 
    // cluster_w_cutoff[roc_dist_index] = rep_matrix(0, num_clusters, num_treatments);
    // cluster_takeup_prop[roc_dist_index] = rep_matrix(0, num_clusters, num_treatments);
    // cluster_social_multiplier[roc_dist_index] = rep_matrix(0, num_clusters, num_treatments);
    // cluster_rep_return[roc_dist_index] = rep_matrix(0, num_clusters, num_treatments);
    // cluster_roc[roc_dist_index] = rep_matrix(0, num_clusters, num_treatments);
    // cluster_roc_diff[roc_dist_index] = rep_matrix(0, num_clusters, num_treatments - 1);
    
    matrix[num_clusters, 2] curr_cluster_mu_rep_control = calculate_mu_rep_deriv(
      roc_compare_treatment_id_right, roc_cluster_dist,
      base_mu_rep, 1, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord);
    
    for (treatment_index in 1:num_treatments) {
      matrix[num_clusters, 2] curr_cluster_mu_rep; 
      
      if (treatment_index > 1) {
        curr_cluster_mu_rep = calculate_mu_rep_deriv(
          treatment_index, roc_cluster_dist,
          base_mu_rep, 1, beliefs_treatment_map_design_matrix, centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord);
      } else {
        curr_cluster_mu_rep = curr_cluster_mu_rep_control;
      }
        
      matrix[num_clusters, 5] roc_results = map_calculate_roc(
        curr_net_benefit[, treatment_index],
        curr_net_benefit[, roc_compare_treatment_id_right],
        
        rep_vector(total_error_sd[treatment_index], num_clusters),
        rep_vector(total_error_sd[roc_compare_treatment_id_right], num_clusters),
        rep_vector(u_sd[treatment_index], num_clusters),
        rep_vector(u_sd[roc_compare_treatment_id_right], num_clusters),
        
        cluster_linear_dist_cost[, treatment_index],
        
        curr_cluster_mu_rep[, 1], 
        curr_cluster_mu_rep_control[, 1], 
        curr_cluster_mu_rep[, 2],
        
        use_u_in_delta,
        alg_sol_rel_tol,
        alg_sol_f_tol,
        alg_sol_max_steps
      );
      
      cluster_w_cutoff[roc_dist_index, , treatment_index] = roc_results[, 1];
      cluster_takeup_prop[roc_dist_index, , treatment_index] = 1 - Phi_approx(roc_results[, 1] / total_error_sd[treatment_index]);
      cluster_social_multiplier[roc_dist_index, , treatment_index] = - roc_results[, 3];
      cluster_rep_return[roc_dist_index, , treatment_index] = roc_results[, 3] .* curr_cluster_mu_rep[, 1];
      cluster_rep_return_dist[roc_dist_index, , treatment_index] = cluster_rep_return[roc_dist_index, , treatment_index] / dist_beta_v[1]; 
      cluster_roc[roc_dist_index, , treatment_index] = roc_results[, 5];
      
      if (treatment_index > roc_compare_treatment_id_right) {
        cluster_roc_diff[roc_dist_index, , diff_index] = cluster_roc[roc_dist_index, , treatment_index] - cluster_roc[roc_dist_index, , roc_compare_treatment_id_right]; 
        diff_index += 1;
      }
    }
  }
}