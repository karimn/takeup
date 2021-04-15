for (cluster_index in 1:num_clusters) {
  int curr_assigned_dist_group = cluster_treatment_map[cluster_assigned_dist_group_treatment[cluster_index], 2];
  
  for (dist_group_index in 1:num_discrete_dist) {
    vector[num_dist_group_mix] missing_group_dist_mix_lp;
    
    for (group_dist_mix_index in 1:num_dist_group_mix) {
      if (lognormal_dist_model) {
        missing_cluster_standard_dist[cluster_index, dist_group_index] = 
          group_dist_mix[dist_group_index, group_dist_mix_index] * 
          lognormal_rng(group_dist_mean[dist_group_index, group_dist_mix_index] + county_dist_effect[dist_group_index, cluster_county_id[cluster_index]] + cluster_dist_effect[dist_group_index, cluster_index], 
                        group_dist_sd[dist_group_index, group_dist_mix_index]);
      } else {
        missing_cluster_standard_dist[cluster_index, dist_group_index] = 
          group_dist_mix[dist_group_index, group_dist_mix_index] * 
          max({ 0.0, normal_rng(group_dist_mean[dist_group_index, group_dist_mix_index] + county_dist_effect[dist_group_index, cluster_county_id[cluster_index]] + cluster_dist_effect[dist_group_index, cluster_index],
          group_dist_sd[dist_group_index, group_dist_mix_index]) });
      }
    }
    
    all_cluster_standard_dist[cluster_index, dist_group_index] = 
      cluster_dist_treatment_id[cluster_index] == dist_group_index ? cluster_standard_dist[cluster_index] : missing_cluster_standard_dist[cluster_index, dist_group_index]; 
  }
}