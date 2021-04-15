for (dist_group_index in 1:num_discrete_dist) { 
  if (use_dist_county_effects) { 
    county_dist_effect[dist_group_index] = county_dist_effect_raw[dist_group_index] * county_dist_effect_sd[dist_group_index];  
  }
  
  if (use_dist_cluster_effects) { 
    cluster_dist_effect[dist_group_index] = cluster_dist_effect_raw[dist_group_index] * cluster_dist_effect_sd[dist_group_index];  
  }
  
  group_dist_mean[dist_group_index] = hyper_dist_mean_mean + hyper_dist_mean_sd * group_dist_mean_raw[dist_group_index]; 
}