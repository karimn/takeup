int num_dist_group_treatments = num_treatments * num_discrete_dist; 
array[num_clusters] int<lower = 1> cluster_size = count(num_clusters, obs_cluster_id);

array[num_clusters] int<lower = 1, upper = num_treatments> cluster_incentive_treatment_id = cluster_treatment_map[cluster_assigned_dist_group_treatment, 1];
array[num_clusters] int<lower = 1, upper = num_dist_group_treatments> cluster_dist_treatment_id = cluster_treatment_map[cluster_assigned_dist_group_treatment, 2];
      