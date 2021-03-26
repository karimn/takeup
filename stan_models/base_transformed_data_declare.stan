int num_dist_group_treatments = num_treatments * num_discrete_dist; 
int<lower = 1> cluster_size[num_clusters] = count(num_clusters, obs_cluster_id);

int<lower = 1, upper = num_treatments> cluster_incentive_treatment_id[num_clusters] = cluster_treatment_map[cluster_assigned_dist_group_treatment, 1];
      