simplex[num_dist_group_mix] group_dist_mix[num_discrete_dist];

ordered[num_dist_group_mix] group_dist_mean[num_discrete_dist];
vector<lower = 0>[num_dist_group_mix] group_dist_sd[num_discrete_dist];

matrix<lower = 0>[num_clusters, num_discrete_dist - 1] missing_cluster_standard_dist; 
  