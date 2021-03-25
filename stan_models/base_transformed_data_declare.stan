int num_dist_group_treatments = num_treatments * num_discrete_dist; 
int<lower = 1> cluster_size[num_clusters] = count(num_clusters, obs_cluster_id);