ordered[num_dist_group_mix] group_dist_mean[num_discrete_dist];
matrix[num_discrete_dist, num_counties] county_dist_effect = rep_matrix(0, num_discrete_dist, num_counties);
matrix[num_discrete_dist, num_clusters] cluster_dist_effect = rep_matrix(0, num_discrete_dist, num_clusters);
