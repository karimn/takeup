array[num_discrete_dist] simplex[num_dist_group_mix] group_dist_mix;

array[num_discrete_dist] ordered[num_dist_group_mix] group_dist_mean_raw;
array[num_discrete_dist] vector<lower = 0>[num_dist_group_mix] group_dist_sd;

matrix[use_dist_county_effects ? num_discrete_dist : 0, num_counties] county_dist_effect_raw;
vector<lower = 0>[use_dist_county_effects ? num_discrete_dist : 0] county_dist_effect_sd;

matrix[use_dist_cluster_effects ? num_discrete_dist : 0, num_clusters] cluster_dist_effect_raw;
vector<lower = 0>[use_dist_county_effects ? num_discrete_dist : 0] cluster_dist_effect_sd;
  