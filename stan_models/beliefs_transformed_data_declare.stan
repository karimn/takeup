int<lower = 1, upper = num_clusters> beliefs_cluster_index[num_beliefs_obs] = obs_cluster_id[beliefs_obs_index];
int<lower = 1, upper = num_counties> beliefs_stratum_index[num_beliefs_obs] = obs_county_id[beliefs_obs_index];

int<lower = 1, upper = num_treatments * num_discrete_dist> beliefs_assigned_dist_group_treatment_id[num_beliefs_obs] = cluster_assigned_dist_group_treatment[beliefs_cluster_index];

matrix[num_beliefs_obs, num_dist_group_treatments] beliefs_treatment_design_matrix = beliefs_treatment_map_design_matrix[beliefs_assigned_dist_group_treatment_id];