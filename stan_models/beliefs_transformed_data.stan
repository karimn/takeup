array[num_beliefs_obs] int<lower = 1, upper = num_clusters> beliefs_cluster_index = obs_cluster_id[beliefs_obs_index];
array[num_beliefs_obs] int<lower = 1, upper = num_counties> beliefs_stratum_index = obs_county_id[beliefs_obs_index];

// int<lower = 1, upper = num_treatments * num_discrete_dist> beliefs_assigned_dist_group_treatment_id[num_beliefs_obs] = cluster_assigned_dist_group_treatment[beliefs_cluster_index];
array[num_beliefs_obs] int<lower = 1, upper = num_treatments> obs_beliefs_treatment_id = cluster_treatment_map[cluster_assigned_dist_group_treatment[beliefs_cluster_index], 1];
array[num_beliefs_obs] int<lower = 1, upper = num_discrete_dist> obs_beliefs_assigned_dist_id = cluster_treatment_map[cluster_assigned_dist_group_treatment[beliefs_cluster_index], 2];

// matrix[num_beliefs_obs, num_dist_group_treatments] beliefs_treatment_design_matrix = beliefs_treatment_map_design_matrix[beliefs_assigned_dist_group_treatment_id];
matrix[num_beliefs_obs, num_treatments] beliefs_treatment_design_matrix = beliefs_treatment_map_design_matrix[obs_beliefs_treatment_id];