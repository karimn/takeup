matrix<lower = 0, upper = 1>[num_dist_group_treatments, num_dist_group_treatments] treatment_map_design_matrix = rep_matrix(0, num_dist_group_treatments, num_dist_group_treatments);
matrix<lower = 0, upper = 1>[num_dist_group_treatments, num_dist_group_treatments] restricted_treatment_map_design_matrix;
matrix<lower = 0, upper = 1>[num_clusters, num_dist_group_treatments] cluster_treatment_design_matrix = rep_matrix(0, num_clusters, num_dist_group_treatments);

// array[num_clusters] int<lower = 0> cluster_takeup_count = count_by_group_test(takeup, obs_cluster_id, { 1 }, 1);
array[num_clusters, num_age_groups] int<lower = 0> cluster_takeup_count = count_by_group_test(takeup, obs_cluster_id, obs_age_group, { 1 }, 1);

// array[num_clusters, num_age_groups] int<lower = 0> cluster_size = count(num_clusters, num_age_groups, obs_cluster_id, obs_age_group);
matrix<lower = 0, upper = 1>[num_clusters, num_age_groups] cluster_age_group_prop;

array[num_obs] int<lower = 1, upper = num_dist_group_treatments> assigned_treatment = cluster_assigned_dist_group_treatment[obs_cluster_id]; 
  
int<lower = 0, upper = num_obs> num_name_matched = sum(is_name_matched);
int<lower = 0, upper = num_obs> num_monitored = num_obs - sum(is_name_matched);

array[num_monitored] int<lower = 1, upper = num_obs> monitored_obs = which(is_name_matched, { 0 }, 1);
array[num_name_matched] int<lower = 1, upper = num_obs> name_matched_obs = which(is_name_matched, { 1 }, 1);

int<lower = 0, upper = num_clusters> num_included_clusters = num_clusters - num_excluded_clusters;
array[num_included_clusters] int<lower = 1, upper = num_clusters> included_clusters = num_excluded_clusters > 0 ? which(seq(1, num_clusters, 1), excluded_clusters, 0) : seq(1, num_clusters, 1);
int<lower = 0, upper = num_obs> num_included_obs = sum(to_array_1d(cluster_size[included_clusters]));
int<lower = 0, upper = num_obs> num_excluded_obs = num_obs - num_included_obs;
array[num_included_obs] int<lower = 1, upper = num_obs> included_obs = num_included_obs != num_obs ? which(obs_cluster_id, included_clusters, 1) : seq(1, num_obs, 1);
array[num_excluded_obs] int<lower = 1, upper = num_obs> excluded_obs; 

int<lower = 0, upper = num_obs> num_included_monitored_obs = num_equals(monitored_obs, included_obs);
int<lower = 0, upper = num_obs> num_excluded_monitored_obs = num_monitored - num_included_monitored_obs;
int<lower = 0, upper = num_obs> num_included_name_matched_obs = num_equals(name_matched_obs, included_obs);
int<lower = 0, upper = num_obs> num_excluded_name_matched_obs = num_name_matched - num_included_name_matched_obs;
array[num_included_monitored_obs] int<lower = 1, upper = num_obs> included_monitored_obs = which(monitored_obs, included_obs, 1);
array[num_excluded_monitored_obs] int<lower = 1, upper = num_obs> excluded_monitored_obs; 
array[num_included_name_matched_obs] int<lower = 1, upper = num_obs> included_name_matched_obs = which(name_matched_obs, included_obs, 1);
array[num_excluded_name_matched_obs] int<lower = 1, upper = num_obs> excluded_name_matched_obs;

array[num_grid_obs] vector[2] sim_grid_mu;

array[num_clusters] int<lower = 1, upper = num_clusters * num_dist_group_treatments> long_cluster_by_treatment_index =
  array_add(array_product(array_subtract(cluster_assigned_dist_group_treatment, { 1 }), { num_clusters }), seq(1, num_clusters, 1));
  
array[num_clusters] int<lower = 1, upper = num_clusters * num_treatments> long_cluster_by_incentive_treatment_index =
  array_add(array_product(array_subtract(cluster_treatment_map[cluster_assigned_dist_group_treatment, 1], { 1 }), { num_clusters }), seq(1, num_clusters, 1));
  
cluster_age_group_prop = to_matrix(cluster_size);
cluster_age_group_prop = cluster_age_group_prop ./ rep_matrix(row_sum(cluster_age_group_prop), num_age_groups);

for (grid_index in 1:num_grid_obs) {
  sim_grid_mu[grid_index] = rep_vector(0, 2);
}

if (num_excluded_obs > 0) {
  excluded_obs = which(obs_cluster_id, included_clusters, 0);
  excluded_monitored_obs = which(monitored_obs, excluded_obs, 1); 
  excluded_name_matched_obs = which(name_matched_obs, excluded_obs, 1);
}

{
  int max_close_index = num_dist_group_treatments %/% num_discrete_dist;
  
  treatment_map_design_matrix[, 1] = rep_vector(1, num_dist_group_treatments);
  
  for(treatment_index in 1:num_dist_group_treatments) {
    treatment_map_design_matrix[treatment_index, treatment_index] = 1;
    
    if (treatment_index > max_close_index) {
      treatment_map_design_matrix[treatment_index, max_close_index + 1] = 1;
      treatment_map_design_matrix[treatment_index, treatment_index - max_close_index] = 1;
    }
  }
}

restricted_treatment_map_design_matrix = treatment_map_design_matrix;

if (use_private_incentive_restrictions) {
  // restricted_treatment_map_design_matrix[3,4] = 1; // Calendar effect is > then bracelet
  // restricted_treatment_map_design_matrix[7,8] = 1; // Calendar effect is > then bracelet
}

cluster_treatment_design_matrix = 
  treatment_map_design_matrix[cluster_assigned_dist_group_treatment]; 