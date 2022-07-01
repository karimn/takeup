int<lower = 0, upper = 1> use_age_groups = num_age_groups > 1;

int num_dist_group_treatments = num_treatments * num_discrete_dist; 
// array[num_clusters] int<lower = 1> cluster_size = count(num_clusters, obs_cluster_id);
array[num_clusters, num_age_groups] int<lower = 0> cluster_size = count(num_clusters, num_age_groups, obs_cluster_id, obs_age_group);

array[num_clusters] int<lower = 1, upper = num_treatments> cluster_incentive_treatment_id = cluster_treatment_map[cluster_assigned_dist_group_treatment, 1];
array[num_clusters] int<lower = 1, upper = num_dist_group_treatments> cluster_dist_treatment_id = cluster_treatment_map[cluster_assigned_dist_group_treatment, 2];
      