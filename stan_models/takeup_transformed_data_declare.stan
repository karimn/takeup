int num_dist_group_treatments = num_treatments * num_discrete_dist; 
  
matrix<lower = 0, upper = 1>[num_dist_group_treatments, num_dist_group_treatments] treatment_map_design_matrix = rep_matrix(0, num_dist_group_treatments, num_dist_group_treatments);
matrix<lower = 0, upper = 1>[num_dist_group_treatments, num_dist_group_treatments] restricted_treatment_map_design_matrix;
matrix<lower = 0, upper = 1>[num_clusters, num_dist_group_treatments] cluster_treatment_design_matrix = rep_matrix(0, num_clusters, num_dist_group_treatments);

int<lower = 1> cluster_size[num_clusters] = count(num_clusters, obs_cluster_id);
int<lower = 0> cluster_takeup_count[num_clusters] = count_by_group_test(takeup, obs_cluster_id, { 1 }, 1);

int<lower = 1, upper = num_dist_group_treatments> assigned_treatment[num_obs] = cluster_assigned_dist_group_treatment[obs_cluster_id]; 
  
int<lower = 0, upper = num_obs> num_name_matched = sum(is_name_matched);
int<lower = 0, upper = num_obs> num_monitored = num_obs - sum(is_name_matched);

int<lower = 1, upper = num_obs> monitored_obs[num_monitored] = which(is_name_matched, { 0 }, 1);
int<lower = 1, upper = num_obs> name_matched_obs[num_name_matched] = which(is_name_matched, { 1 }, 1);

int<lower = 0, upper = num_clusters> num_included_clusters = num_clusters - num_excluded_clusters;
int<lower = 1, upper = num_clusters> included_clusters[num_included_clusters] = num_excluded_clusters > 0 ? which(seq(1, num_clusters, 1), excluded_clusters, 0) : seq(1, num_clusters, 1);
int<lower = 0, upper = num_obs> num_included_obs = sum(cluster_size[included_clusters]);
int<lower = 0, upper = num_obs> num_excluded_obs = num_obs - num_included_obs;
int<lower = 1, upper = num_obs> included_obs[num_included_obs] = num_included_obs != num_obs ? which(obs_cluster_id, included_clusters, 1) : seq(1, num_obs, 1);
int<lower = 1, upper = num_obs> excluded_obs[num_excluded_obs]; 

int<lower = 0, upper = num_obs> num_included_monitored_obs = num_equals(monitored_obs, included_obs);
int<lower = 0, upper = num_obs> num_excluded_monitored_obs = num_monitored - num_included_monitored_obs;
int<lower = 0, upper = num_obs> num_included_name_matched_obs = num_equals(name_matched_obs, included_obs);
int<lower = 0, upper = num_obs> num_excluded_name_matched_obs = num_name_matched - num_included_name_matched_obs;
int<lower = 1, upper = num_obs> included_monitored_obs[num_included_monitored_obs] = which(monitored_obs, included_obs, 1);
int<lower = 1, upper = num_obs> excluded_monitored_obs[num_excluded_monitored_obs]; 
int<lower = 1, upper = num_obs> included_name_matched_obs[num_included_name_matched_obs] = which(name_matched_obs, included_obs, 1);
int<lower = 1, upper = num_obs> excluded_name_matched_obs[num_excluded_name_matched_obs];

int<lower = 0, upper = 1> use_dist_salience = in_array(use_cost_model, 
                                                       { COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE, COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE, COST_MODEL_TYPE_SEMIPARAM_SALIENCE });
int<lower = 0, upper = 1> use_semiparam = in_array(use_cost_model, { COST_MODEL_TYPE_SEMIPARAM, COST_MODEL_TYPE_SEMIPARAM_SALIENCE }); 

int num_treatments_param_kappa = use_cost_model == COST_MODEL_TYPE_PARAM_KAPPA ? num_treatments : 0;
int num_treatments_param = 0;
int num_treatments_param_quadratic = 0; 
int num_treatments_semiparam = 0;

vector[2] sim_grid_mu[num_grid_obs];

int<lower = 1, upper = num_clusters * num_dist_group_treatments> long_cluster_by_treatment_index[num_clusters] =
  array_add(array_product(array_subtract(cluster_assigned_dist_group_treatment, { 1 }), { num_clusters }), seq(1, num_clusters, 1));
  
int<lower = 1, upper = num_clusters * num_treatments> long_cluster_by_incentive_treatment_index[num_clusters] =
  array_add(array_product(array_subtract(cluster_treatment_map[cluster_assigned_dist_group_treatment, 1], { 1 }), { num_clusters }), seq(1, num_clusters, 1));

int<lower = 1> num_dist_group_mix = 2;
