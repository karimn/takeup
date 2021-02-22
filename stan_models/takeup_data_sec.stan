int MIN_COST_MODEL_TYPE_VALUE;
int MAX_COST_MODEL_TYPE_VALUE;
int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_KAPPA;
int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_LINEAR;
int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_QUADRATIC;
int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_SEMIPARAM;
int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_LINEAR_SALIENCE;
int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_PARAM_QUADRATIC_SALIENCE;
int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_SEMIPARAM_SALIENCE;
int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> COST_MODEL_TYPE_DISCRETE;

int<lower = 1> num_v_mix;

// Multilevel Configuration 
int<lower = 0, upper = 1> use_cluster_effects;
int<lower = 0, upper = 1> use_county_effects;
int<lower = 0, upper = 1> use_mu_cluster_effects;
int<lower = 0, upper = 1> use_mu_county_effects;
int<lower = 0, upper = 1> use_param_dist_cluster_effects; // These are used for parameteric (linear, quadratic) distance cost models only
int<lower = 0, upper = 1> use_param_dist_county_effects;

int<lower = 0, upper = 1> use_binomial;
int<lower = 0, upper = 1> use_cost_k_restrictions;
int<lower = 0, upper = 1> use_private_incentive_restrictions;
int<lower = 0, upper = 1> use_salience_effect;
int<lower = 0, upper = 1> use_single_cost_model;
int<lower = 0, upper = 1> use_name_matched_obs;
int<lower = 0, upper = 1> use_shifting_v_dist;
int<lower = MIN_COST_MODEL_TYPE_VALUE, upper = MAX_COST_MODEL_TYPE_VALUE> use_cost_model;
int<lower = 0, upper = 1> suppress_reputation;
int<lower = 0, upper = 1> use_u_in_delta;
int<lower = 0, upper = 1> multithreaded;
int<lower = 0, upper = 1> generate_rep;
int<lower = 0, upper = 1> generate_sim;
int<lower = 0, upper = 1> fit_model_to_data;
int<lower = 0, upper = 1> cross_validate;

real<lower = 0> alg_sol_f_tol;
real<lower = 0> alg_sol_rel_tol;
int<lower = 0> alg_sol_max_steps;

int<lower = 1> num_obs; // Actual observations
int<lower = 1> num_clusters;
int<lower = 1> num_counties;
int<lower = 1> num_treatments;
int<lower = 1> num_discrete_dist;

int<lower = 1, upper = num_treatments> CALENDAR_TREATMENT_INDEX;
int<lower = 1, upper = num_treatments> BRACELET_TREATMENT_INDEX;

int<lower = 1, upper = num_clusters> obs_cluster_id[num_obs];
int<lower = 1, upper = num_counties> obs_county_id[num_obs];
int<lower = 1, upper = num_counties> cluster_county_id[num_clusters];

int<lower = 0, upper = 1> takeup[num_obs]; // Observed outcome variable
int<lower = 0, upper = 1> is_name_matched[num_obs];
vector[num_clusters] cluster_standard_dist; // Standardized distance to treatment

int<lower = 1> cluster_treatment_map[num_treatments * num_discrete_dist, 2];
int<lower = 1, upper = num_treatments * num_discrete_dist> cluster_assigned_dist_group_treatment[num_clusters]; 

// Reputation

// Semiparametric Cost Model (Splines)

int<lower = 3> num_knots_v;
matrix[num_clusters, num_knots_v] Z_splines_v; 

real<lower = 0> u_splines_v_sigma_sd;

// K-fold CV 

int<lower = 0, upper =1> cluster_log_lik; 
int<lower = 0> num_excluded_clusters;
int<lower = 1, upper = num_clusters> excluded_clusters[num_excluded_clusters];

// Simulation

int<lower = 1> num_grid_obs; // Simulation observations
int<lower = 1> num_small_grid_obs; // Simulation observations
vector[num_grid_obs] grid_dist; // Simulation distances
vector[num_small_grid_obs] small_grid_dist; // Simulation distances
matrix[num_grid_obs, num_knots_v] Z_grid_v;

// int<lower = 0> num_sim_sm_v;
// vector[num_sim_sm_v] sim_sm_v;

// Prior hyperparameters

real<lower = 0> mu_rep_sd;
real<lower = 0> beta_control_sd;
real<lower = 0> beta_far_effect_sd;
real<lower = 0> beta_ink_effect_sd;
real<lower = 0> beta_calendar_effect_sd;
real<lower = 0> beta_bracelet_effect_sd;

real<lower = 0> dist_beta_v_sd;

real<lower = 0> structural_beta_county_sd_sd;
real<lower = 0> structural_beta_cluster_sd_sd;