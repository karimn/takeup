// Multilevel Configuration 
int<lower = 0, upper = 1> use_cluster_effects;
int<lower = 0, upper = 1> use_county_effects;

int<lower = 0, upper = 1> use_binomial;
int<lower = 0, upper = 1> use_private_incentive_restrictions;
int<lower = 0, upper = 1> use_name_matched_obs;
int<lower = 0, upper = 1> multithreaded;
int<lower = 0, upper = 1> generate_rep;
int<lower = 0, upper = 1> generate_sim;
int<lower = 0, upper = 1> fit_model_to_data;
int<lower = 0, upper = 1> cross_validate;

real<lower = 0> alg_sol_f_tol;
real<lower = 0> alg_sol_rel_tol;
int<lower = 0> alg_sol_max_steps;

int<lower = 1, upper = num_treatments> CALENDAR_TREATMENT_INDEX;
int<lower = 1, upper = num_treatments> BRACELET_TREATMENT_INDEX;

array[num_obs] int<lower = 0, upper = 1> takeup; // Observed outcome variable
array[num_obs] int<lower = 0, upper = 1> is_name_matched;

// K-fold CV 

int<lower = 0, upper =1> cluster_log_lik; 
int<lower = 0> num_excluded_clusters;
array[num_excluded_clusters] int<lower = 1, upper = num_clusters> excluded_clusters;

// Prior hyperparameters

real<lower = 0> beta_control_sd;
real<lower = 0> beta_ink_effect_sd;
real<lower = 0> beta_calendar_effect_sd;
real<lower = 0> beta_bracelet_effect_sd;

real<lower = 0> dist_beta_v_sd;

real<lower = 0> structural_beta_county_sd_sd;
real<lower = 0> structural_beta_cluster_sd_sd;