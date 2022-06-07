int<lower = 0, upper = 1> fit_wtp_model_to_data;
int<lower = 0, upper = 1> use_strata_levels;

int<lower = 0> num_wtp_obs;
int<lower = 0> num_strata;
array[num_strata] int<lower = 1, upper = num_wtp_obs> wtp_strata_sizes; 

// Calendar = 1
// Bracelet = -1
vector<lower = -1, upper = 1>[num_wtp_obs] gift_choice;

// Switch = 1
// Keep = -1
vector<lower = -1, upper = 1>[num_wtp_obs] wtp_response;

vector<lower = 0>[num_wtp_obs] wtp_offer;

int<lower = 0> num_preference_value_diff;
vector[num_preference_value_diff] preference_value_diff;