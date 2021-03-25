functions {
#include takeup_functions.stan
}

data { 
#include base_data_sec.stan
#include beliefs_data_sec.stan
}

transformed data {
#include base_transformed_data_declare.stan
#include beliefs_transformed_data_declare.stan
}

parameters {
#include beliefs_parameters_sec.stan
}

transformed parameters {
#include beliefs_transformed_parameters_declare.stan
#include beliefs_transformed_parameters_define.stan
}

model {
#include beliefs_model_sec.stan
}

generated quantities {
#include beliefs_generated_quantities_declare.stan
#include beliefs_generated_quantities_define.stan
}
