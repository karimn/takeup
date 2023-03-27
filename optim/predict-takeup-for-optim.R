#!/usr/bin/Rscript
script_options = docopt::docopt(
    stringr::str_glue("Usage:
    predict-takeup-for-optim.R <fit-version> <private-benefit-z> <visibility-z> [options] [ --from-csv | --to-csv ]

    Takes posterior fits and calculates estimated takeup for a given distance.


    Options:
        --from-csv  Load parameter draws from a csv file
        --to-csv   Load stanfit object and write draws to csv
        --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
        --output-path=<path>  Path to find results [default: {file.path('optim', 'data')}]
        --output-name=<output-name>  Prepended to output file
        --pred-distance  Estimate takeup on a grid from 0 to 5km
        --dist-cutoff=<dist-cutoff>  Don't estimate takeup past this cutoff (set takeup to 0) [default: 2500]
        --rep-cutoff=<rep-cutoff>  Don't let reputational returns increase past this point. [default: Inf]
        --num-post-draws=<num-post-draws>  Number of posterior draws to use [default: 200]
        --num-cores=<num-cores>  Number of cores to use [default: 8]
        --type-lb=<type-lb>  Lower bound of type distribution [default: -Inf]
        --type-ub=<type-ub>  Upper bound of type distribution [default: Inf]
        --model=<model>  Which model to use [default: STRUCTURAL_LINEAR_U_SHOCKS]
        --fit-rf  Estimate a simple brms reduced form model with distance entering continuously 
        --run-estimation  Whether th run estimation stuff
        --num-extra-pots=<num-extra-pots>  How many extra PoTs to sample per village above experiment [default: 2]
        --data-input-path=<data-input-path>  Where village and PoT data is stored [default: {file.path('optim', 'data')}]
        --data-input-name=<data-input-name>  Filename of village and PoT data [default: full-experiment-4-extra-pots.rds]
        --suppress-reputation  Suppress reputational returns
        --single-chain  Only use first chain for draws (useful for debugging) 
    "),
    args = if (interactive()) "
                            85
                            control
                            control
                            --output-name=cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
                            --to-csv
                            --num-post-draws=200
                            --rep-cutoff=Inf
                            --dist-cutoff=3500
                            --type-lb=-Inf
                            --type-ub=Inf
                            --num-cores=12
                            --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
                            --num-extra-pots=4
                            --data-input-name=full-experiment.rds
                            --single-chain
                            --run-estimation
                              " 
           else commandArgs(trailingOnly = TRUE)
)


# Loading functions
run_estimation = script_options$run_estimation
set.seed(19484)

library(posterior)
library(tidyverse)
library(tidybayes)
library(broom)
library(rstan)
library(sf)
library(nleqslv)
library(cmdstanr)
library(econometr)
library(furrr)

source(file.path("optim", "optim-functions.R"))
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))


source(file.path("multilvlr", "multilvlr_util.R"))
source(file.path("stan_models", "stan_owen_t.R"))

stan_owen_t = exposed_stan_func$stan_owen_t

script_options = script_options %>%
    modify_at(c("dist_cutoff", 
                "rep_cutoff",
                "num_post_draws", 
                "num_cores",
                "num_extra_pots",
                "type_lb",
                "type_ub"), as.numeric)
fit_version = script_options$fit_version

script_options$bounds = c(script_options$type_lb, script_options$type_ub)

append_output = if (!is.null(script_options$output_name)) {
    paste0("-", script_options$output_name)
} else {
    ""
}

mu_rep_type = switch(
    script_options$model,
    STRUCTURAL_LINEAR_U_SHOCKS = 0,
    STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP = 1,
    STRUCTURAL_LINEAR_U_SHOCKS_LINEAR_MU_REP = 2,
    STRUCTURAL_LINEAR_U_SHOCKS_NO_REP = 3,
    STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP = 4
)


models_we_want = c(
    script_options$model
)


# sd_of_dist = sd(rf_analysis_data$cluster.dist.to.pot)


if (script_options$to_csv) {
    if (script_options$single_chain) {
        struct_model_files = fs::dir_ls(
            script_options$input_path, 
            regex = str_glue("dist_fit{fit_version}_{script_options$model}-1\\.csv"))
    } else {
        struct_model_files = fs::dir_ls(
            script_options$input_path, 
            regex = str_glue("dist_fit{fit_version}_{script_options$model}-.*csv"))
    }
    
    struct_model_fit = as_cmdstan_fit(struct_model_files)
    # N.B. this is very hardcoded for vanilla structural model
    struct_param_draws = struct_model_fit %>%
        gather_draws(
            beta[k],
            dist_beta_v[j],
            centered_cluster_beta_1ord[j, k],
            centered_cluster_dist_beta_1ord[j, k],
            base_mu_rep, 
            total_error_sd[k],
            u_sd[k]
        )
    # model_fit is huge so get rid of it asap
    rm(struct_model_fit) 
    gc()
    # N.B. we use input path for `to_csv` output and reserve output_path for 
    # output later in script
    struct_param_draws %>%
        write_csv(file.path(script_options$input_path, str_interp("param_posterior_draws_dist_fit${fit_version}_${script_options$model}.csv") ))
    ## Now reduced form  - we actually fit rf model continuous in distance here
    # since main analysis RF just uses close/far
    if (script_options$fit_rf) {

        # Fit Loading
        load(file.path("temp-data", str_interp("processed_dist_fit${fit_version}_lite.RData")))


        stan_data = (dist_fit_data %>%
            slice(1) %>% # just take first model as only using for analysis_data
            pull(stan_data))[[1]]
        rf_analysis_data = stan_data$analysis_data

        library(brms)
        brms_rf_analysis_data = rf_analysis_data %>%
            select(
                dewormed,
                assigned.treatment, 
                cluster.id, 
                county, 
                cluster.dist.to.pot
            ) %>%
            mutate(sd_dist = cluster.dist.to.pot/sd(cluster.dist.to.pot))


        options(mc.cores = script_options$num_cores)
        rf_model_fit = brm(
            data = brms_rf_analysis_data,
            dewormed ~ sd_dist*assigned.treatment, 
            bernoulli(link = "probit")
        )

        rf_model_fit %>%
            saveRDS(
                file.path(
                    script_options$input_path, 
                    str_interp("param_posterior_draws_dist_fit${fit_version}_REDUCED_FORM_NO_RESTRICT.rds") ))

    } 
} 

if (script_options$from_csv) {
    struct_param_draws = read_csv(
        file.path(
            script_options$input_path, 
            str_interp(
                "param_posterior_draws_dist_fit${fit_version}_${script_options$model}.csv"
            )
        )
    )
    if (!script_options$fit_rf) {
        rf_model_fit = read_rds(
            file.path(
                script_options$input_path, 
                str_interp(
                    "param_posterior_draws_dist_fit${fit_version}_REDUCED_FORM_NO_RESTRICT.rds"
                )
            )
        )
    }
}





## B(z,d):
# \beta is treatment effect
# dist_beta_v is distance cost




## Create estimated demand functions
max_draw = max(struct_param_draws$.draw)
draw_treat_grid = expand.grid(
    draw = 1:min(max_draw, script_options$num_post_draws)
    # treatment = script_options$treatment
)

dist_data = read_rds(
    file.path(
        script_options$data_input_path,
        script_options$data_input_name
    )
)

brms_long_distance_mat = dist_data$brms_long_distance_mat
long_distance_mat = dist_data$long_distance_mat
sd_of_dist = dist_data$sd_of_dist
sd_of_dist %>%
    saveRDS("temp-data/sd_of_dist.rds")
if (run_estimation == TRUE){



pred_functions = map(
    draw_treat_grid$draw,
    ~extract_params(
        param_draws = struct_param_draws,
        private_benefit_treatment = script_options$private_benefit_z,
        visibility_treatment = script_options$visibility_z,
        draw_id = .x,
        dist_sd = sd_of_dist,
        j_id = 1,
        rep_cutoff = script_options$rep_cutoff,
        dist_cutoff = script_options$dist_cutoff, 
        bounds = script_options$bounds,
        mu_rep_type = mu_rep_type,
        suppress_reputation = script_options$suppress_reputation
    ) %>% find_pred_takeup()
)



if (script_options$fit_rf) {

    rf_pred_df = map_dfr(
        script_options$private_benefit_z,
        ~{
            add_epred_draws(
                brms_long_distance_mat %>% mutate(assigned.treatment = .x),
                rf_model_fit,
                ndraws = script_options$num_post_draws
            )
        }
    ) %>%
        rename(
            treatment = assigned.treatment, 
            draw = .draw, 
            demand = .epred
        )  %>%
        ungroup() 

}


subset_long_distance_mat = long_distance_mat %>%
    filter(dist < script_options$dist_cutoff)

tictoc::tic()
plan(multicore, workers = script_options$num_cores)
pred_df = future_imap_dfr(
    pred_functions,
    ~{
        subset_long_distance_mat %>%
            mutate( 
                as_tibble(.x(dist)),
                draw = draw_treat_grid[.y, "draw"]
            )
    },
    .options = furrr_options(
        seed = TRUE
    ),
    .progress = TRUE
)
tictoc::toc()

cutoff_pred_df = future_map_dfr(
    draw_treat_grid$draw,
    ~{
        long_distance_mat %>%
            filter(dist >= script_options$dist_cutoff) %>%
            mutate(
                pred_takeup = NA,
                linear_pred = NA,
                b = NA, 
                mu_rep = NA, 
                delta_v_star = NA, 
                v_star = NA, 
                draw = .x
            )
    }
)

pred_df = bind_rows(
    pred_df,
    cutoff_pred_df
)


structural_demand_df = long_distance_mat %>%
    left_join(
        pred_df %>% 
            select(-dist, -dist_km),
        by = c("index_i", "index_j")
    ) %>%
    mutate(
        across(
            pred_takeup, replace_na, 0.0
        )) %>%
  rename(
    village_i = index_i, 
    pot_j = index_j, 
    demand = pred_takeup) %>%
  group_by(village_i) %>%
  mutate(closest_pot = dist == min(dist)) %>%
  ungroup() %>%
  mutate( 
    model = script_options$model 
  )

if (script_options$fit_rf) {

    rf_demand_df = rf_pred_df %>%
        rename(
            village_i = index_i, 
            pot_j = index_j, 
        ) %>%
        group_by(village_i, closest_pot = dist == min(dist)) %>%
        ungroup() %>%
        mutate(
            model = "REDUCED_FORM_NO_RESTRICT"
        ) %>%
        mutate(
            demand = if_else(
                dist > script_options$dist_cutoff,
                0,
                demand
            )
        ) %>%
        select(-`.row`, -`.chain`, -`.iteration`)
    rf_demand_df$private_benefit_z = script_options$private_benefit_z
    rf_demand_df$visibility_z = script_options$visibility_z
}




structural_demand_df$private_benefit_z = script_options$private_benefit_z
structural_demand_df$visibility_z = script_options$visibility_z






structural_demand_df %>%
    write_csv(
        file.path(
            script_options$output_path, 
            str_interp("pred-demand-dist-fit${fit_version}${append_output}.csv")
        )
    )

if (script_options$fit_rf) {

    rf_demand_df %>%
        write_csv(
            file.path(
                script_options$output_path, 
                str_interp("pred-demand-dist-fit${fit_version}-REDUCED_FORM_NO_RESTRICT.csv")
            )
        )
}

        

}

if (script_options$pred_distance) {

plan(multicore, workers = script_options$num_cores)
pred_functions = map(
    draw_treat_grid$draw,
    ~extract_params(
        param_draws = struct_param_draws,
        private_benefit_treatment = script_options$private_benefit_z,
        visibility_treatment = script_options$visibility_z,
        draw_id = .x,
        dist_sd = sd_of_dist,
        j_id = 1,
        rep_cutoff = script_options$rep_cutoff,
        dist_cutoff = script_options$dist_cutoff, 
        bounds = script_options$bounds,
        mu_rep_type = mu_rep_type,
        suppress_reputation = script_options$suppress_reputation
    ) %>% find_pred_takeup()
)



dist_df = tibble(
        dist = seq(from = 0, to = 2500, by = 100)
) 

sm_df = future_map_dfr(
    pred_functions,
    ~{
        dist_df %>%
            mutate(
                as_tibble(.x(dist))
            ) 
    }, 
    .id = "draw", 
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

sm_df = sm_df %>%
    mutate(
        sm = (-delta + mu_rep_deriv*delta_v_star) / (1 + mu_rep*delta_v_star_deriv)
    ) 

# sm_df %>%
#     group_by(dist) %>%
#     summarise(
#         estimate = mean(sm), 
#         conf.low = quantile(sm, 0.05), 
#         conf.high = quantile(sm, 0.95)
#     ) %>%
#     ggplot(aes(
#         x = dist, 
#         y = estimate, 
#         ymin = conf.low, 
#         ymax = conf.high
#     )) +
#     geom_pointrange()

# sm_df %>%
#     mutate(
#         roc = 100*dnorm(v_star, sd = sqrt(total_error_sd))*sm
#     ) %>%
#     group_by(dist) %>%
#     summarise(
#         estimate = mean(roc), 
#         conf.low = quantile(roc, 0.05), 
#         conf.high = quantile(roc, 0.95)
#     ) %>%
#     ggplot(aes(
#         x = dist, 
#         y = estimate, 
#         ymin = conf.low, 
#         ymax = conf.high
#     )) +
#     geom_pointrange()

sm_df %>%
    write_csv(
        file.path(
            script_options$output_path, 
            str_interp("pred-social-multiplier-fit${fit_version}${append_output}.csv")
        )
    )

# dist_df %>%
#     mutate(
#         sm = (-delta + mu_rep_deriv*delta_v_star) / (1 + mu_rep*delta_v_star_deriv)
#     ) %>%
#     ggplot(aes(
#         x = dist, 
#         y = sm
#     )) +
#     geom_point()


 }
