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
    "),
    args = if (interactive()) "
                            71
                            control
                            control
                            --output-name=cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS
                            --from-csv
                            --num-post-draws=200
                            --rep-cutoff=Inf
                            --dist-cutoff=2500
                            --type-lb=-Inf
                            --type-ub=Inf
                            --num-cores=12
                            --model=STRUCTURAL_LINEAR_U_SHOCKS
                            --num-extra-pots=4
                            --pred-distance
                            --data-input-name=full-experiment.rds
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

mu_rep_type = switch(
    script_options$model,
    STRUCTURAL_LINEAR_U_SHOCKS = 0,
    STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP = 1,
    STRUCTURAL_LINEAR_U_SHOCKS_LINEAR_MU_REP = 2,
    STRUCTURAL_LINEAR_U_SHOCKS_NO_REP = 3
)


models_we_want = c(
    script_options$model
)


# sd_of_dist = sd(rf_analysis_data$cluster.dist.to.pot)




if (script_options$to_csv) {
    struct_model_files = fs::dir_ls(
        script_options$input_path, 
        regex = str_glue("dist_fit{fit_version}_{script_options$model}-.*csv"))
    
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
    if (script_options$fit_rf) {
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

if (run_estimation == TRUE){

dist_data = read_rds(
    file.path(
        script_options$data_input_path,
        script_options$data_input_name
    )
)

brms_long_distance_mat = dist_data$brms_long_distance_mat
long_distance_mat = dist_data$long_distance_mat
sd_of_dist = dist_data$sd_of_dist


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
        mu_rep_type = mu_rep_type
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


append_output = if (!is.null(script_options$output_name)) {
    paste0("-", script_options$output_name)
} else {
    ""
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
df = expand.grid(
    w = seq(from = -2, to = 2, length.out = 100), 
    u_sd = c(0.1, 0.5, 1)
) %>% 
    as_tibble()

df = df %>% 
    mutate(
        ed = analytical_delta(w, u_sd), 
        stan = map2_dbl(
            w, 
            u_sd, 
            ~exposed_funcs$expected_delta(
        w = .x,
        u_sd = .y, 
        total_error_sd = sqrt(1 + .y^2),
        x_r = 1,
        x_i = 1
    )
        )
    )


df %>%
    gather(variable, 
    value, ed:stan) %>%
    ggplot(aes(
        x = w, 
        y = value, 
        colour = variable
    )) +
    geom_point()


w = 0.2
u_sd = 0.3
b = 0.2
mu = 0.4
total_error_sd = sqrt(u_sd^2 + 1)
    analytical_delta(
        w = w,
        u_sd = u_sd 
    )

    exposed_funcs$expected_delta(
        w = w,
        u_sd = u_sd, 
        total_error_sd = sqrt(1 + u_sd^2),
        x_r = 1,
        x_i = 1

    )

   g_f = generate_v_cutoff_fixedpoint(
    b = b, 
    mu = mu,
    total_error_sd = total_error_sd, 
    u_sd = u_sd
   ) 

    stan_soln
    nleqslv(fn = g_f,x =  -b)

    stan_soln = exposed_funcs$find_fixedpoint_solution(
        benefit_cost = b,
        mu_rep = mu,
        total_error_sd = total_error_sd,
        u_sd = u_sd,

        use_u_in_delta = 1,
        alg_sol_f_tol = 0.001,
        alg_sol_max_steps = 1e9L,
        alg_sol_rel_tol = 0.0000001
    )




  stan_mu_df %>%
    ggplot(aes(
        x = cluster.dist.to.pot, 
        y = median, 
        colour = assigned.treatment
    )) +
    geom_point() +
    geom_line()



comp_df %>%
    filter(variable == "pred_takeup") %>%
    ggplot(aes(
        x = dist, 
        y = mean_est, 
        colour = type
    )) +
    geom_point() +
    facet_wrap(~treatment)


comp_df %>% 
    filter(variable == "v_star") %>%
    filter(dist < 3500) %>%
    ggplot(aes(
        x = dist, 
        y = mean_est, 
        colour = type
    )) +
    geom_point(size = 2) +
    facet_wrap(~treatment, scales = "free")

comp_df %>%
    filter(variable == "mu_rep") %>%
    filter(dist < 3500) %>%
    ggplot(aes(
        x = dist, 
        y = median_est, 
        colour = type
    )) +
    geom_point(size = 2) +
    facet_wrap(~treatment, scales = "free")

    summ_df %>%
        filter(
            variable == "mu_rep"
        ) %>%
        filter(dist < 2500) %>%
        ggplot(aes(
            x = dist, 
            y = median_est, 
            colour = treatment
        ))  +
        geom_point() +
        geom_line() +
        geom_hline(
            yintercept = 0.1
        )




comp_df %>%
    filter(variable == "pred_takeup") 


    pred_dist_df %>%
        write_csv(
            file.path(
                script_options$output_path, 
                str_interp("all-treat-demand-dist-fit${fit_version}-${script_options$model}.csv")
            )
        )



comp_df %>%
    filter(
        variable %in% c(
            "pred_takeup", 
            "v_star"
        )
    ) %>%
    ggplot(aes(
        x = dist, 
        y = mean_est, 
        colour = type
    )) +
    geom_point() +
    facet_grid(treatment~variable)


# high_bound = Inf
# low_bound = 2

# n_posts = 10

# agg_data = function(
#     high_bound, 
#     low_bound, 
#     n_posts, 
#     treatment, 
#     mu_rep_type){
#     high_bound_pred_fs = map(
#         1:n_posts,
#         ~extract_params(
#             param_draws = struct_param_draws,
#             treatment = treatment,
#             draw_id = .x,
#             dist_sd = dist_sd,
#             j_id = 1,
#             rep_cutoff = script_options$rep_cutoff,
#             dist_cutoff = script_options$dist_cutoff, 
#             bounds = c(-high_bound, high_bound) ,
#             mu_rep_type = mu_rep_type
#         ) %>% find_pred_takeup()
#     )

#     low_bound_pred_fs = map(
#         1:n_posts,
#         ~extract_params(
#                 param_draws = struct_param_draws,
#                 treatment = treatment,
#                 draw_id = .x,
#                 dist_sd = dist_sd,
#                 j_id = 1,
#                 rep_cutoff = script_options$rep_cutoff,
#                 dist_cutoff = script_options$dist_cutoff, 
#                 bounds = c(-low_bound, low_bound),
#                 mu_rep_type = mu_rep_type
#             ) %>% find_pred_takeup()
#     )
#     from_d = 10
#     to_d = 10000
#     len_out = 100

#     dist_df = tibble(
#         distance = seq(from_d, to_d, length.out = len_out)
#     )



#     high_dfs = imap_dfr(
#         high_bound_pred_fs,
#         ~ {
#             dist_df %>%
#                 mutate(as_tibble(.x(distance))) %>%
#                 mutate(draw = .y) 
#         }
#     ) %>%
#         mutate(bound = high_bound)
#     low_dfs = imap_dfr(
#         low_bound_pred_fs,
#         ~ {
#             dist_df %>%
#                 mutate(as_tibble(.x(distance))) %>%
#                 mutate(draw = .y) 
#         }
#     ) %>%
#         mutate(bound = low_bound)


#     wide_comp_dfs = bind_rows(
#         high_dfs,
#         low_dfs
#     ) %>%
#     mutate(bound = factor(bound)) %>%
#     mutate(treatment = treatment)
#     return(wide_comp_dfs)
# }


# wide_comp_bracelet_rep_dfs = agg_data(
#     high_bound = Inf,
#     low_bound = 2,
#     n_posts = 20,
#     treatment = "bracelet",
#     mu_rep_type = mu_rep_type,
# ) 

# wide_comp_bracelet_no_rep_dfs = agg_data(
#     high_bound = Inf,
#     low_bound = 2,
#     n_posts = 20,
#     treatment = "bracelet",
#     mu_rep_type = mu_rep_type
# )  


# wide_comp_bracelet_dfs = bind_rows(
#     wide_comp_bracelet_rep_dfs,
#     wide_comp_bracelet_no_rep_dfs
# )



# wide_comp_control_no_rep_dfs = agg_data(
#     high_bound = Inf,
#     low_bound = 2,
#     n_posts = 20,
#     treatment = "control",
#     mu_rep_type = mu_rep_type,
# )  

# wide_comp_control_rep_dfs = agg_data(
#     high_bound = Inf,
#     low_bound = 2,
#     n_posts = 20,
#     treatment = "control",
#     mu_rep_type = mu_rep_type
# )  

# wide_comp_control_dfs = bind_rows(
#     wide_comp_control_no_rep_dfs,
#     wide_comp_control_rep_dfs
# )

# long_comp_bracelet_df = wide_comp_bracelet_dfs %>%
#     gather(
#         variable, value, 
#         -distance, 
#         -bound, 
#         -draw, 
#         -treatment
#     )

# long_comp_control_df = wide_comp_control_dfs %>%
#     gather(
#         variable, value, 
#         -distance, 
#         -bound, 
#         -draw, 
#         -treatment
#     )

# long_comp_df = bind_rows(
#     long_comp_bracelet_df,
#     long_comp_control_df
# )

# long_comp_df

# long_comp_df %>%
#     filter(treatment == "bracelet") %>%
#     ggplot(aes(
#         x = distance, 
#         y = value, 
#         colour = bound, 
#         group = interaction(bound, draw)
#     )) +
#     geom_line() +
#     facet_wrap(~variable, scales = 'free') +
#     labs(
#         title = str_glue("Bracelet Posterior Draws Over Distance"), 
#         subtitle = str_glue("{script_options$model}")
#     ) +
#     theme_bw() +
#     theme(legend.position = "bottom") +
#     labs(colour = "V Bound:")

# ggsave(
#     str_glue(
#         "temp-plots/bracelet-model-predictions-dist-fit{fit_version}-{script_options$model}.png"
#     ),
#     width = 8,
#     height = 6,
#     dpi = 500
# )

# long_comp_df %>%
#     filter(treatment == "control") %>%
#     ggplot(aes(
#         x = distance, 
#         y = value, 
#         colour = bound, 
#         group = interaction(bound, draw)
#     )) +
#     geom_line() +
#     facet_wrap(~variable, scales = 'free') +
#     labs(
#         title = str_glue("Control Posterior Draws Over Distance"), 
#         subtitle = str_glue("{script_options$model}")
#     ) +
#     theme_bw() +
#     theme(legend.position = "bottom") +
#     labs(colour = "V Bound:")
# ggsave(
#     str_glue(
#         "temp-plots/control-model-predictions-dist-fit{fit_version}-{script_options$model}.png"
#     ),
#     width = 8,
#     height = 6,
#     dpi = 500
# )


# long_comp_df %>%
#     filter(fct_match(bound, high_bound)) %>%
#     ggplot(aes(
#         x = distance, 
#         y = value, 
#         colour = treatment, 
#         group = interaction(treatment, draw)
#     )) +
#     geom_line() +
#     facet_wrap(~variable, scales = 'free') +
#     labs(
#         title = str_glue("Comparison Posterior Draws Over Distance"), 
#         subtitle = str_glue("{script_options$model}, Unbounded")
#     ) +
#     theme_bw() +
#     theme(legend.position = "bottom") +
#     labs(colour = "Treatment:")


# ggsave(
#     str_glue(
#         "temp-plots/comp-model-unbounded-predictions-dist-fit{fit_version}-{script_options$model}.png"
#     ),
#     width = 8,
#     height = 6,
#     dpi = 500
# )


# long_comp_df %>%
#     filter(fct_match(bound, high_bound)) %>%
#     filter(variable == "pred_takeup") %>%
#     filter(distance <= 2500) %>%
#     ggplot(aes(
#         x = distance, 
#         y = value, 
#         colour = interaction(treatment, suppress_rep), 
#         group = interaction(treatment, suppress_rep, draw)
#     )) +
#     geom_line() +
#     # facet_wrap(
#     #     ~suppress_rep, 
#     #     scales = 'free', 
#     #     ncol = 1
#     #     ) +
#     labs(
#         title = str_glue("Comparison Posterior Draws Over Distance"), 
#         subtitle = str_glue("{script_options$model}, Unbounded")
#     ) +
#     theme_bw() +
#     theme(legend.position = "bottom") +
#     labs(colour = "Treatment:")

# long_comp_df %>%
#     filter(fct_match(bound, high_bound)) %>%
#     filter(variable == "pred_takeup") %>%
#     filter(distance <= 2500) %>%
#     ggplot(aes(
#         x = distance, 
#         y = value, 
#         colour = treatment, 
#         group = interaction(treatment, draw)
#     )) +
#     geom_line() +
#     facet_wrap(
#         ~suppress_rep, 
#         # scales = 'free', 
#         ncol = 2
#         ) +
#     labs(
#         title = str_glue("Comparison Posterior Draws Over Distance"), 
#         subtitle = str_glue("{script_options$model}, 2.5km cutoff - Suppress Reputation FALSE/TRUE")
#     ) +
#     theme_bw() +
#     theme(legend.position = "bottom") +
#     labs(colour = "Treatment:")

# ggsave(
#     str_glue(
#         "temp-plots/close-comp-model-bounded-predictions-dist-fit{fit_version}-{script_options$model}.png"
#     ),
#     width = 8,
#     height = 6,
#     dpi = 500
# )


# long_comp_df %>%
#     filter(fct_match(bound, low_bound)) %>%
#     filter(suppress_rep == FALSE) %>%
#     ggplot(aes(
#         x = distance, 
#         y = value, 
#         colour = treatment, 
#         group = interaction(treatment, draw)
#     )) +
#     geom_line() +
#     facet_wrap(~variable, scales = 'free') +
#     labs(
#         title = str_glue("Comparison Posterior Draws Over Distance"), 
#         subtitle = str_glue("{script_options$model}, Bounded")
#     ) +
#     theme_bw() +
#     theme(legend.position = "bottom") +
#     labs(colour = "Treatment:")


# ggsave(
#     str_glue(
#         "temp-plots/comp-model-bounded-predictions-dist-fit{fit_version}-{script_options$model}.png"
#     ),
#     width = 8,
#     height = 6,
#     dpi = 500
# )


# long_comp_df %>%
#     filter(
#         treatment == "bracelet"
#     ) %>%
#     filter(fct_match(bound, high_bound)) %>%
#     ggplot(aes(
#         x = distance,
#         y = value, 
#         colour = suppress_rep, 
#         group = interaction(suppress_rep, draw)
#     )) +
#     geom_line() +
#     facet_wrap(~variable, scales = 'free') +
#     labs(
#         title = str_glue("Comparison Posterior Draws Over Distance"), 
#         subtitle = str_glue("{script_options$model}, Unbounded - Comparing Rep vs No Rep")
#     ) +
#     theme_bw() +
#     theme(legend.position = "bottom") +
#     labs(colour = "Suppress Reputational Returns:")
# ggsave(
#     str_glue(
#         "temp-plots/comp-model-unbounded-rep-comp-predictions-dist-fit{fit_version}-{script_options$model}.png"
#     ),
#     width = 8,
#     height = 6,
#     dpi = 500
# )
# # long_comp_df %>%
# #     filter(fct_match(bound, high_bound)) %>%
# #     filter(distance < 5000) %>%
# #     ggplot(aes(
# #         x = distance, 
# #         y = value, 
# #         colour = bound, 
# #         group = interaction(bound, draw)
# #     )) +
# #     geom_point() +
# #     geom_line() +
# #     facet_wrap(~variable, scales = 'free') +
# #     theme_bw()

# # wide_comp_dfs %>%
# #     filter(fct_match(bound, low_bound)) %>%
# #     ggplot(aes(
# #         x = v_star, 
# #         y = delta_v_star, 
# #         group = draw, 
# #         colour = distance
# #     )) +
# #     geom_line() +
# #     geom_point()

# # test_df_1 = tibble(
# #     x = seq(from = from_d, to = to_d, length.out = len_out),
# # ) %>%
# # mutate(as_tibble(high_bound_pred_f(x))) %>%
# # gather(variable, value, -x) %>%
# # mutate(bound = high_bound)


# #     test_df_2 = tibble(
# #         x = seq(from = from_d, to = to_d, length.out = len_out),
# #     ) %>%
# #     mutate(as_tibble(low_bound_pred_f(x))) %>%
# #     gather(variable, value, -x) %>%
# #     mutate(bound = low_bound)

# #     comp_test_df = bind_rows(
# #         test_df_1,
# #         test_df_2
# #     ) %>%
# #     mutate(bound = factor(bound))



# #     comp_test_df %>%
# #         filter(x < 4000) %>%
# #         ggplot(aes(
# #             x = x, 
# #             y = value, 
# #             colour = bound
# #         )) +
# #         geom_line() +
# #         geom_point() +
# #         facet_wrap(~variable, scales = "free")  

# #     comp_test_df %>%
# #         filter(fct_match(bound, low_bound)) %>%
# #         ggplot(aes(
# #             x = x, 
# #             y = value, 
# #             colour = bound
# #         )) +
# #         geom_line() +
# #         facet_wrap(~variable, scales = "free")  
# # comp_test_df %>%
# #    spread(variable, value)  %>%
# #    ggplot(aes(
# #     x = x, 
# #     y = v_star, 
# #     colour = bound
# #    )) +
# #    geom_line() +
# #    geom_point() +
# #    theme_bw() +
# #    facet_wrap(~bound) 


# # inf_x = comp_test_df %>%
# #    spread(variable, value)  %>%
# #    filter(!is.finite(delta_v_star)) %>%
# #    pull(x)


# # dist_x = seq(from = from_d, to = to_d, length.out = len_out)
# # pred_manual = map_dfr(dist_x, ~low_bound_pred_f(.x) %>% as_tibble())

# # manual_df = tibble(
# #     x = dist_x,
# #     as_tibble(pred_manual)
# # )

# # manual_df %>%
# #     ggplot(aes(
# #         x = x, 
# #         y = delta_v_star
# #     )) +
# #     geom_point()


# # manual_df %>%
# #     filter(abs(v_star) < low_bound) %>%
# #     gather(variable, value, -x) %>%
# #     ggplot(aes(
# #         x = x, 
# #         y = value
# #     )) +
# #     geom_point() +
# #     geom_line() +
# #     facet_wrap(~variable, scales = "free")


# # manual_df %>%
# #     filter(abs(v_star) < low_bound) %>%
# #     ggplot(aes(
# #         x = v_star,
# #         y = delta_v_star, 
# #         colour = x
# #     )) +
# #     geom_point()

# # manual_df %>%
# #     filter(!is.finite(delta_v_star))




# # low_bound_pred_f(inf_x)



# # low_bound_pred_f(inf_x)



# # comp_test_df %>%
# #    spread(variable, value)  %>%
# #    ggplot(aes(
# #     x = x, 
# #     y = delta_v_star, 
# #     colour = bound
# #    )) +
# #    geom_line() +
# #    geom_point() +
# #    theme_bw() +
# #    facet_wrap(~bound) 



# # comp_test_df %>%
# #    spread(variable, value)  %>%
# #    ggplot(aes(
# #     x = v_star, 
# #     y = delta_v_star, 
# #     colour = bound
# #    )) +
# #    geom_line() +
# #    geom_point() +
# #    theme_bw() +
# #    facet_wrap(~bound) +
# #    xlim(-10, 10)



# # pred_df %>%
# #     filter(treatment == "bracelet") %>%
# #     ggplot(aes(
# #         x = v_star,
# #         y = delta_v_star, 
# #         group = draw
# #     )) +
# #     geom_line()

# # pred_df %>%
# #     ggplot(aes(
# #         x = distance, 
# #         y = pred_takeup, 
# #         group = interaction(draw, treatment), 
# #         colour = treatment
# #     )) +
# #     geom_line()

# # pred_df %>%
# #     gather(variable, value, -distance, -treatment, -draw) %>%
# #     ggplot(aes(
# #         x = distance, 
# #         y = value, 
# #         colour = treatment, 
# #         group = interaction(draw, treatment)
# #     )) +
# #     geom_line() +
# #     facet_wrap(~variable, scales = "free") 

# #     ggplot(aes(
# #         x = distance
# #     ))

# #     pred_functions[[1]](500)

# #     pred_functions[[1]](c(500, 10000))

# #     struct_param_draws
# #     test_params = extract_params(struct_param_draws,
# #             treatment = "bracelet",
# #             draw_id = 1,
# #             dist_sd = dist_sd,
# #             j_id = 1,
# #             rep_cutoff = script_options$rep_cutoff,
# #             dist_cutoff = script_options$dist_cutoff, 
# #             bounds = c(-Inf, Inf) )

# #     test_func = find_pred_takeup(test_params)
# #     test_params

# #     test_func(500)

# #     test_func(10)


# #     test_func(c(500, 10))

# #     comp_test_df %>%
# #         filter(x < 2000) %>%
# #         ggplot(aes(
# #             x = x, 
# #             y = value, 
# #             colour = type
# #         )) +
# #         geom_line() +
# #         facet_wrap(~variable, scales = "free") 

# #     unbound_pred_df = future_imap_dfr(unbound_pred_functions,
# #         ~tibble(
# #             distance = seq(from = 0, to = 10000, length.out = 100)) %>% 
# #                 mutate(
# #                     as_tibble(.x(distance)), 
# #                     treatment = draw_treat_grid[.y, "treatment"],
# #                     draw = draw_treat_grid[.y, "draw"]),
# #         .options = furrr_options(
# #             globals = c(
# #             "draw_treat_grid",
# #             "script_options",
# #             "long_distance_mat",
# #             "extract_params", 
# #             "struct_param_draws", 
# #             "dist_sd", 
# #             "find_pred_takeup",
# #             ".find_pred_takeup",
# #             "calculate_mu_rep",
# #             "calculate_belief_latent_predictor",
# #             "find_v_star",
# #             "generate_v_cutoff_fixedpoint",
# #             "analytical_delta",
# #             "analytical_delta_bounded",
# #             "analytical_conv_Fw",
# #             "analytical_conv_Fw_part",
# #             "stan_owen_t"
# #             ),
# #             packages = c("dplyr", "nleqslv")
# #         ),
# #         .progress = TRUE
# #     )  

# #     pred_df = future_imap_dfr(pred_functions,
# #         ~tibble(
# #             distance = seq(from = 0, to = 10000, length.out = 100)) %>% 
# #                 mutate(
# #                     as_tibble(.x(distance)), 
# #                     treatment = draw_treat_grid[.y, "treatment"],
# #                     draw = draw_treat_grid[.y, "draw"]),
# #         .options = furrr_options(
# #             globals = c(
# #             "draw_treat_grid",
# #             "script_options",
# #             "long_distance_mat",
# #             "extract_params", 
# #             "struct_param_draws", 
# #             "dist_sd", 
# #             "find_pred_takeup",
# #             ".find_pred_takeup",
# #             "calculate_mu_rep",
# #             "calculate_belief_latent_predictor",
# #             "find_v_star",
# #             "generate_v_cutoff_fixedpoint",
# #             "analytical_delta",
# #             "analytical_delta_bounded",
# #             "analytical_conv_Fw"),
# #             packages = c("dplyr", "nleqslv")
# #         ),
# #         .progress = TRUE
# #     )  



# #     treat_levels = c(
# #         "bracelet",
# #         "ink",
# #         "calendar",
# #         "control"
# #     )
# #     pred_df = pred_df %>%
# #         mutate(treatment = factor(treatment, levels = treat_levels))

# #     pred_df %>%
# #         write_csv(
# #             file.path(
# #                 script_options$output_path, 
# #                 str_interp("pred_posterior_draws_dist_fit${fit_version}_STRUCTURAL_LINEAR_U_SHOCKS.csv")
# #             )
# #         )

# #     summ_pred_df = pred_df %>%
# #         group_by(distance, treatment) %>%
# #         summarise_all(median) %>%
# #         mutate(treatment = factor(treatment, levels = treat_levels))

# #     pred_df %>%
# #         ggplot(aes(
# #             x = distance, 
# #             y = pred_takeup, 
# #             group = interaction(draw, treatment), 
# #             colour = treatment
# #         )) +
# #         geom_line(
# #             alpha = 0.1
# #         )   +
# #         geom_line(
# #             data = summ_pred_df, 
# #             inherit.aes = FALSE,
# #             aes(
# #                 x = distance, 
# #                 y = pred_takeup, 
# #                 group = interaction(treatment), 
# #                 colour = treatment), 
# #             size = 1
# #         ) +
# #         theme_minimal() +
# #         theme(legend.position = "bottom") +
# #         labs(
# #             title = "Predicted Takeup", 
# #             y = "Pred Takeup", 
# #             x = "Distance (m)"
# #         )

# #     ggsave(
# #         "temp-plots/pred-takeup.png", 
# #         width = 10,
# #         height = 10,
# #         dpi = 500
# #         )

# #     pred_df %>%
# #         ggplot(aes(
# #             x = distance, 
# #             y = pred_takeup, 
# #             group = interaction(draw, treatment),
# #             colour = treatment
# #         )) +
# #         geom_line(
# #             alpha = 0.3
# #         )  +
# #         geom_line(
# #             data = summ_pred_df, 
# #             inherit.aes = FALSE,
# #             aes(
# #                 x = distance, 
# #                 y = pred_takeup, 
# #                 group = interaction(treatment)), 
# #             colour = "black",
# #             size = 1
# #         ) +
# #         facet_wrap(~treatment) +
# #         theme_bw() +
# #         guides(color = "none") +
# #         labs(
# #             title = "Pred Takeup", 
# #             y = "Takeup", 
# #             x = "Distance (m)"
# #         )
        
# #     ggsave(
# #         "temp-plots/pred-takeup-facet.png", 
# #         width = 10,
# #         height = 10,
# #         dpi = 500
# #         )

        





# #     filter_structural_fit = function(data){
# #         data %>%
# #             filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) 
# #     }

# #     treat_levels = c(
# #         "bracelet",
# #         "ink",
# #         "calendar", 
# #         "control"
# #     )

# #     dist_fit_data %>%
# #         filter_structural_fit() %>%
# #         select(cluster_takeup_prop) %>%
# #         unnest() %>%
# #         mutate(
# #             treatment = factor(assigned_treatment, levels = treat_levels), 
# #             type = "stan", 
# #             pred_takeup = mean_est, 
# #             distance = roc_distance
# #         ) 

# #     dist_fit_data %>%
# #         filter_structural_fit() %>%
# #         select(cluster_takeup_prop) %>%
# #         unnest() %>%
# #         mutate(
# #             treatment = factor(assigned_treatment, levels = treat_levels), 
# #             type = "stan", 
# #             pred_takeup = per_0.5, 
# #             distance = roc_distance
# #         ) %>%
# #         bind_rows(
# #             summ_pred_df %>%
# #                 filter(distance < 2500) %>%
# #                 mutate(type = "extracted")
# #         ) %>%
# #         select( 
# #             treatment,
# #             distance, 
# #             pred_takeup, 
# #             type
# #         ) %>%
# #         ggplot(aes(
# #             x = distance, 
# #             y = pred_takeup, 
# #             colour = type 
# #         )) + 
# #         geom_point() +
# #         geom_line() +
# #         facet_wrap(~treatment) +
# #         theme_bw() +
# #         theme(legend.position = "bottom") +
# #         labs(
# #             title = "Comparing Stan vs R Posterior Medians"
# #         )

# #     ggsave(
# #         "temp-plots/temp-p-comp.png",
# #         width = 10,
# #         height = 10,
# #         dpi = 500
# #     )
        
        
# #        pred_df  %>%
# #         ggplot(aes(
# #             x = distance, 
# #             colour = treatment, 
# #             y = delta_v_star, 
# #             group = interaction(draw, treatment)
# #         )) +
# #         geom_line()
        
# #     pred_df    
        
# #     summ_pred_df  %>%
# #         mutate(
# #             r_by_d = mu_rep*delta_v_star
# #         ) %>%
# #         ggplot(aes(
# #             x = distance, 
# #             colour = treatment, 
# #             y = pred_takeup
# #         )) +
# #         geom_line() +
# #         geom_point()

# #     cowplot::plot_grid(
# #         dist_fit_data %>%
# #             filter_structural_fit() %>%
# #             select(cluster_takeup_prop) %>%
# #             unnest() %>%
# #             mutate(
# #                 assigned_treatment = factor(assigned_treatment, levels = treat_levels)
# #             ) %>%
# #             ggplot(aes(
# #                 x = roc_distance, 
# #                 y = per_0.5, 
# #                 colour = assigned_treatment
# #             )) +
# #             geom_point() +
# #             geom_line(),

# #         summ_pred_df  %>%
# #             filter(distance <= 2500) %>%
# #             mutate(
# #                 r_by_d = mu_rep*delta_v_star
# #             ) %>%
# #             ggplot(aes(
# #                 x = distance, 
# #                 colour = treatment, 
# #                 y = pred_takeup
# #             )) +
# #             geom_line() +
# #             geom_point() 
# #     ) 

# #     # cowplot::plot_grid(
# #     #     dist_fit_data %>%
# #     #         filter_structural_fit() %>%
# #     #         select(cluster_rep_return) %>%
# #     #         unnest() %>%
# #     #         mutate(
# #     #             assigned_treatment = factor(assigned_treatment, levels = treat_levels)
# #     #         ) %>%
# #     #         ggplot(aes(
# #     #             x = roc_distance, 
# #     #             y = per_0.5, 
# #     #             colour = assigned_treatment
# #     #         )) +
# #     #         geom_point() +
# #     #         geom_line(),
# #     #     summ_pred_df  %>%
# #     #         filter(distance <= 2500) %>%
# #     #         mutate(
# #     #             r_by_d = mu_rep*delta_v_star
# #     #         ) %>%
# #     #         ggplot(aes(
# #     #             x = distance, 
# #     #             colour = treatment, 
# #     #             y = r_by_d
# #     #         )) +
# #     #         geom_line() +
# #     #         geom_point()
# #     # )

# #     pred_df  %>%
# #         # filter(v_star > -5) %>%
# #         ggplot(aes(
# #             x = v_star, 
# #             group = interaction(draw, treatment),
# #             colour = treatment, 
# #             y = delta_v_star
# #         )) +
# #         geom_line(alpha = 0.2) +
# #         geom_point(alpha = 0.2) +
# #         theme_minimal() +
# #         theme(legend.position = "bottom") +
# #         labs(
# #             title = "Delta(v^*) vs v^*"
# #         )

# #     ggsave(
# #         "temp-plots/all-delta-v-vs-v.png",
# #         width = 10,
# #         height = 10,
# #         dpi = 500
# #     )



# #     pred_df  %>%
# #         filter(abs(v_star) < 2) %>%
# #         # filter(v_star > -5) %>%
# #         ggplot(aes(
# #             x = v_star, 
# #             group = interaction(draw, treatment),
# #             colour = treatment, 
# #             y = delta_v_star
# #         )) +
# #         geom_line(alpha = 0.2) +
# #         geom_point(alpha = 0.2) +
# #         theme_minimal() +
# #         theme(legend.position = "bottom") +
# #         labs(
# #             title = "Delta(v^*) vs v^*"
# #         )

# #     ggsave(
# #         "temp-plots/subset-delta-v-vs-v.png",
# #         width = 10,
# #         height = 10,
# #         dpi = 500
# #     )

# #     summ_pred_df  %>%
# #         # filter(v_star > -5) %>%
# #         ggplot(aes(
# #             x = v_star, 
# #             colour = treatment, 
# #             y = delta_v_star
# #         )) +
# #         geom_line() +
# #         geom_point() +
# #         theme_minimal() +
# #         theme(legend.position = "bottom") +
# #         labs(
# #             title = "Delta(v^*) vs v^*"
# #         )


# #     ggsave(
# #         "temp-plots/delta-v-vs-v.png",
# #         width = 10,
# #         height = 10,
# #         dpi = 500
# #     )


# #     summ_pred_df  %>%
# #         ggplot(aes(
# #             x = distance, 
# #             colour = treatment, 
# #             y = v_star
# #         )) +
# #         geom_line() +
# #         geom_point()


# #     dist_fit_data %>%
# #         filter_structural_fit()  %>%
# #         colnames()

# #     dist_fit_data %>%
# #         filter_structural_fit()  %>%
# #         select(cluster_w_cutoff) %>%
# #         unnest() %>%
# #         ggplot(aes(
# #             x = roc_distance, 
# #             y = mean_est, 
# #             colour = assigned_treatment
# #         )) +
# #         geom_point() +
# #         geom_line()


# #     pred_df  %>%
# #         filter(v_star > -5) %>%
# #         ggplot(aes(
# #             x = distance, 
# #             colour = treatment, 
# #             y = delta_v_star
# #         )) +
# #         geom_line() +
# #         geom_point()


# #     pred_df  %>%
# #         filter(v_star > -5) %>%
# #         filter(distance <= 2.5) %>%
# #         # filter(treatment %in% c("bracelet", "control")) %>%
# #         ggplot(aes(
# #             x = distance, 
# #             colour = treatment, 
# #             y = pred_takeup
# #         )) +
# #         geom_line() +
# #         # geom_point() +
# #         geom_hline( 
# #             yintercept = c(0.5, 0.2), 
# #             linetype = "longdash"
# #         ) 


# #     dist_prob_plot <- dist_fit_data %>% 
# #         filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural"))  %>%
# #         select(cluster_takeup_prop) %>%
# #         unnest() %>%
# #         ggplot(aes(
# #             x = roc_distance, 
# #             y = mean_est, 
# #             colour = assigned_treatment
# #         )) +
# #         geom_line()


# #     dist_prob_plot  



# #     pred_df  %>%
# #         filter(v_star > -5) %>%
# #         ggplot(aes(
# #             x = distance, 
# #             colour = treatment, 
# #             y = v_star
# #         )) +
# #         geom_line() +
# #         geom_point()

# #     pred_df %>%
# #         filter(distance == max(distance))



# #     pred_df %>%
# #         ggplot(aes(
# #             x = distance,
# #             y = y_hat,
# #             colour = treatment
# #         )) +
# #         geom_line()


# #     analysis_data %>%
# #         select(standard_cluster.dist.to.pot) %>%
# #         summary()



# #     ed = find_pred_takeup(
# #         distance = 0.5, 
# #         beta = betas[[1]], 
# #         mu_beta = mu_betas[[1]], 
# #         base_mu_rep = 0.121, 
# #         total_error_sd = 1.2,
# #         u_sd =0)

# #     pred_df = tibble(
# #         dist = seq(from = -2, to = 2, length.out = 100)
# #     ) %>%
# #         mutate(y_hat = ed(dist))


# #     pred_df %>%
# #         ggplot(aes(
# #             x = dist,
# #             y = y_hat
# #         )) +
# #         geom_point()



# #     find_v_star(
# #         distance = 0.5, 
# #         beta = betas[[1]], 
# #         mu_beta = mu_betas[[1]], 
# #         base_mu_rep = 0.121, 
# #         total_error_sd = 1.2,
# #         u_sd =0)








# #     param_summary %>%
# #         head(11)



# #     variables = fitted_draw_df %>%
# #         colnames()

# #     variables[!str_detect(variables, "cf")]

# #     model_fit$summary("beta_dist_v")

# #     variables(model_fit)
# #     model_fit$variables()


# #     model_fit$metadata() %>% head(20)

# #     str(model_fit)


# #     dist_df = imap_dfr(
# #         mu_betas,
# #         ~expand.grid(
# #             dist = seq(from = 0, to = 4, length.out = 100)
# #         ) %>%
# #         tibble() %>%
# #         mutate(
# #             mu_rep = calculate_mu_rep(dist = dist, 
# #                                     base_mu_rep = clean_param_summary %>%  
# #                                         filter(variable == "base_mu_rep") %>%
# #                                         pull(mean), 
# #                                     mu_beliefs_effect = 1, 
# #                                     beta = .x[1], 
# #                                     dist_beta = .x[2], 
# #                                     beta_control = mu_betas[[1]][1]),
# #             treatment = treatments[.y]
# #         )
# #     ) %>%
# #         mutate(orig_dist = dist*dist_sd)


# #     dist_df %>%
# #         # filter(treatment %in% c("bracelet", "control")) %>%
# #         ggplot(aes(
# #             x = orig_dist, 
# #             y = mu_rep, 
# #             colour = treatment
# #         )) +
# #         geom_point() +
# #         labs(title = "Mu(z, d) from model") +
# #         theme_bw()


 }
