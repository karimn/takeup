
#!/usr/bin/Rscript

script_options = docopt::docopt(
    stringr::str_glue("Usage:
    predict-takeup-for-optim.R <fit-version> [options] [ --from-csv | --to-csv ]

    Takes posterior fits and calculates estimated takeup for a given distance.


    Options:
        --from-csv  Load parameter draws from a csv file
        --to-csv   Load stanfit object and write draws to csv
        --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
        --output-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
    "),
    args = if (interactive()) "66 --from-csv"
)


library(posterior)
library(tidyverse)
library(tidybayes)
library(broom)
library(rstan)

library(nleqslv)
library(cmdstanr)
library(econometr)
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))


source(file.path("multilvlr", "multilvlr_util.R"))

fit_version = script_options$fit_version

models_we_want = c(
  "STRUCTURAL_LINEAR_U_SHOCKS"
)

wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
kenya.proj4 <- "+proj=utm +zone=36 +south +ellps=clrk80 +units=m +no_defs"

rct.schools.data <- read_rds(file.path("data", "takeup_rct_schools.rds"))
rct.cluster.selection <- read_rds(file.path("data", "rct_cluster_selection_2.0.rds"))
cluster.strat.data <- read_rds(file.path("data", "takeup_processed_cluster_strat.rds"))

load(file.path("data", "takeup_village_pot_dist.RData"))

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)

nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

analysis_data <- monitored_nosms_data

dist_sd <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
    summarise(ed_sd = sd(cluster.dist.to.pot)) %>%
    pull()

## Fit Loading
# always load processed fit as we can grab analysis_data used for model and 
# calculate sd to unstandardise
load(file.path("temp-data", str_interp("processed_dist_fit${fit_version}_lite.RData")))


stan_data = (dist_fit_data %>%
    slice(1) %>% # just take first model as only using for analysis_data
    pull(stan_data))[[1]]
rf_analysis_data = stan_data$analysis_data
if (script_options$to_csv) {
    model_files = fs::dir_ls(script_options$input_path, regex = str_glue("dist_fit{fit_version}_STRUCTURAL_LINEAR_U_SHOCKS.*\\.csv"))
    model_fit = as_cmdstan_fit(model_files)
    # N.B. this is very hardcoded for vanilla structural model
    param_draws = model_fit %>%
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
    rm(model_fit) 
    gc()
    # N.B. we use input path for `to_csv` output and reserve output_path for 
    # output later in script
    param_draws %>%
        write_csv(file.path(script_options$input_path, str_interp("param_posterior_draws_dist_fit${fit_version}_STRUCTURAL_LINEAR_U_SHOCKS.csv") ))
} 

if (script_options$from_csv) {
    param_draws = read_csv(
        file.path(
            script_options$input_path, 
            str_interp(
                "param_posterior_draws_dist_fit${fit_version}_STRUCTURAL_LINEAR_U_SHOCKS.csv"
            )
        )
    )
}



## B(z,d):
# \beta is treatment effect
# dist_beta_v is distance cost





## Calculating \mu(z,d)
# centered_cluster_beta_1ord, centered_cluster_dist_beta_1ord

# beta_1ord and dist_beta_1ord go into calculate mu rep

# vector calculate_mu_rep(array[] int treatment_ids, vector dist,
#                         real base_mu_rep, real mu_beliefs_effect,
#                         matrix design_matrix,
#                         matrix beta, matrix dist_beta) {
#   vector[rows(beta)] beliefs_latent = calculate_beliefs_latent_predictor(design_matrix[treatment_ids], beta, dist_beta, dist);
  
#   return base_mu_rep * exp(mu_beliefs_effect * (beliefs_latent - beta[, 1])); // Remove intercept 
# }

# vector calculate_beliefs_latent_predictor(matrix design_matrix, matrix beta, matrix dist_beta, vector dist) {
#   if (rows(design_matrix) == 1) {
#     return (beta * design_matrix[1]') + ((dist_beta * design_matrix[1]') .* dist);
#   } else {
#     return rows_dot_product(design_matrix, beta) + (rows_dot_product(design_matrix, dist_beta) .* dist);
#   }
# }




treatments = c(
    "control",
    "ink", 
    "calendar", 
    "bracelet"
)
#' Calculate latent predictor for beliefs model
#'
#' Typically takes as input the distance, the coef on beliefs for each treatment 
#' and the coef on beliefs for each distance x treatment
calculate_belief_latent_predictor = function(beta, dist_beta, dist) {
    val =  beta + dist*dist_beta 
    return(val)
}

#' Calculate Visibility \mu(z,d)
#' 
#' Takes in distance, base_mu_rep which act as control, mu_beliefs_effect (=1),
#' betas from beliefs model, dist_beta from beliefs model and the control coef
#' from the beliefs model.
#' (The exponential model defines \mu_0 exp(\beta_dist d)) as the control effect
calculate_mu_rep = function(dist, 
                            base_mu_rep, 
                            mu_beliefs_effect, 
                            beta, 
                            dist_beta, 
                            beta_control) {
    beliefs_latent = calculate_belief_latent_predictor(beta = beta, dist_beta = dist_beta, dist = dist)
    mu_rep = base_mu_rep * exp(mu_beliefs_effect * (beliefs_latent - beta_control))
    return(mu_rep)
}



#' Find Cutoff Type
#' 
#'  Given net benefit b, visibility \mu, total error sd and u_sd  we calculate
#' the cutoff type.
#' 
find_v_star = function(distance, b, mu_rep, total_error_sd, u_sd){
    v_f = generate_v_cutoff_fixedpoint(
        b = b,
        mu = mu_rep,
        total_error_sd = total_error_sd,
        u_sd = u_sd
    )

    soln = nleqslv(x = -b, fn = v_f)
    v_star = soln$x
    delta_v_star = calculate_delta(v_star, total_error_sd, u_sd)
    return(lst(
        v_star,
        delta_v_star
    ))
}


#' Find Predicted Takeup
#'
#' @param beta_b_z The coefficients on private benefit - treatment effect. 
#' @param beta_b_d The coefficients on private benefit - distance effect.
#' @param mu_beta_z The coefficients on beliefs that feed into \mu(z,d) model - treatment effect
#' @param mu_beta_d The coefficients on beliefs that feed into \mu(z,d) model - distance effect
#' @param total_error_sd Total sd of W
#' @param u_sd Variance of idiosyncratic shock.
#' @param dist_sd standard deviation of distance in the study. Stan uses 
#'  standardised distance to estimate models so we need to give it standardised 
#'  distance for our prediction exercise.
#' @param mu_beta_z_control The coefficient on mu_beta in the control arm - needed to 
#'  renormalise \mu(z,d) given \mu_0 exp(.) setup.
.find_pred_takeup = function(beta_b_z, 
                            beta_b_d,
                            mu_beta_z, 
                            mu_beta_d, 
                            base_mu_rep, 
                            total_error_sd, 
                            u_sd, 
                            dist_sd, 
                            mu_beta_z_control) {
    function(distance){
        distance = distance/dist_sd
        b = beta_b_z - beta_b_d*distance
        mu_rep = calculate_mu_rep(
            dist = distance,
            base_mu_rep = base_mu_rep,
            mu_beliefs_effect = 1,
            beta = mu_beta_z,
            dist_beta = mu_beta_d,
            beta_control = mu_beta_z_control)
        v_star_soln = find_v_star(
            distance = distance,
            b = b,
            mu_rep = mu_rep,
            total_error_sd = total_error_sd,
            u_sd = u_sd
        )

        delta_v_star = v_star_soln$delta_v_star
        v_star = v_star_soln$v_star

        linear_pred = b + mu_rep*delta_v_star
        pred_takeup = 1 - pnorm(v_star/(total_error_sd))
        return(lst(pred_takeup, linear_pred, b, mu_rep, delta_v_star, v_star = v_star_soln$v_star))
    }
}

find_pred_takeup = function(params) {
    .find_pred_takeup(
        beta_b_z = params$beta,
        beta_b_d = params$dist_beta_v,
        mu_beta_z = params$centered_cluster_beta_1ord,
        mu_beta_d = params$centered_cluster_dist_beta_1ord,
        base_mu_rep = params$base_mu_rep,
        total_error_sd = params$total_error_sd,
        u_sd = params$u_sd,
        dist_sd = params$dist_sd,
        mu_beta_z_control = params$mu_beta_z_control
    )
}


extract_params = function(param_draws, treatment, draw_id, j_id, dist_sd) {

    treatments = c(
        "control",
        "ink",
        "calendar",
        "bracelet"
    )
    treatment_id  = which(treatments == treatment)
    draw_df = param_draws %>%
        filter(.draw == draw_id) 

    params = draw_df %>%
        filter((k == treatment_id | is.na(k)) & (j == j_id | is.na(j)) ) %>%
        pull(.value, name = .variable)
    
    mu_beta_z_control = draw_df %>%
        filter(.variable == "centered_cluster_beta_1ord" & k == 1 & (j == j_id | is.na(j)) ) %>%
        pull(.value)
    params = c(params, "mu_beta_z_control" = mu_beta_z_control, "dist_sd" = dist_sd) %>%
        as.list()
    return(params)
}


draw_treat_grid =  expand.grid(
    draw = 1:800,
    treatment = treatments
)

pred_functions = map2(
    draw_treat_grid$draw,
    draw_treat_grid$treatment,
    ~extract_params(
        param_draws = param_draws,
        treatment = .y,
        draw_id = .x,
        dist_sd = dist_sd,
        j_id = 1
    ) %>% find_pred_takeup()
)

library(furrr)
plan(multisession, workers = 12)

pred_df = future_imap_dfr(pred_functions,
    ~tibble(
        distance = seq(from = 0, to = 5000, length.out = 40)) %>% 
            mutate(
                as_tibble(.x(distance)), 
                treatment = draw_treat_grid[.y, "treatment"],
                draw = draw_treat_grid[.y, "draw"]),
    .options = furrr_options(
        globals = c(
            "draw_treat_grid",
            "extract_params", 
            "param_draws", 
            "dist_sd", 
            "find_pred_takeup",
            ".find_pred_takeup",
            "calculate_mu_rep",
            "calculate_belief_latent_predictor",
            "find_v_star",
            "generate_v_cutoff_fixedpoint",
            "calculate_delta",
            "integrate_delta_part"),
        packages = c("dplyr", "nleqslv")
    ),
    .progress = TRUE
)  

treat_levels = c(
    "bracelet",
    "ink",
    "calendar",
    "control"
)
pred_df = pred_df %>%
    mutate(treatment = factor(treatment, levels = treat_levels))

pred_df %>%
    write_csv(
        file.path(
            script_options$output_path, 
            str_interp("pred_posterior_draws_dist_fit${fit_version}_STRUCTURAL_LINEAR_U_SHOCKS.csv")
        )
    )

summ_pred_df = pred_df %>%
    group_by(distance, treatment) %>%
    summarise_all(median) %>%
    mutate(treatment = factor(treatment, levels = treat_levels))

pred_df %>%
    ggplot(aes(
        x = distance, 
        y = pred_takeup, 
        group = interaction(draw, treatment), 
        colour = treatment
    )) +
    geom_line(
        alpha = 0.1
    )   +
    geom_line(
        data = summ_pred_df, 
        inherit.aes = FALSE,
        aes(
            x = distance, 
            y = pred_takeup, 
            group = interaction(treatment), 
            colour = treatment), 
        size = 1
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        title = "Predicted Takeup", 
        y = "Pred Takeup", 
        x = "Distance (m)"
    )

ggsave(
    "temp-plots/pred-takeup.png", 
    width = 10,
    height = 10,
    dpi = 500
    )

pred_df %>%
    ggplot(aes(
        x = distance, 
        y = pred_takeup, 
        group = interaction(draw, treatment),
        colour = treatment
    )) +
    geom_line(
        alpha = 0.3
    )  +
    geom_line(
        data = summ_pred_df, 
        inherit.aes = FALSE,
        aes(
            x = distance, 
            y = pred_takeup, 
            group = interaction(treatment)), 
        colour = "black",
        size = 1
    ) +
    facet_wrap(~treatment) +
    theme_bw() +
    guides(color = "none") +
    labs(
        title = "Pred Takeup", 
        y = "Takeup", 
        x = "Distance (m)"
    )
    
ggsave(
    "temp-plots/pred-takeup-facet.png", 
    width = 10,
    height = 10,
    dpi = 500
    )

    





filter_structural_fit = function(data){
    data %>%
        filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural")) 
}

treat_levels = c(
    "bracelet",
    "ink",
    "calendar", 
    "control"
)

dist_fit_data %>%
    filter_structural_fit() %>%
    select(cluster_takeup_prop) %>%
    unnest() %>%
    mutate(
        treatment = factor(assigned_treatment, levels = treat_levels), 
        type = "stan", 
        pred_takeup = mean_est, 
        distance = roc_distance
    ) 

dist_fit_data %>%
    filter_structural_fit() %>%
    select(cluster_takeup_prop) %>%
    unnest() %>%
    mutate(
        treatment = factor(assigned_treatment, levels = treat_levels), 
        type = "stan", 
        pred_takeup = per_0.5, 
        distance = roc_distance
    ) %>%
    bind_rows(
        summ_pred_df %>%
            filter(distance < 2500) %>%
            mutate(type = "extracted")
    ) %>%
    select( 
        treatment,
        distance, 
        pred_takeup, 
        type
    ) %>%
    ggplot(aes(
        x = distance, 
        y = pred_takeup, 
        colour = type 
    )) + 
    geom_point() +
    geom_line() +
    facet_wrap(~treatment) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(
        title = "Comparing Stan vs R Posterior Medians"
    )

ggsave(
    "temp-plots/temp-p-comp.png",
    width = 10,
    height = 10,
    dpi = 500
)
    
    
    
    
pred_df    
    
summ_pred_df  %>%
    mutate(
        r_by_d = mu_rep*delta_v_star
    ) %>%
    ggplot(aes(
        x = distance, 
        colour = treatment, 
        y = pred_takeup
    )) +
    geom_line() +
    geom_point()

cowplot::plot_grid(
    dist_fit_data %>%
        filter_structural_fit() %>%
        select(cluster_takeup_prop) %>%
        unnest() %>%
        mutate(
            assigned_treatment = factor(assigned_treatment, levels = treat_levels)
        ) %>%
        ggplot(aes(
            x = roc_distance, 
            y = per_0.5, 
            colour = assigned_treatment
        )) +
        geom_point() +
        geom_line(),

    summ_pred_df  %>%
        filter(distance <= 2500) %>%
        mutate(
            r_by_d = mu_rep*delta_v_star
        ) %>%
        ggplot(aes(
            x = distance, 
            colour = treatment, 
            y = pred_takeup
        )) +
        geom_line() +
        geom_point() 
) 

# cowplot::plot_grid(
#     dist_fit_data %>%
#         filter_structural_fit() %>%
#         select(cluster_rep_return) %>%
#         unnest() %>%
#         mutate(
#             assigned_treatment = factor(assigned_treatment, levels = treat_levels)
#         ) %>%
#         ggplot(aes(
#             x = roc_distance, 
#             y = per_0.5, 
#             colour = assigned_treatment
#         )) +
#         geom_point() +
#         geom_line(),
#     summ_pred_df  %>%
#         filter(distance <= 2500) %>%
#         mutate(
#             r_by_d = mu_rep*delta_v_star
#         ) %>%
#         ggplot(aes(
#             x = distance, 
#             colour = treatment, 
#             y = r_by_d
#         )) +
#         geom_line() +
#         geom_point()
# )

pred_df  %>%
    # filter(v_star > -5) %>%
    ggplot(aes(
        x = v_star, 
        group = interaction(draw, treatment),
        colour = treatment, 
        y = delta_v_star
    )) +
    geom_line(alpha = 0.2) +
    geom_point(alpha = 0.2) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        title = "Delta(v^*) vs v^*"
    )

ggsave(
    "temp-plots/all-delta-v-vs-v.png",
    width = 10,
    height = 10,
    dpi = 500
)



pred_df  %>%
    filter(abs(v_star) < 2) %>%
    # filter(v_star > -5) %>%
    ggplot(aes(
        x = v_star, 
        group = interaction(draw, treatment),
        colour = treatment, 
        y = delta_v_star
    )) +
    geom_line(alpha = 0.2) +
    geom_point(alpha = 0.2) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        title = "Delta(v^*) vs v^*"
    )

ggsave(
    "temp-plots/subset-delta-v-vs-v.png",
    width = 10,
    height = 10,
    dpi = 500
)

summ_pred_df  %>%
    # filter(v_star > -5) %>%
    ggplot(aes(
        x = v_star, 
        colour = treatment, 
        y = delta_v_star
    )) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        title = "Delta(v^*) vs v^*"
    )


ggsave(
    "temp-plots/delta-v-vs-v.png",
    width = 10,
    height = 10,
    dpi = 500
)


summ_pred_df  %>%
    ggplot(aes(
        x = distance, 
        colour = treatment, 
        y = v_star
    )) +
    geom_line() +
    geom_point()


dist_fit_data %>%
    filter_structural_fit()  %>%
    colnames()

dist_fit_data %>%
    filter_structural_fit()  %>%
    select(cluster_w_cutoff) %>%
    unnest() %>%
    ggplot(aes(
        x = roc_distance, 
        y = mean_est, 
        colour = assigned_treatment
    )) +
    geom_point() +
    geom_line()


pred_df  %>%
    filter(v_star > -5) %>%
    ggplot(aes(
        x = distance, 
        colour = treatment, 
        y = delta_v_star
    )) +
    geom_line() +
    geom_point()


pred_df  %>%
    filter(v_star > -5) %>%
    filter(distance <= 2.5) %>%
    # filter(treatment %in% c("bracelet", "control")) %>%
    ggplot(aes(
        x = distance, 
        colour = treatment, 
        y = pred_takeup
    )) +
    geom_line() +
    # geom_point() +
    geom_hline( 
        yintercept = c(0.5, 0.2), 
        linetype = "longdash"
    ) 


dist_prob_plot <- dist_fit_data %>% 
    filter(fct_match(fit_type, "fit"), fct_match(model_type, "structural"))  %>%
    select(cluster_takeup_prop) %>%
    unnest() %>%
    ggplot(aes(
        x = roc_distance, 
        y = mean_est, 
        colour = assigned_treatment
    )) +
    geom_line()


dist_prob_plot  



pred_df  %>%
    filter(v_star > -5) %>%
    ggplot(aes(
        x = distance, 
        colour = treatment, 
        y = v_star
    )) +
    geom_line() +
    geom_point()

pred_df %>%
    filter(distance == max(distance))



pred_df %>%
    ggplot(aes(
        x = distance,
        y = y_hat,
        colour = treatment
    )) +
    geom_line()


analysis_data %>%
    select(standard_cluster.dist.to.pot) %>%
    summary()



ed = find_pred_takeup(
    distance = 0.5, 
    beta = betas[[1]], 
    mu_beta = mu_betas[[1]], 
    base_mu_rep = 0.121, 
    total_error_sd = 1.2,
    u_sd =0)

pred_df = tibble(
    dist = seq(from = -2, to = 2, length.out = 100)
) %>%
    mutate(y_hat = ed(dist))


pred_df %>%
    ggplot(aes(
        x = dist,
        y = y_hat
    )) +
    geom_point()



find_v_star(
    distance = 0.5, 
    beta = betas[[1]], 
    mu_beta = mu_betas[[1]], 
    base_mu_rep = 0.121, 
    total_error_sd = 1.2,
    u_sd =0)



ed(1)





param_summary %>%
    head(11)



variables = fitted_draw_df %>%
    colnames()

variables[!str_detect(variables, "cf")]

model_fit$summary("beta_dist_v")

variables(model_fit)
model_fit$variables()


model_fit$metadata() %>% head(20)

str(model_fit)


dist_df = imap_dfr(
    mu_betas,
    ~expand.grid(
        dist = seq(from = 0, to = 4, length.out = 100)
    ) %>%
    tibble() %>%
    mutate(
        mu_rep = calculate_mu_rep(dist = dist, 
                                  base_mu_rep = clean_param_summary %>%  
                                    filter(variable == "base_mu_rep") %>%
                                    pull(mean), 
                                  mu_beliefs_effect = 1, 
                                  beta = .x[1], 
                                  dist_beta = .x[2], 
                                  beta_control = mu_betas[[1]][1]),
        treatment = treatments[.y]
    )
) %>%
    mutate(orig_dist = dist*dist_sd)


dist_df %>%
    # filter(treatment %in% c("bracelet", "control")) %>%
    ggplot(aes(
        x = orig_dist, 
        y = mu_rep, 
        colour = treatment
    )) +
    geom_point() +
    labs(title = "Mu(z, d) from model") +
    theme_bw()


