
#!/usr/bin/Rscript

script_options = docopt::docopt(
    stringr::str_glue("Usage:
    predict-takeup-for-optim.R <fit-version> [options] [ --from-csv | --to-csv ] [<treatment>...]

    Takes posterior fits and calculates estimated takeup for a given distance.


    Options:
        --from-csv  Load parameter draws from a csv file
        --to-csv   Load stanfit object and write draws to csv
        --input-path=<path>  Path to find results [default: {file.path('data', 'stan_analysis_data')}]
        --output-path=<path>  Path to find results [default: {file.path('optim', 'data')}]
        --output-name=<output-name>  Prepended to output file
        --pred-distance  Estimate takeup on a grid from 0 to 5km
        --dist-cutoff=<dist-cutoff>  Don't estimate takeup past this cutoff (set takeup to 0) [default: 2500]
        --rep-cutoff=<rep-cutoff>  Don't let reputational returns increase past this point. [default: 2500]
        --num-post-draws=<num-post-draws>  Number of posterior draws to use [default: 1600]
        --num-cores=<num-cores>  Number of cores to use [default: 8]
    "),
    args = if (interactive()) "
                            66 
                            --output-name=rep
                            --from-csv
                            --num-post-draws=80
                            --rep-cutoff=2500
                            --dist-cutoff=5000
                            --num-cores=12
                            bracelet
                            control
                              " 
           else commandArgs(trailingOnly = TRUE)
)
if (length(script_options$treatment) == 0) {
    script_options$treatments = c(
        "control",
        "ink", 
        "calendar", 
        "bracelet"
    )
}


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
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))


source(file.path("multilvlr", "multilvlr_util.R"))

script_options = script_options %>%
    modify_at(c("dist_cutoff", 
                "rep_cutoff",
                "num_post_draws", 
                "num_cores"), as.numeric)
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

sd_of_dist = sd(rf_analysis_data$cluster.dist.to.pot)

if (script_options$to_csv) {
    struct_model_files = fs::dir_ls(script_options$input_path, regex = str_glue("dist_fit{fit_version}_STRUCTURAL_LINEAR_U_SHOCKS.*\\.csv"))
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
        write_csv(file.path(script_options$input_path, str_interp("param_posterior_draws_dist_fit${fit_version}_STRUCTURAL_LINEAR_U_SHOCKS.csv") ))
    ## Now reduced form  - we actually fit rf model continuous in distance here
    # since main analysis RF just uses close/far

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
        dewormed ~ (1 | cluster.id) + sd_dist*assigned.treatment, 
        bernoulli(link = "probit")
    )

    rf_model_fit %>%
        saveRDS(
            file.path(
                script_options$input_path, 
                str_interp("param_posterior_draws_dist_fit${fit_version}_REDUCED_FORM_NO_RESTRICT.rds") ))
    rm(rf_model_fit)
    gc() # not as big but get rid of it anyway

} 

if (script_options$from_csv) {
    struct_param_draws = read_csv(
        file.path(
            script_options$input_path, 
            str_interp(
                "param_posterior_draws_dist_fit${fit_version}_STRUCTURAL_LINEAR_U_SHOCKS.csv"
            )
        )
    )
    rf_model_fit = read_rds(
        file.path(
            script_options$input_path, 
            str_interp(
                "param_posterior_draws_dist_fit${fit_version}_REDUCED_FORM_NO_RESTRICT.rds"
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
    dim_distance = length(distance)
    n_iters = 5

    v_f = generate_v_cutoff_fixedpoint(
        b = b,
        mu = mu_rep,
        total_error_sd = total_error_sd,
        u_sd = u_sd
    )
    starting_points = map(1:n_iters, ~rnorm(dim_distance))
    solns = map(starting_points, ~nleqslv(x = 0 + .x, fn = v_f))
    soln_matrix = matrix(
        map(solns, "x") %>% 
            unlist(), 
        nrow = n_iters, 
        ncol = dim_distance,
        byrow = TRUE )
    # Take median from nleqslv with 10 different starting points
    v_star = apply(soln_matrix, 2, median)
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
                             mu_beta_z_control,
                             rep_cutoff = Inf) {
    function(distance){
        over_cutoff = distance > rep_cutoff # note rep_cutoff not standardised
        distance = distance/dist_sd
        b = beta_b_z - beta_b_d*distance
        mu_rep = calculate_mu_rep(
            dist = distance,
            base_mu_rep = base_mu_rep,
            mu_beliefs_effect = 1,
            beta = mu_beta_z,
            dist_beta = mu_beta_d,
            beta_control = mu_beta_z_control)
        # if distance greater than cutoff, set mu_rep to cutoff mu_rep within distance
        cutoff_mu_rep = calculate_mu_rep(
            dist = rep_cutoff/dist_sd,
            base_mu_rep = base_mu_rep,
            mu_beliefs_effect = 1,
            beta = mu_beta_z,
            dist_beta = mu_beta_d,
            beta_control = mu_beta_z_control)
        mu_rep[which(over_cutoff)] = cutoff_mu_rep

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
        mu_beta_z_control = params$mu_beta_z_control,
        rep_cutoff = params$rep_cutoff
    )
}


extract_params = function(param_draws, 
                          treatment, 
                          draw_id, 
                          j_id, 
                          dist_sd,
                          dist_cutoff = Inf,
                          rep_cutoff = Inf) {

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
    

    params = c(
        params, 
        "mu_beta_z_control" = mu_beta_z_control, 
        "dist_sd" = dist_sd,
        "dist_cutoff" = dist_cutoff,
        "rep_cutoff" = rep_cutoff
        ) %>%
        as.list()

    return(params)
}

## Create estimated demand functions
max_draw = max(struct_param_draws$.draw)
draw_treat_grid = expand.grid(
    draw = 1:min(max_draw, script_options$num_post_draws),
    treatment = script_options$treatment
)


pred_functions = map2(
    draw_treat_grid$draw,
    draw_treat_grid$treatment,
    ~extract_params(
        param_draws = struct_param_draws,
        treatment = .y,
        draw_id = .x,
        dist_sd = dist_sd,
        j_id = 1,
        rep_cutoff = 2500,
        dist_cutoff = Inf
    ) %>% find_pred_takeup()
)
## Now extract distance data
## Grabbing Data
cluster_ids_actually_used = analysis_data %>%
  select(cluster.id) %>%
  pull() %>%
  unique()

rct_cluster_df = st_as_sf(rct.cluster.selection) 

rct_pot_cluster_ids = rct_cluster_df %>%
  filter(selected == TRUE) %>%
  filter(cluster.id %in% cluster_ids_actually_used) %>%
  pull(cluster.id)

n_orig_pots = length(rct_pot_cluster_ids)

rct_school_df = st_as_sf(rct.schools.data)

unselected_pots = setdiff(
  rct_school_df$cluster.id,
  rct_pot_cluster_ids)

# additional_pots = sample(unselected_pots, round(n_orig_pots*1.1), replace = FALSE)
additional_pots = NULL



pot_df = rct_school_df %>%
  filter(cluster.id %in% c(additional_pots, rct_pot_cluster_ids)) %>%
  mutate(id = 1:n())


village_df = village.centers %>%
  st_as_sf(coords = c("lon", "lat"), crs = wgs.84) %>%
  mutate(id = 1:n())



distance_matrix = st_distance(
  village_df,
  pot_df
)

long_distance_mat = distance_matrix %>%
    as_tibble() %>%
    mutate(index_i = 1:n()) %>%
    gather(variable, dist, -index_i) %>%
    mutate(index_j = str_extract(variable, "\\d+") %>% as.numeric())  %>%
    select(index_i, index_j, dist) %>%
    mutate(
        dist = as.numeric(dist),
        dist_km = as.numeric(dist/1000))


brms_long_distance_mat = long_distance_mat %>%
    left_join(
        village_df %>% as_tibble() %>% select(cluster.id, id),
        by = c("index_i" = "id")
    )  %>%
    mutate(
        sd_dist = dist / sd_of_dist, 
        assigned.treatment = "bracelet"
    ) 


rf_pred_df = map_dfr(
    script_options$treatment,
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

tictoc::tic()
plan(multisession, workers = script_options$num_cores)
pred_df = future_imap_dfr(
    pred_functions,
    ~{
        long_distance_mat %>%
        filter(dist < script_options$dist_cutoff) %>%
            mutate( 
                as_tibble(.x(dist)),
                treatment = draw_treat_grid[.y, "treatment"], 
                draw = draw_treat_grid[.y, "draw"]
            )
    },
    .options = furrr_options(
        globals = c(
            "draw_treat_grid",
            "script_options",
            "long_distance_mat",
            "extract_params", 
            "struct_param_draws", 
            "dist_sd", 
            "find_pred_takeup",
            ".find_pred_takeup",
            "calculate_mu_rep",
            "calculate_belief_latent_predictor",
            "find_v_star",
            "generate_v_cutoff_fixedpoint",
            "calculate_delta",
            "integrate_delta_part"),
        packages = c("dplyr", "nleqslv", "purrr"),
        seed = TRUE
    ),
    .progress = TRUE
)
tictoc::toc()


cutoff_pred_df = future_map2_dfr(
    draw_treat_grid$draw,
    draw_treat_grid$treatment,
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
                treatment = .y, 
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
    model = "STRUCTURAL_LINEAR_U_SHOCK"
  )


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


append_output = if (!is.null(script_options$output_name)) {
    paste0("_", script_options$output_name)
} else {
    ""
}

bind_rows(
    structural_demand_df,
    rf_demand_df
) %>%
    write_csv(
        file.path(
            script_options$output_path, 
            str_interp("pred_demand_dist_fit${fit_version}{append_output}.csv")
        )
    )
        

pot_df %>%
    mutate(
        lon = sf::st_coordinates(.)[, 1],
        lat = sf::st_coordinates(.)[, 2]) %>%
    as_tibble() %>%
    write_csv(
        file.path(
            script_options$output_path,
            "pot-df.csv"
        ))

village_df %>%
    mutate(
        lon = sf::st_coordinates(.)[, 1],
        lat = sf::st_coordinates(.)[, 2]) %>%
    as_tibble() %>%
    write_csv(
        file.path(
            script_options$output_path,
            "village-df.csv"
        ))



if (script_options$pred_distance) {


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
                "struct_param_draws", 
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


}
