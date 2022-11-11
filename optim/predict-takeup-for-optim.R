library(magrittr)
library(posterior)
library(tidyverse)
library(broom)
library(ggrepel)
library(ggmap)
library(ggstance)
library(gridExtra)
library(cowplot)
library(rgeos)
library(sp)
library(knitr)
library(modelr)
library(car)
library(rstan)
library(latex2exp)
library(ggthemes)

library(nleqslv)
library(cmdstanr)
library(econometr)
source(file.path("rct-design-fieldwork", "takeup_rct_assign_clusters.R"))
source(file.path("analysis_util.R"))
source(file.path( "dist_structural_util.R"))


source(file.path("multilvlr", "multilvlr_util.R"))

fit_version <- 66

# 66 ed fit
# 60 Karim fit
# 62 also Karim fit


model_fit_by = if_else(fit_version %in% c(60, 62), "Karim", "Ed")


models_we_want = c(
  "STRUCTURAL_LINEAR_U_SHOCKS"
)


quant_probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

output_basepath = str_glue("temp-data/output_dist_fit{fit_version}")

canva_palette_vibrant <- "Primary colors with a vibrant twist"

theme_set(theme_minimal() +
            theme(legend.position = "bottom"))

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

load(file.path("temp-data", str_interp("processed_dist_fit${fit_version}_lite.RData")))

dist_fit_data %<>%
  left_join(tribble(
    ~ fit_type,        ~ model_color,
      "fit",           "black", 
      "prior-predict", "darkgrey",
  ), by = "fit_type")


stan_data = (dist_fit_data %>%
    slice(1) %>%
    pull(stan_data))[[1]]

rf_analysis_data = stan_data$analysis_data

model_files = fs::dir_ls("data/stan_analysis_data", regex = str_glue("dist_fit{fit_version}_STRUCTURAL_LINEAR_U_SHOCKS.*\\.csv"))
model_fit = as_cmdstan_fit(model_files)

stop()
#  Levels: control ink calendar bracelet
param_summary = model_fit$summary(c(
    "beta", 
    "dist_beta_v",
    "centered_cluster_beta_1ord",
    "centered_cluster_dist_beta_1ord",
    "base_mu_rep", 
    "total_error_sd",
    "u_sd"))
#  saveRDS(param_summary, "temp-data/ed-temp-param-summary.rds")

# param_summary = read_rds("temp-data/ed-temp-param-summary.rds")

param_summary


#  rm(model_fit)
#  gc()

## B(z,d):
# \beta is treatment effect
# dist_beta_v is distance cost


clean_param_summary = param_summary %>%
    mutate(
        index_1 = as.numeric(str_extract(variable, "(?<=\\[)\\d+(?=\\,|\\])")),
        index_2 = as.numeric(str_extract(variable, "(?<=,)\\d+(?=\\])"))) %>%
    mutate(variable = str_remove_all(variable, "\\[.*$"))   %>%
    mutate(index_1 = if_else(str_detect(variable, "centered"), index_2, index_1)) %>% # should have used tidybayes
    group_by(variable, index_1) %>%
    select(mean, median) %>%
    unique()




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

# matrix[num_clusters, num_treatments] centered_cluster_dist_beta_1ord; 



# very hardcoded
calculate_belief_latent_predictor = function(beta, dist_beta, dist) {
    # intercept + treatment - intercept
    val =  beta + dist*dist_beta 
    return(val)
}

calculate_mu_rep = function(dist, base_mu_rep, mu_beliefs_effect, beta, dist_beta, beta_control) {

    beliefs_latent = calculate_belief_latent_predictor(beta = beta, dist_beta = dist_beta, dist = dist)
    return(
        base_mu_rep * exp(mu_beliefs_effect * (beliefs_latent - beta_control))
    )
}

extract_mu_betas = function(summ_df, index){
    summ_df %>%
        filter(index_1 == index) %>%
        filter(str_detect(variable, "cluster")) %>%
        pull(mean, name = variable)
}

extract_betas = function(summ_df, index){
    summ_df %>%
        filter((variable == "beta" & index_1 == index) | (variable == "dist_beta_v"))  %>%
        pull(mean, name = variable)
}

mu_betas = map(1:4, ~extract_mu_betas(clean_param_summary, .x))
betas = map(1:4, ~extract_betas(clean_param_summary, .x))
# mu_betas[[1]][1] = 0 # control have to net off intercept and use \mu_0 pre-mult exp


betas



treatments = c(
    "control",
    "ink", 
    "calendar", 
    "bracelet"
)


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


find_pred_takeup = function(beta, mu_beta, base_mu_rep, total_error_sd, u_sd, dist_sd, beta_control) {
    function(distance){
        distance = distance/dist_sd
        b = beta[1] - beta[2]*distance
        mu_rep = calculate_mu_rep(
            dist = distance,
            base_mu_rep = base_mu_rep,
            mu_beliefs_effect = 1,
            beta = mu_beta[1],
            dist_beta = mu_beta[2],
            beta_control = beta_control)
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

clean_param_summary
hat_total_error_sd = clean_param_summary %>%
    ungroup() %>%
    filter(variable == "total_error_sd") %>%
    select(mean) %>%
    unique() %>%
    pull()
hat_u_sd = clean_param_summary %>%
    ungroup() %>%
    filter(variable == "u_sd") %>%
    select(mean) %>%
    unique() %>%
    pull()
hat_base_mu_rep = clean_param_summary %>%
    ungroup() %>%
    filter(variable == "base_mu_rep") %>%
    select(mean) %>%
    unique() %>%
    pull()
hat_dist_beta_v = clean_param_summary %>%
    ungroup() %>%
    filter(variable == "dist_beta_v") %>%
    select(mean) %>%
    unique() %>%
    pull()
hat_base_mu_rep
hat_total_error_sd
hat_u_sd


pred_funcs = map(1:4,
            ~find_pred_takeup(
                beta = betas[[.x]],
                mu_beta = mu_betas[[.x]],
                base_mu_rep = hat_base_mu_rep,
                total_error_sd = hat_total_error_sd,
                u_sd = hat_u_sd,
                dist_sd = dist_sd,
                beta_control = mu_betas[[1]][1]
            )
)


treat_levels = c(
    "bracelet",
    "ink",
    "calendar",
    "control"
)
pred_df = imap_dfr(pred_funcs,
    ~tibble(distance = seq(from = 0, to = 2500, length.out = 40)) %>% mutate(as_tibble(.x(distance)), treatment = treatments[.y])
)  %>%
    mutate(
        treatment = factor(treatment, levels = treat_levels)
    ) 


pred_df  %>%
    ggplot(aes(
        x = distance, 
        colour = treatment, 
        y = mu_rep
    )) +
    geom_line() +
    geom_point()

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

pred_df

pred_df  %>%
    mutate(
        r_by_d = mu_rep*delta_v_star/hat_dist_beta_v
    ) %>%
    ggplot(aes(
        x = distance, 
        colour = treatment, 
        y = pred_takeup
    )) +
    geom_line() +
    geom_point()

dist_fit_data %>%
    filter_structural_fit() %>%
    select(cluster_takeup_prop) %>%
    unnest() %>%
    mutate(
        treatment = factor(assigned_treatment, levels = treat_levels), 
        type = "stan", 
        pred_takeup = mean_est, 
        distance = roc_distance
    ) %>%
    bind_rows(
        pred_df %>%
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
        title = "Comparing Stan Posterior Means and R Extractions"
    )
ggsave(
    "temp-p-comp.png",
    width = 10,
    height = 10,
    dpi = 500
)
    
    
    
    
    
    
    %>%
    ggplot(aes(
        x = roc_distance, 
        y = mean_est, 
        colour = assigned_treatment
    )) +
    geom_point() +
    geom_line(),

pred_df  %>%
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
            y = mean_est, 
            colour = assigned_treatment
        )) +
        geom_point() +
        geom_line(),

    pred_df  %>%
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
cowplot::plot_grid(
    dist_fit_data %>%
        filter_structural_fit() %>%
        select(cluster_rep_return) %>%
        unnest() %>%
        mutate(
            assigned_treatment = factor(assigned_treatment, levels = treat_levels)
        ) %>%
        ggplot(aes(
            x = roc_distance, 
            y = mean_est, 
            colour = assigned_treatment
        )) +
        geom_point() +
        geom_line(),

    pred_df  %>%
        mutate(
            r_by_d = mu_rep*delta_v_star
        ) %>%
        ggplot(aes(
            x = distance, 
            colour = treatment, 
            y = r_by_d
        )) +
        geom_line() +
        geom_point()
)

pred_df  %>%
    filter(v_star > -5) %>%
    ggplot(aes(
        x = v_star, 
        colour = treatment, 
        y = delta_v_star
    )) +
    geom_line() +
    geom_point()

pred_df  %>%
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
