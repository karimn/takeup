library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(posterior)


fit_file = "data/stan_analysis_data/dist_fit87_STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-1.csv"
fit = as_cmdstan_fit(fit_file)


# # To find fixed point need:
# # benefit_cost
# # mu_rep
# # total_error_sd
# # u_sd


bc_draws = spread_rvars(
    fit,
    # cluster_rep_return[i,j,k],
    # cluster_rep_return_dist[i,j,k],
    # cluster_w_cutoff[i,j,k], 
    # cluster_w_control_cutoff[i,j,k], 
    # total_error_sd[1],
    # dist_beta_v[k]
    cluster_roc[i,j,k], 
    cluster_roc_no_vis[i,j,k]
    )

rm(fit)
gc()

summ_bc_draws = bc_draws %>%
    filter(j == 1) %>%
    pivot_longer(
        cols = c(cluster_roc, cluster_roc_no_vis)
    ) 

summ_bc_draws = summ_bc_draws %>%
    median_qi(value, .width = c(0.80, 0.50))

dist_df = tibble(
    i = 1:51,
    dist = seq(from = 0, to = 5000, by = 100)
)
summ_bc_draws = summ_bc_draws %>%
    left_join(dist_df)

treat_df = tibble(
    k = 1:4, 
    treatment = c("control", "ink", "calendar", "bracelet")
)
summ_bc_draws = summ_bc_draws %>%
    left_join(treat_df)

summ_bc_draws = summ_bc_draws %>%
    mutate(
        treatment = factor(treatment, levels = c("control", "ink", "calendar", "bracelet"))
    )

sd_of_dist = read_rds("temp-data/sd_of_dist.rds")
plot_summ_bc_draws = summ_bc_draws %>%
    filter(dist < 2500) %>%
    mutate(dist = dist/1000) %>%
    mutate(
        across(c(value, .lower, .upper), magrittr::divide_by, sd_of_dist), 
        across(c(value, .lower, .upper), magrittr::multiply_by, 1000),
        across(c(value, .lower, .upper), magrittr::multiply_by, 100))


load(file.path("temp-data", str_interp("processed_dist_fit86_lite.RData")))
library(magrittr)
rf_analysis_data <- dist_fit_data %>% 
  filter(
    fct_match(model_type, "reduced form"),
    fct_match(fit_type, "fit"),
  ) %$% 
  stan_data[[1]]$analysis_data 

plot_summ_bc_draws %>%
    filter(name == "cluster_roc") %>%
    filter(treatment %in% c("control", "bracelet")) %>%
    mutate(k = factor(k)) %>%
        ggplot(aes(
            x = dist
        )) +
        geom_line(aes(
            color = treatment,
            y = value
        )) +
        geom_ribbon(
            data = . %>%
                filter(.width == 0.5),
            aes(
                ymax = .upper, 
                ymin = .lower,
                fill = treatment), alpha = 0.3) +
        geom_ribbon(
            data = . %>%
                filter(.width == 0.8),
            aes(
                ymax = .upper, 
                ymin = .lower,
                fill = treatment), alpha = 0.3) +
        theme_minimal() +
        theme( 
            legend.position = "bottom",
            legend.title = element_blank()
        )  +
        labs(
            x = "Distance to Treatment (d) [km]", 
            # y = "Estimated Takeup", 
            colour = ""
        ) +
        guides(
            linetype = "none"
        ) +
        labs(
            caption = "Line: Median. Outer ribbon: 80% credible interval. Inner ribbon: 50% credible interval." 
        ) +
        geom_line(
            inherit.aes = FALSE, 
            aes(x = dist, y = value),
            data = plot_summ_bc_draws %>%
                filter(name == "cluster_roc_no_vis"),
            linetype = "longdash"
        )  +
        geom_rug(
            inherit.aes = FALSE,
            aes(dist),
            alpha = 0.75,
            data = rf_analysis_data %>%
            filter(fct_match(assigned.treatment, c("control",  "bracelet"))) %>%
            distinct(cluster_id, assigned_treatment = assigned.treatment, dist = cluster.dist.to.pot / 1000)) +
        annotate(
            "text", 
            x = 0.3, 
            y = -12, 
            label = "Amplification", 
            size = 4, 
            alpha = 0.7
        ) +
        annotate(
            "text", 
            x = 2.3, 
            y = -7.5, 
            label = "Mitigation", 
            size = 4, 
            alpha = 0.7
        ) +
        labs(
            y = "Rate of Change [pp/km]"
        ) +
        ggthemes::scale_color_canva( 
            "",
            palette = "Primary colors with a vibrant twist", 
            labels = str_to_title
        ) +
        ggthemes::scale_fill_canva( 
            "",
            palette = "Primary colors with a vibrant twist", 
            labels = str_to_title
        ) 

ggsave("temp-plots/last-minute-roc-plot.pdf", width = 8, height = 6)
    
    # %>%




summ_bc_draws = summ_bc_draws %>%
    left_join(dist_df)


summ_bc_draws %>%
    filter(dist <= 2500) %>%
    ggplot(aes(
        x = dist, 
        y = cluster_roc, 
        colour = factor(k )
    )) +
    geom_point() +
    geom_line() +
    geom_line(
        aes(y = cluster_roc_no_vis), 
        colour = "black"
    )

bc_draws %>%
    select(cluster_w_control_cutoff, total_error_sd) 

r_dnorm = rfun(dnorm)

bc_draws = bc_draws %>%
    mutate(
        delta_roc =  -1*dist_beta_v*r_dnorm(cluster_w_control_cutoff, sd = total_error_sd)
    ) %>%
    left_join(
        dist_df
    )


summ_bc_draws = bc_draws %>%
    mutate(across(where(is_rvar), mean))

summ_bc_draws %>%
    filter(i == 1) %>%
   select(i, j, k, delta_roc) %>%
   tail()

summ_bc_draws %>%
    filter(dist < 2500) %>%
    ggplot(aes(
        x = dist, 
        y = delta_roc, 
    )) +
    geom_line() +
    geom_point()



    
summ_bc_draws = bc_draws %>%
    mutate(
        across( 
            where(is_rvar), mean
        )
    )

summ_bc_draws %>%
    ggplot(aes(
        x = 
    ))


# bc_draws = fit$draws(
#   variables = c(
#     "structural_cluster_benefit_cost", 
#     "cluster_dist_cost", 
#     "obs_cluster_mu_rep",
#     "total_error_sd[1]", 
#     "u_sd[1]", 
#     "beta", 
#     "dist_beta_v")
# )

# rm(fit)
# gc()

# control_mu = 114
# bracelet_mu = 93
# calendar_mu = 110

# bracelet_b = 106
# ink_mu = 94
# ink_b = 138

# pct_change = function(x, y){100*(x/y - 1)}


# pct_change(bracelet_mu, control_mu)
# pct_change(calendar_mu, control_mu)
# pct_change(bracelet_b, control_mu)
# pct_change(ink_mu, control_mu)
# pct_change(ink_b, control_mu)





# betas = spread_rvars(
#     bc_draws, 
#     beta[k], 
#     dist_beta_v[k]
#     )

# betas


# optim_files = fs::dir_ls(
#     "optim/data/",
#     regexp = "*optimal-allocation.rds"
# )


# pred_demand_df = read_csv(
#     "optim/data/all-treat-demand-dist-fit71-STRUCTURAL_LINEAR_U_SHOCKS.csv"
# )

# demand_files = fs::dir_ls(
#     "optim/data/agg-log-full",
#     regexp = "pred"
# )

# extract_and_harmonise_variables = function(stan_df){
#     v_star_df = stan_df %>%
#         select(cluster_w_cutoff) %>%
#         unnest(cols = c(cluster_w_cutoff))  %>%
#         select(roc_distance, assigned_treatment, mean_est, per_0.5) %>%
#         rename(
#             dist = roc_distance,
#             treatment = assigned_treatment, 
#             mean_est = mean_est, 
#             median_est = per_0.5
#         ) %>%
#         mutate(variable = "v_star")

#     pred_takeup_df = stan_df %>%
#         select(obs_cluster_takeup_level) %>%
#         unnest(cols = c(obs_cluster_takeup_level)) %>%
#         select(
#             dist = assigned_dist_obs, 
#             treatment = assigned_treatment_obs, 
#             median_est = per_0.5
#         ) %>% 
#         mutate(variable = "actual_takeup")

#     rep_return_df = stan_df %>%
#         select(cluster_rep_return) %>%
#         unnest(cols = c(cluster_rep_return)) %>%
#         select( 
#             dist = roc_distance, 
#             mean_est = mean_est, 
#             median_est = per_0.5, 
#             treatment = assigned_treatment
#         ) %>% 
#         mutate(variable = "rep_return")

#     prop_df = stan_df %>%
#         select(cluster_takeup_prop) %>%
#         unnest(c(cluster_takeup_prop)) %>%
#         select(
#             dist = roc_distance, 
#             mean_est = mean_est, 
#             median_est = per_0.5, 
#             treatment = assigned_treatment
#         ) %>%
#         mutate(variable = "pred_takeup")

#     df = bind_rows(
#         v_star_df, 
#         pred_takeup_df, 
#         rep_return_df, 
#         prop_df
#     )
#     return(df)

# }




# fit_version = 71
# load(file.path("temp-data", str_interp("processed_dist_fit${fit_version}_lite.rdata")))

# stan_fit_df = dist_fit_data %>%
#     filter(model == "structural_linear_u_shocks") %>%
#     filter(fct_match(fit_type, "fit"))

# stan_clean_df = extract_and_harmonise_variables(
#     stan_fit_df
# )


# pred_demand_df


# summ_df %>%
#     ggplot(aes(
#         x = dist, 
#         y = mean_est, 
#         colour = treatment
#     )) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~variable, scales = "free")
# stop()


# summ_df %>%
#     filter(variable == "pred_takeup") %>%
#     ggplot(aes(
#         x = dist, 
#         y = mean_est, 
#         colour = treatment
#     )) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~variable, scales = "free")


# comp_df = bind_rows(
#     summ_df %>% mutate(type = "ed"), 
#     stan_clean_df %>% mutate(type = 'stan')
# )

# stan_clean_df %>%
#     filter(variable == "pred_takeup")

# comp_df %>%
#     filter(variable == "pred_takeup") %>%
#     ggplot(aes(
#         x = dist, 
#         y = mean_est, 
#         colour = type
#     )) +
#     geom_point() +
#     facet_wrap(~treatment)

#     ggsave("temp-plots/uh-oh.png", width = 8, height = 6, dpi = 500)



# if (use_map_demand) {

# demand_files = fs::dir_ls(
#     "optim/data/agg-log-full",
#     regexp = "pred"
# )
# demand_files = demand_files[!str_detect(demand_files, "no-cutoff") & str_detect(demand_files, "STRUCTURAL_LINEAR_U_SHOCKS.csv")]



# demand_df = tibble(
#     demand_file = demand_files
# ) %>%
#     mutate(demand_df = map(demand_file, read_csv))



# demand_df = demand_df  %>%
#     mutate(
#         cutoff_type = if_else(
#             str_detect(demand_file, "no-cutoff"),
#             "no_cutoff", 
#             "cutoff"
#         ),
#         b_type = str_extract(demand_file, "(?<=b-)(.*?)(?=-)"), 
#         mu_type = str_extract(demand_file, "(?<=mu-)(.*?)(?=-)"), 
#         model_type = str_extract(demand_file, str_glue("(?<=mu-{mu_type}-)(.*?)(?=\\.csv)"))
#     )   


# subset_demand_df = demand_df %>%
#     select(-demand_file) %>%
#     mutate(
#         subset_df = map(
#             demand_df, 
#             ~{
#                 filter(.x, dist <= 10000) %>%
#                 filter(model != "REDUCED_FORM_NO_RESTRICT") %>%
#                     group_by(village_i, pot_j) %>%
#                     summarise(across(where(is.numeric), mean, na.rm = TRUE))
#             }
#         )
#     )  %>%
#     select(-demand_df)



# rm(demand_df)

# long_subset_demand_df = subset_demand_df %>%
#     unnest(subset_df)
# }


# long_subset_demand_df %>%
#     filter(b_type == mu_type) %>%
#     filter(cutoff_type == "cutoff") %>%
#     filter(dist < 2500) %>%

#     ggplot(aes(
#         x = dist, 
#         y = demand, 
#         colour = b_type
#     )) +
#     geom_point()

# long_subset_demand_df %>%
#     filter(b_type == mu_type) %>%
#     filter(cutoff_type == "cutoff") %>%
#     filter(dist < 2500) %>%
#     mutate(close = dist < 1250)  %>%
#     select( 
#         b_type, 
#         mu_type, 
#         dist, 
#         demand, 
#         close
#     )  %>%
#     group_by(
#         b_type, 
#         close
#     ) %>%
#     summarise(
#         mean_demand = mean(demand)
#     )



# plot_demand_panels = function(long_data, cutoff, model, b = "control"){
#     p = long_data %>%
#         filter(cutoff_type == cutoff) %>%
#         filter(model_type == model) %>%
#         filter(b_type == b)  %>%
#         filter(dist <= 5000) %>%
#         gather(variable, value, demand:v_star) %>%
#         mutate(variable = factor(variable,
#                                 levels = c("b", "v_star", "delta_v_star", "mu_rep", "linear_pred", "demand"))) %>%
#         ggplot(aes(
#             x = dist, 
#             y = value, 
#             colour = mu_type
#         )) +
#         facet_wrap(~variable, scales = "free") +
#         geom_point(alpha = 0.1) +
#         geom_line() +
#         theme_bw() +
#         theme(legend.position = "bottom") +
#         labs(
#             colour = "Visibility Type", 
#             x = "Distance (m)", 
#             y = "Value", 
#             title = str_glue("Fixing Private Incentive at {b}, Varying Visibility"), 
#             subtitle = str_glue("model: {model}, Cutoff at {cutoff_distance}km")
#         )

#     return(p)
# }


# long_subset_demand_df %>%
#         filter(cutoff_type == "cutoff") %>%
#         filter(b_type == mu_type)  %>%
#         filter(dist <= 3500) %>%
#         gather(variable, value, demand:v_star) %>%
#         mutate(variable = factor(variable,
#                                 levels = c("b", "v_star", "delta_v_star", "mu_rep", "linear_pred", "demand"))) %>%
#         ggplot(aes(
#             x = dist, 
#             y = value, 
#             colour = mu_type
#         )) +
#         facet_wrap(~variable, scales = "free") +
#         geom_point(alpha = 0.1) +
#         geom_line() +
#         theme_bw() +
#         theme(legend.position = "bottom")



# stan_fit_df %>%
#     select(cluster_w_cutoff) %>%
#     unnest()  %>%
#     filter(assigned_treatment == "bracelet") %>%
#     ggplot(aes(
#         x = roc_distance, 
#         y =  mean_est, 
#         colour = assigned_treatment
#     ))  +
#     geom_point() +
#     geom_line()




# stan_clean_df %>%
#     filter(
#         variable == "takeup_prop"
#     )  


# long_subset_demand_df %>%
#         filter(cutoff_type == "cutoff") %>%
#         filter(b_type == mu_type)  %>%
#         filter(dist <= 3500)  %>%
#         select(dist, demand, b_type) %>%
#         filter(dist < 100)




# stan_clean_df %>%
#     ggplot(aes(
#         x = dist, 
#         y = median_est, 
#         colour = treatment
#     )) +
#     facet_wrap(~variable, scales = "free") +
#     geom_point() +
#     geom_line()


# stan_clean_df %>%
#     filter(variable == "pred_takeup")

# long_subset_demand_df %>%
#         filter(cutoff_type == "cutoff") %>%
#         filter(b_type == mu_type)  %>%
#         filter(dist <= 3500) %>%
#         gather(variable, value, demand:v_star) %>%
#         filter(variable == "v_star") %>%
#         filter(b_type == "bracelet") %>%
#         mutate(variable = factor(variable,
#                                 levels = c("b", "v_star", "delta_v_star", "mu_rep", "linear_pred", "demand"))) %>%
#         ggplot(aes(
#             x = dist, 
#             y = value, 
#             colour = mu_type
#         )) +
#         geom_point(alpha = 0.1) +
#         geom_line() +
#         theme_bw() +
#         theme(legend.position = "bottom")

# stan_fit_df %>%
#     select(obs_cluster_takeup_level) %>%
#     unnest() %>%
#     filter(assigned_dist_obs < 200) %>%
#     filter(assigned_treatment == "bracelet") %>%
#     select(assigned_dist_obs, contains("per"))



# long_subset_demand_df %>%
#         filter(cutoff_type == "cutoff") %>%
#         filter(b_type == "control")  %>%
#         filter(dist <= 3500) %>%
#         gather(variable, value, demand:v_star) %>%
#         mutate(variable = factor(variable,
#                                 levels = c("b", "v_star", "delta_v_star", "mu_rep", "linear_pred", "demand"))) %>%
#         ggplot(aes(
#             x = dist, 
#             y = value, 
#             colour = mu_type
#         )) +
#         facet_wrap(~variable, scales = "free") +
#         geom_point(alpha = 0.1) +
#         geom_line() +
#         theme_bw() +
#         theme(legend.position = "bottom")
# ggsave(
#     "temp-plots/b-control-mu-vary.png",
#     width = 10,
#     height = 10
# )

# long_subset_demand_df %>%
#         filter(cutoff_type == "cutoff") %>%
#         filter(mu_type == "control")  %>%
#         filter(dist <= 3500) %>%
#         gather(variable, value, demand:v_star) %>%
#         mutate(variable = factor(variable,
#                                 levels = c("b", "v_star", "delta_v_star", "mu_rep", "linear_pred", "demand"))) %>%
#         ggplot(aes(
#             x = dist, 
#             y = value, 
#             colour = b_type
#         )) +
#         facet_wrap(~variable, scales = "free") +
#         geom_point(alpha = 0.1) +
#         geom_line() +
#         theme_bw() +
#         theme(legend.position = "bottom")

# ggsave(
#     "temp-plots/b-vary-mu-control.png",
#     width = 10,
#     height = 10
# )

# long_subset_demand_df %>%
#         filter(cutoff_type == "cutoff") %>%
#         filter(b_type == mu_type) %>%
#         filter(dist <= 2500) %>%
#         ggplot(aes(
#             x = dist, 
#             y = demand, 
#             colour = b_type
#         )) +
#         geom_point() +
#         geom_vline(
#             xintercept = 2500
#         )

# long_subset_demand_df %>%
#         filter(cutoff_type == "cutoff") %>%
#         filter(b_type == mu_type) %>%
#         filter(dist <= 3500) %>%
#         gather(variable, value, demand:v_star) %>%
#         mutate(
#             variable = factor(variable,
#                                 levels = c(
#                                     "b", 
#                                     "v_star", 
#                                     "delta_v_star", 
#                                     "mu_rep", 
#                                     "linear_pred", 
#                                     "demand"))
#         )  %>%
#         ggplot(aes(
#             x = dist, 
#             y = value, 
#             colour = mu_type
#         )) +
#         facet_wrap(~variable, scales = "free") +
#         geom_point(alpha = 0.1) +
#         geom_line() +
#         theme_bw() +
#         theme(legend.position = "bottom")



# p_panels_normal = plot_demand_panels(long_subset_demand_df, "cutoff", "STRUCTURAL_LINEAR_U_SHOCKS", b = "calendar")

# p_panels_normal

# p_panels_normal = plot_demand_panels(long_subset_demand_df, "no_cutoff", "STRUCTURAL_LINEAR_U_SHOCKS")
# p_panels_log = plot_demand_panels(long_subset_demand_df, "no_cutoff", "STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP")

# p_panels_log %>%
#     ggsave(
#         plot = .,
#         filename = "temp-plots/panel-plot-b-control-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP.png", 
#         width = 8,
#         height = 6, 
#         dpi = 500)

# p_panels_normal %>%
#     ggsave(
#         plot = .,
#         filename = "temp-plots/panel-plot-b-control-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS.png", 
#         width = 8,
#         height = 6, 
#         dpi = 500)

# long_subset_demand_df %>%
#     filter(model_type == "STRUCTURAL_LINEAR_U_SHOCKS") %>%
#     filter(b_type == mu_type) %>%
#     filter(dist < 2500) %>%
#     mutate(mu_type = factor(mu_type, levels = c("control", "ink", "calendar", "bracelet"))) %>%
#     ggplot(aes(
#         x = dist, 
#         y = mu_rep, 
#         colour = mu_type
#     )) +
#     geom_point() +
#     geom_hline(yintercept = 0.1) +
#     scale_y_continuous(breaks = seq(from = 0, to = 0.6, by = 0.1))

# ggsave("temp-plots/tmp.png", width = 8, height = 6, dpi = 500)

# long_subset_demand_df %>%
#     filter(
#         cutoff_type == "cutoff", 
#         b_type == "control", 
#         model_type == "STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP"
#     ) %>%
#     filter(dist <= 2500) %>%
#     ggplot(aes(
#         x = dist, 
#         y = mu_rep, 
#         colour = mu_type
#     )) +
#     geom_point()

# stop()



# optim_df = tibble(
#     optim_file = optim_files
# ) %>%
#     mutate(
#         optim_df = map(optim_files, read_rds)
#     )


# optim_df = optim_df %>%
#     mutate(
#         cutoff_type = if_else(
#             str_detect(optim_file, "\\/cutoff"),
#             "cutoff", 
#             "no_cutoff"
#         )
#         # b_type = str_extract(optim_file, "(?<=b-)(.*?)(?=-)"), 
#         # mu_type = str_extract(optim_file, "(?<=mu-)(.*?)(?=-)"), 
#         # model_type = str_extract(optim_file, str_glue("(?<=mu-{mu_type}-)(.*?)(?=-median)"))
#     )   %>%
#     select(-optim_file)



# # optim_no_rep_data = read_rds("optim/data/archive/cutoff-no-rep-median-optimal-allocation.rds") %>%
# #     as_tibble()
# # optim_rep_data = read_rds("optim/data/archive/cutoff-rep-median-optimal-allocation.rds") %>%
# #     as_tibble()
# village_data = read_csv("optim/data/village-df.csv")
# pot_data = read_csv("optim/data/pot-df.csv")

# ## Analysis data
# load(file.path("data", "analysis.RData"))

# standardize <- as_mapper(~ (.) / sd(.))
# unstandardize <- function(standardized, original) standardized * sd(original)
# # stick to monitored sms.treatment group
# # remove sms.treatment.2
# monitored_nosms_data <- analysis.data %>% 
#   filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
#   left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
#             by = "cluster.id") %>% 
#   mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
#   group_by(cluster.id) %>% 
#   mutate(cluster_id = cur_group_id()) %>% 
#   ungroup()

# analysis_data <- monitored_nosms_data %>% 
# mutate(
#     assigned_treatment = assigned.treatment, 
#     assigned_dist_group = dist.pot.group,
#     sms_treatment = factor(sms.treatment.2))

# ##



# optim_data = optim_df %>%
#     unnest(optim_df) %>%
#     rename(
#         b_type = private_benefit_z, 
#         mu_type = visibility_z, 
#         model_type = model)




# optim_data %>%
#     filter(model_type == "STRUCTURAL_LINEAR_U_SHOCKS") %>%
#     filter(b_type == "calendar", mu_type == "calendar") %>%
#     filter(cutoff_type == "cutoff") %>%
#     unnest(model_output) %>%
#     summarise(
#         n_pot = n_distinct(j), 
#         takeup = mean(demand)
#     )

# optim_data %>%
#     filter(model_type == "STRUCTURAL_LINEAR_U_SHOCKS") %>%
#     filter(
#         b_type == "control" , cutoff_type == "cutoff"
#     ) %>%
#     unnest(model_output) %>%
#     select( 
#         i, j, 
#         b_type, 
#         mu_type, 
#         demand, 
#         dist
#     ) %>%
#     spread(mu_type, demand)   %>%
#     summarise(
#         mean_diff = 100*mean(bracelet - control, na.rm = TRUE)
#     )

# optim_data %>%
#     filter(model_type == "STRUCTURAL_LINEAR_U_SHOCKS") %>%
#     filter(
#         b_type == "control" , cutoff_type == "cutoff"
#     ) %>%
#     unnest(model_output) %>%
#     filter(demand != 0) %>%
#     select( 
#         i, j, 
#         b_type, 
#         mu_type, 
#         demand, 
#         dist
#     ) %>%
#     spread(mu_type, demand)  %>%
#     ggplot(aes(
#         x = control, 
#         y = bracelet
#     )) +
#     geom_point() +
#     geom_abline(linetype = "longdash")



# optim_data = optim_data %>%
#     mutate(
#         n_pot = map_dbl(
#             model_output, 
#             ~{
#                 n_distinct(.x$j)

#             }
#         )
#     )


# optim_df

# optim_data

# optim_data %>%
#     filter(model_type == "STRUCTURAL_LINEAR_U_SHOCKS") %>%
#     select( 
#         b_type, 
#         mu_type, 
#         model_type, 
#         n_pot
#     )

# comp_df = bind_rows(
#     subset_demand_df %>%
#         filter(
#             model_type == "STRUCTURAL_LINEAR_U_SHOCKS"
#         ) %>%
#         filter(b_type == "control") %>%
#         filter(cutoff_type == "cutoff") %>%
#         unnest(subset_df) %>%
#         filter(dist < 2500) %>%
#         mutate(df_type = "demand"), 
#     optim_data %>%
#         filter(
#             model_type == "STRUCTURAL_LINEAR_U_SHOCKS"
#         ) %>%
#         filter(b_type == "control") %>%
#         filter(cutoff_type == "cutoff") %>%
#         unnest(demand_data) %>%
#         filter(dist < 2500) %>%
#         mutate(df_type = "optim") 
# )

# comp_df %>%
#     ggplot(aes(
#         x = dist, 
#         y = demand, 
#         colour = mu_type
#     )) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~df_type) +
#     labs(
#         title = "Control Private Incentive, Varying Visibility", 
#         y = "Takeup",
#         x = "Distance, Meters"
#     ) +
#     theme(legend.position = "bottom")

# subset_demand_df %>%
#     filter(
#         model_type == "STRUCTURAL_LINEAR_U_SHOCKS"
#     ) %>%
#     filter(b_type == "control") %>%
#     filter(cutoff_type == "cutoff") %>%
#     unnest(subset_df) %>%
#     filter(dist < 2500) %>%
#     ggplot(aes(
#         x = dist, 
#         y = demand, 
#         colour = mu_type
#     )) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~cutoff_type) +
#     labs(
#         title = "Control Private Incentive, Varying Visibility", 
#         y = "Takeup",
#         x = "Distance, Meters"
#     ) +
#     theme(legend.position = "bottom")

# optim_data %>%
#     filter(
#         model_type == "STRUCTURAL_LINEAR_U_SHOCKS"
#     ) %>%
#     filter(b_type == "control") %>%
#     filter(cutoff_type == "cutoff") %>%
#     unnest(demand_data) %>%
#     filter(dist < 2500) %>%
#     ggplot(aes(
#         x = dist, 
#         y = demand, 
#         colour = mu_type
#     )) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~cutoff_type) +
#     labs(
#         title = "Control Private Incentive, Varying Visibility", 
#         y = "Takeup",
#         x = "Distance, Meters"
#     ) +
#     theme(legend.position = "bottom")


# optim_data %>%
#     head(1) %>%
#     select(model_output) %>%
#     unnest() %>%
#     colnames()




# optim_data %>%
#     filter( 
#         model_type == "STRUCTURAL_LINEAR_U_SHOCKS"
#     ) %>%
#     # filter(
#     #     treatment %in% c("bracelet", "control")
#     # ) %>%
#     unnest(demand_data) %>%
#     filter(dist < 2500) %>%
#     ggplot(aes(
#         x = dist, 
#         y = demand,
#         color = mu_type
#     )) +
#     geom_point() +
#     facet_wrap(~b_type) +
#     geom_hline(
#         yintercept = 0.33, 
#         linetype = "longdash"
#     ) +
#     theme_bw()






# summ_analysis_data = analysis_data %>%
#     group_by(cluster.id) %>%
#     summarise(
#         demand = mean(dewormed), 
#         dist = unique(cluster.dist.to.pot), 
#         treatment = unique(assigned_treatment)
#     )   
# summ_analysis_fit = lm(
#     data = analysis_data %>% rename(dist = cluster.dist.to.pot, treatment = assigned_treatment),
#     formula = dewormed ~ dist + treatment
# )

# optim_data %>%
#     filter( 
#         model == "STRUCTURAL_LINEAR_U_SHOCKS"
#     ) %>%
#     filter(
#         treatment %in% c("bracelet", "control")
#     ) %>%
#     unnest(demand_data) %>%
#     filter(dist < 2500) %>%
#     modelr::add_predictions(summ_analysis_fit) %>%
#     ggplot(aes(
#         x = dist, 
#         y = demand,
#         color = treatment
#     )) +
#     geom_point() +
#     facet_wrap(~rep) +
#     geom_hline(
#         yintercept = 0.33, 
#         linetype = "longdash"
#     ) +
#     theme_bw() +
#     geom_line(aes(y = pred)) +
#     labs(
#         title = "Median Demand vs Distance", 
#         subtitle = "Bracelet vs Control, Fixing Rep at Control vs Treatment Rep", 
#         caption = str_wrap(
#             "Points show output from structural model. Line shows frequentist OLS on 'raw' data", 
#             width = 140
#         )    
#     ) +
#     theme(legend.position = "bottom")

# ggsave(
#     "temp-plots/control-treat-vis-with-ols-too.png",
#     width = 8,
#     height = 6
# )
# summ_analysis_fit %>%
#     tidy()


# 0.460 + 0.0821

# optim_data %>%
#     filter( 
#         model == "STRUCTURAL_LINEAR_U_SHOCKS"
#     ) %>%
#     filter(
#         treatment %in% c("bracelet", "control")
#     ) %>%
#     unnest(demand_data) %>%
#     filter(dist < 2500)  %>%
#     filter(demand > 0.30 & demand < 0.36) %>%
#     group_by(
#         rep, treatment
#     ) %>%
#     summarise(
#         mean_dist = mean(dist)
#     )
