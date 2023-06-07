#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        misc-optim-plots.R  [options] 

        Options:
        --model=<model>
        --output-path=<output-path>
        --fit-version=<fit-version>
        --welfare-function=<welfare-function> Which utility function to use [default: log]
        --distance-constraint=<distance-constraint>  Distance constraint, in meters [default: 3500]
"),
  args = if (interactive()) "
                            --output-path=optim/plots/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots \
                            --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP \
                            --fit-version=86 \
                            --welfare-function=identity
                             " else commandArgs(trailingOnly = TRUE)
) 
library(tidyverse)
library(tidybayes)
library(posterior)
library(cmdstanr)
library(ggthemes)

# Amplification/Mitigation demand curve plots

demand_files = fs::dir_ls(
    str_glue("optim/data/{script_options$model}/agg-full-many-pots"),
    regexp = "pred"
)


demand_file_df = tibble(
    file = demand_files
) %>%
    mutate(
        b_type = str_extract(file, "(?<=-b-)\\w+(?=-mu)"),
        mu_type = str_extract(file, "(?<=-mu-)\\w+(?=-)"),
    ) %>%
    filter(str_detect(file, paste0(script_options$model, ".csv")))



subset_demand_df = demand_file_df %>%
    filter(b_type == "control") %>%
    mutate(
        rep_type = if_else(
            str_detect(file, "suppress-rep"), 
            "suppress_rep", 
            "rep"
        ), 
        v_star_type = if_else(
            str_detect(file, "static"), 
            "static", 
            "v_star"
        )
    )  %>%
    filter(rep_type == "rep") %>%
    mutate(
        cutoff_type = if_else(
                str_detect(file, "no-cutoff"), 
                "no_cutoff", 
                "cutoff"
        )
    )


subset_demand_df = subset_demand_df %>%
    filter(
        (mu_type == "bracelet" & v_star_type == "v_star") |
        (mu_type == "control" & v_star_type == "v_star") |
        (mu_type == "control" & v_star_type == "static") 
    )  %>%
    filter(cutoff_type == "cutoff")


subset_demand_df = subset_demand_df %>%
    mutate(
        demand_data = map(file, read_csv)
    )

subset_demand_df = subset_demand_df %>%
    mutate(
        demand_data = map(demand_data, ~filter(.x, dist_km <= 2.5))
    )

summ_subset_demand_df = subset_demand_df %>%
    mutate(
        summ_demand_data = map(
            demand_data, 
            ~ { .x %>%
                group_by(dist_km) %>%
                median_qi(demand, .width = c(0.80, 0.50))
            }
        )
    )
    



summ_subset_demand_df = summ_subset_demand_df %>%
    unnest(summ_demand_data) %>%
    mutate(
        v_star_type = if_else(v_star_type == "static", "Static", "W*")
    ) 


plot_amp_mit_fun = function(data) {
    data %>%
        ggplot(aes(
            x = dist_km, 
            y = demand,
            ymin = .lower, 
            ymax = .upper, 

        )) +
        geom_line(aes(
            colour = type, 
            group = type, 
            linetype = v_star_type
        )) +
        geom_ribbon(
            data = . %>%
                filter(.width == 0.5),
            aes(fill = type), alpha = 0.3) +
        geom_ribbon(
            data = . %>%
                filter(.width == 0.8),
            aes(fill = type), alpha = 0.3) +
        theme_minimal() +
        theme( 
            legend.position = "bottom",
            legend.title = element_blank()
        )  +
        labs(
            x = "Distance (km)", 
            y = "Estimated Takeup", 
            colour = ""
        ) +
        guides(
            fill = "none", 
            linetype = "none"
        ) +
        ggthemes::scale_color_canva( 
            palette = "Primary colors with a vibrant twist"
        ) +
        ggthemes::scale_fill_canva( 
            palette = "Primary colors with a vibrant twist"
        ) +
        labs(
            caption = "Line: Median. Outer ribbon: 80% credible interval. Inner ribbon: 50% credible interval." 
        ) +
        scale_linetype_manual(
            values = c("longdash", "solid")
        )

}

plot_summ_subset_demand_df = summ_subset_demand_df %>%
    mutate(mu_type = str_to_title(mu_type)) %>%
    mutate(
        type = paste0(v_star_type, ": ", mu_type)
    ) 

full_p_amp_mit = plot_summ_subset_demand_df %>%
    plot_amp_mit_fun() 

ylim_p_amp_mit = ggplot_build(full_p_amp_mit)$layout$panel_scales_y[[1]]$range$range
full_p_amp_mit = plot_summ_subset_demand_df %>%
    plot_amp_mit_fun()  +
    ylim(c(0, ylim_p_amp_mit[2]))

first_p_amp_mit = plot_summ_subset_demand_df %>%
    filter(v_star_type == "Static") %>%
    plot_amp_mit_fun() +
    ylim(c(0, ylim_p_amp_mit[2]))

second_p_amp_mit = plot_summ_subset_demand_df %>%
    filter(v_star_type == "Static" | mu_type == "Bracelet") %>%
    plot_amp_mit_fun() +
    ylim(c(0, ylim_p_amp_mit[2]))


imap(
    list(
        first_p_amp_mit, 
        second_p_amp_mit,
        full_p_amp_mit
    ), 
    ~ggsave(
        plot = .x,
    file.path(
        script_options$output_path,
        str_glue(
            "plot{.y}-{script_options$model}-agg-{script_options$welfare_function}-full-many-pots-pred-demand-vstar-comp.pdf"
        )
    ),
    width = 8, 
    height = 6, dpi = 500
)
)





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

plot_summ_bc_draws %>%
    write_csv("temp-data/plot-summ-draws-roc-vis-no-vis.csv")

