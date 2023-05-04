#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        misc-optim-plots.R  [options] 

        Options:
        --model=<model>
        --output-path=<output-path>
        --fit-version=<fit-version>
"),
  args = if (interactive()) "
                            --output-path=optim/plots/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-log-full-many-pots \
                            --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP \
                            --fit-version=86
                             " else commandArgs(trailingOnly = TRUE)
) 

library(tidyverse)
library(data.table)
library(sf)
library(tidybayes)

control_oa_files = str_glue("optim/data/{script_options$model}/agg-log-full-many-pots/target-rep-cutoff-b-control-mu-control-{script_options$model}-median-optimal-allocation.rds")
bracelet_oa_files = str_glue("optim/data/{script_options$model}/agg-log-full-many-pots/target-rep-cutoff-b-control-mu-bracelet-{script_options$model}-median-optimal-allocation.rds")


control_oa_df = read_rds(control_oa_files)
bracelet_oa_df = read_rds(bracelet_oa_files)



dist_data = read_rds(
  file.path(
    "optim/data",
    "full-many-pots-experiment.rds"
  )
)

demand_df = read_csv(
    file.path(
        str_glue("optim/data/{script_options$model}/agg-log-full-many-pots"),
        str_glue("pred-demand-dist-fit{script_options$fit_version}-cutoff-b-control-mu-control-{script_options$model}.csv")
    )
)

median_demand_df = demand_df %>%
    group_by(
        pot_j, 
        village_i
    ) %>%
    summarise(
        demand = median(demand), 
        dist = unique(dist)
    )

## Loading Data
wgs.84 = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
village_data = dist_data$village_df %>%
    st_as_sf(coords = c("lon", "lat"), crs = wgs.84)
pot_data = dist_data$pot_df %>%
    st_as_sf(coords = c("lon", "lat"), crs = wgs.84)

n = nrow(village_data)
m = nrow(pot_data)



pots_used_in_experiment = unique(village_data$cluster.id)


library(scales)
#extract hex color codes for a plot with three elements in ggplot2 
hex <- ggthemes::canva_pal("Primary colors with a vibrant twist")(2)

village_data %>%
        ggplot() +
        geom_sf() +
        geom_sf(
            data = pot_data,
            colour = hex[1],
            shape = 17,
            # size = 2, 
            alpha = 0.7) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        labs(
            title = "Optimal PoT Allocation Problem",
            subtitle = "Black dots indicate villages. Triangles indicate potential clinic locations."
        )  +
        labs(x = "", y = "") + 
        NULL

optimal_data = pot_data %>%
    filter(cluster.id %in% village_data$cluster.id)

optimal_data = inner_join(
    pot_data %>%
     dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
        select(
            cluster.id, 
            pot_id = id,
            pot_lon = lon, 
            pot_lat = lat,
            ) %>%
            as_tibble() %>%
            mutate(cluster.id = as.numeric(cluster.id)), 
    village_data %>% 
     dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
        select(
            cluster.id, 
            village_id = id,
            village_lon = lon, 
            village_lat = lat) %>%
            as_tibble() %>%
            mutate(cluster.id = as.numeric(cluster.id)),
    by = "cluster.id"
) 


optimal_data = optimal_data %>%
    left_join(
        median_demand_df, 
        by = c(
            "village_id" = "village_i", 
            "pot_id" = "pot_j"
    ))



comp_df = bind_rows(
    control_oa_df %>%
        pull(model_output) %>%
        first() %>%
        select(
            pot_j = j, 
            village_i = i, 
            pot_lon, pot_lat, 
            village_lon, village_lat,
            dist, 
            demand
        ) %>% mutate(type = "Optimal: Control"), 
    bracelet_oa_df %>%
        pull(model_output) %>%
        first() %>%
        select(
            pot_j = j, 
            village_i = i, 
            pot_lon, pot_lat, 
            village_lon, village_lat,
            dist, 
            demand
        ) %>% mutate(type = "Optimal: Bracelet"), 
    optimal_data %>%
        select(
            pot_j = pot_id, 
            village_i = village_id,
            pot_lon, pot_lat, 
            village_lon, village_lat,
            dist, 
            demand
        ) %>%
        mutate(type = "Experimental")
)


comp_df = comp_df %>%
    mutate(
        type = factor(type, levels = c("Experimental", "Optimal: Control", "Optimal: Bracelet"))
    )




plot_fun = function(data) {
    data %>%
    ggplot(aes(
        x = dist/1000, 
        fill = type
    )) +
    geom_density(
        colour = "black", 
        alpha = 0.5
        ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        fill = "", 
        y = "Density",
        x = "Distance Walked (km)"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    ylim(0, 0.75) +
    xlim(0, 3.5)
}


c1 = comp_df %>%
    mutate(dist = if_else(type != "Experimental", Inf, dist)) %>%
    plot_fun()
c2 = comp_df %>%
    mutate(dist = if_else(!(type %in% c("Experimental", "Optimal: Control")), Inf, dist)) %>%
    plot_fun() +
    annotate(
        "text", 
        x = 0.5, 
        y = 0.7, 
        label = "Amplification",
        size = 5, 
        alpha = 0.7
    )

c3 = comp_df %>%
    plot_fun() +
    annotate(
        "text", 
        x = 0.5, 
        y = 0.7, 
        label = "Amplification",
        size = 5, 
        alpha = 0.7
    ) +
    annotate(
        "text", 
        x = 3, 
        y = 0.25, 
        label = "Mitigation",
        size = 5, 
        alpha = 0.7
    )


imap(
    list(c1, c2, c3), 
    ~ggsave(
        plot = .x,
        filename = file.path(
            script_options$output_path,
            str_glue(
                "comp-dist-plot{.y}-fit{script_options$fit_version}-{script_options$model}.pdf"
            )
        ),
        width = 8, 
        height = 6
        )
)


rm(list = c(
    "comp_df", 
    "demand_df",
    "optimal_data",
    "median_demand_df"
))

gc()



# Amplification/Mitigation demand curve plots

demand_files = fs::dir_ls(
    str_glue("optim/data/{script_options$model}/agg-log-full-many-pots"),
    regexp = "pred"
)



demand_file_df = tibble(
    file = demand_files
) %>%
    mutate(
        b_type = str_extract(file, "(?<=-b-)\\w+(?=-mu)"),
        mu_type = str_extract(file, "(?<=-mu-)\\w+(?=-)"),
    )



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
    filter(rep_type == "rep")


subset_demand_df = subset_demand_df %>%
    filter(
        (mu_type == "bracelet" & v_star_type == "v_star") |
        (mu_type == "control" & v_star_type == "v_star") |
        (mu_type == "control" & v_star_type == "static") 
    )

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
            "plot{.y}-{script_options$model}-agg-log-full-many-pots-pred-demand-vstar-comp.pdf"
        )
    ),
    width = 8, 
    height = 6, dpi = 500
)
)




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

plot_summ_bc_draws %>%
    write_csv("temp-data/plot-summ-draws-roc-vis-no-vis.csv")


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
    