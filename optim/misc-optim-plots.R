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
bracelet_oa_files = str_glue("optim/data/{script_options$model}/agg-log-full-many-pots/target-rep-cutoff-b-bracelet-mu-bracelet-{script_options$model}-median-optimal-allocation.rds")


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
    ylim(0, 7.1e-4) +
    xlim(0, 3000)
}


c1 = comp_df %>%
    mutate(dist = if_else(type != "Experimental", Inf, dist)) %>%
    plot_fun()

c2 = comp_df %>%
    mutate(dist = if_else(!(type %in% c("Experimental", "Optimal: Control")), Inf, dist)) %>%
    plot_fun()

c3 = comp_df %>%
    plot_fun()


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
        demand_data = map(demand_data, ~filter(.x, dist_km <= 3.5))
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


summ_subset_demand_df %>%
    mutate(mu_type = str_to_title(mu_type)) %>%
    mutate(
        type = paste0(v_star_type, ": ", mu_type)
    ) %>%
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



ggsave(
    file.path(
        script_options$output_path,
        str_glue(
            "{script_options$model}-agg-log-full-many-pots-pred-demand-vstar-comp.pdf"
        )
    ),
    width = 8, 
    height = 6, dpi = 500
)


