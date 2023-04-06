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

stop()


comp_df %>%
    select(dist, type) %>%
    group_by(type) %>%
    mutate(id = 1:n()) %>%
    spread(
        type, dist
    ) %>%
    ggplot(aes(
        x = Experimental, 
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
        x = "Distance Walked (m)"
    ) +
    scale_fill_brewer(palette = "Dark2")


c1 = comp_df %>%
    mutate(dist = if_else(type != "Experimental", Inf, dist)) %>%
    ggplot(aes(
        x = dist, 
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
        x = "Distance Walked (m)"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    ylim(0, 7.1e-4) +
    xlim(0, 3000)

c2 = comp_df %>%
    mutate(dist = if_else(!(type %in% c("Experimental", "Optimal: Control")), Inf, dist)) %>%
    ggplot(aes(
        x = dist, 
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
        x = "Distance Walked (m)"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    ylim(0, 7.1e-4) +
    xlim(0, 3000)

c3 = comp_df %>%
    ggplot(aes(
        x = dist, 
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
        x = "Distance Walked (m)"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    ylim(0, 7.1e-4) +
    xlim(0, 3000)



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


ggsave(
    plot = c1,
    filename = file.path(
        script_options$output_path,
        str_glue(
            "comp-dist-plot1-fit{script_options$fit_version}-{script_options$model}.pdf"
        )
    ),
    width = 8, 
    height = 6
)




# pot_we_want = 8

# villages_in_581 = comp_df %>%
#     filter(pot_j == pot_we_want) %>%
#     pull(village_i)


# exp_pots_we_want = comp_df %>%
#     filter(type == "experimental") %>%
#     filter(village_i %in% villages_in_581) %>%
#     pull(pot_j)

# comp_df %>%
#     filter(village_i %in% villages_in_581) %>%
#     select(village_i, pot_j, dist, demand, type)  %>%
#     gather(
#         variable, value, -type, -village_i
#     ) %>%
#     spread(variable, value)

# comp_df %>%
#     filter(village_i %in% villages_in_581) 

# pot_data %>%
#     colnames()



# village_data %>%
#         ggplot() +
#         geom_sf(alpha = 0.2) +
#         geom_sf(
#             data = pot_data,
#             colour = hex[1],
#             shape = 17,
#             # size = 2, 
#             alpha = 0.2) +
#         geom_sf(
#             data = pot_data %>%
#                 filter(id %in% exp_pots_we_want ),
#             colour = "red",
#             shape = 17,
#             size = 5, 
#             alpha = 1) +
#         geom_sf(
#             data = pot_data %>%
#                 filter(id %in%  pot_we_want ),
#             colour = "hotpink",
#             shape = 17,
#             size = 5, 
#             alpha = 1) +
#         geom_sf(
#             data = . %>%
#                 filter(id %in% villages_in_581), 
#             size = 2
#         ) +
#         theme_bw() +
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#         labs(
#             title = "Optimal PoT Allocation Problem",
#             subtitle = "Black dots indicate villages. Triangles indicate potential clinic locations."
#         )  +
#         labs(x = "", y = "") + 
#         NULL
