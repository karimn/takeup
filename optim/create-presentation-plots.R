#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        create-presentation-plots.R  [options] 

        Options:
          --min-cost  Flag to minimise programme cost for a given takeup level.
          --max-takeup  Flag to maximise vaccine takeup for a given cost level.
          --comp-output-basename=<comp-output-basename>  Comparison output basename
          --output-path=<output-path>  Path where output should be saved.
          --output-basename=<output-basename>  Output basename.
          --map-plot  Plot map without adding lines for closest villages and active PoTs 
          --posterior-median  If posterior median over takeup demand is used.
          --cutoff-type=<cutoff-type>  Which cutoff type to use if using full posterior draws [default: no-cutoff]
          --constraint-type=<constraint-type>  Aggregate or individual village welfare
          --welfare-function=<welfare-function>  Log or identity utility
          --data-input-path=<data-input-path>  Where village and PoT data is stored [default: {file.path('optim', 'data')}]
          --data-input-name=<data-input-name>  Filename of village and PoT data [default: full-experiment-4-extra-pots.rds]
          --pdf-output-path=<pdf-output-path>  Output path for PDF plots 
          --demand-input-path=<demand-input-path>  Path where demand data is stored.
          --demand-input-filename=<demand-input-filename> Demand input filename. 
"),
  args = if (interactive()) "
                            --constraint-type=agg \
                            --welfare-function=log \
                            --min-cost \
                            --output-path=optim/plots/agg-log-full-many-pots \
                            --output-basename=agg-log-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median \
                            --cutoff-type=cutoff
                            --data-input-name=full-many-pots-experiment.rds \
                            --posterior-median \
                            --pdf-output-path=presentations/takeup-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-fig/
                            --demand-input-path=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-log-full-many-pots \
                            --demand-input-filename=pred-demand-dist-fit86-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv
                             " else commandArgs(trailingOnly = TRUE)
) 

library(tidyverse)
library(sf)
library(data.table)
library(latex2exp)



source('optim/optim-functions.R')

optim_type = if_else(script_options$min_cost, "min_cost", "max_takeup")
stat_type = if_else(script_options$posterior_median, "median", "post-draws")



dist_data = read_rds(
  file.path(
    script_options$data_input_path,
    script_options$data_input_name
  )
)

demand_df = read_csv(
    file.path(
        script_options$demand_input_path,
        script_options$demand_input_filename
    )
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
        demand_df, 
        by = c(
            "village_id" = "village_i", 
            "pot_id" = "pot_j"
    ))


swf_summ_df  = optimal_data %>%
    group_by(draw) %>%
        mutate(
            social_welfare = sum(log(demand))
        ) 




plot_optimal_allocation = function(village_data,
                                    pot_data,
                                    optimal_data){
    summ_df = optimal_data %>%
    group_by(draw) %>%
        mutate(
            swf = sum(log(demand))
        ) %>%
        ungroup() %>%
        summarise(
            mean_demand = round(100*mean(demand), 2), 
            mean_dist = mean(dist), 
            n_pot = n_distinct(pot_id), 
            mean_swf = mean(swf)
        )


    takeup_hit = summ_df$mean_demand
    mean_dist = summ_df$mean_dist
    n_pots_used = summ_df$n_pot
    sw = summ_df$mean_swf
    pot_data$assigned_pot = pot_data$cluster.id %in% village_data$cluster.id

    title_str = str_glue(
        "Takeup: {takeup_hit}%"
    )

    subtitle = str_glue(
        "Average Distance: {round(mean_dist/1000, 2)}km, Social Welfare: {round(sw, 2)}"
    )
    sub_str = str_glue("Assigned PoTs: {n_pots_used}")


    plot_data = bind_rows(
        village_data %>% mutate(type = "village") %>% select(-cluster.id), 
        pot_data %>% mutate(type = paste0("pot_", assigned_pot)) %>% select(-cluster.id)
    )


    optimal_pot_plot = plot_data %>%
        ggplot() +
        geom_sf(aes(
            color = type, 
            shape = type, 
            size = type, 
            alpha = type
        )) +
        theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_segment(
            data = optimal_data, 
            aes(x = village_lon, 
                y = village_lat, 
                xend = pot_lon, 
                yend = pot_lat
                ))  +
        annotate(
            "text",
            x = 34.8, 
            y = 0.75,
            label = sub_str
        ) +
        labs(x = "", y = "", title = title_str, subtitle = subtitle) +
        scale_size_manual(
            values = c("pot_TRUE" = 1.5, "pot_FALSE" = 1 , "village" = 1),
            breaks = c("pot_TRUE", "pot_FALSE", "village"),
            labels = c("PoT Used", "PoT Unused", "Village")
        ) +
        scale_alpha_manual(
            values = c("pot_TRUE" = 1, "pot_FALSE" = 0.3, "village" = 0.8),
            breaks = c("pot_TRUE", "pot_FALSE", "village"),
            labels = c("PoT Used", "PoT Unused", "Village")
        ) +
        scale_color_manual(
            values = c("pot_TRUE" = hex[2], "pot_FALSE" = hex[1], "village" = "black"),
            breaks = c("pot_TRUE", "pot_FALSE", "village"),
            labels = c("PoT Used", "PoT Unused", "Village")
        ) +
        scale_shape_manual(
            values = c("pot_TRUE" = 17, "pot_FALSE" = 17, "village" = 16),
            breaks = c("pot_TRUE", "pot_FALSE", "village"),
            labels = c("PoT Used", "PoT Unused", "Village")
        ) +
        theme(legend.title = element_blank()) +
        theme(legend.position = "bottom") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        guides(alpha = guide_legend(override.aes = list(alpha = 1)))
    return(optimal_pot_plot)
}


control_in_experiment = plot_optimal_allocation(village_data, pot_data, optimal_data)
ggsave(
    file.path(
        script_options$pdf_output_path,
        str_glue(
            "{script_options$output_basename}-control-used-optimal-allocation-plot.pdf"
        )
    ),
    width = 8, 
    height = 6, 
    dpi = 500
)



