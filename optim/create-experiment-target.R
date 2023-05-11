#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        create-experiment-target.R  [options] 

        Options:
          --min-cost  Flag to minimise programme cost for a given takeup level.
          --max-takeup  Flag to maximise vaccine takeup for a given cost level.
          --comp-output-basename=<comp-output-basename>  Comparison output basename
          --output-path=<output-path>  Path where output should be saved.
          --output-basename=<output-basename>  Output basename.
          --posterior-median  If posterior median over takeup demand is used.
          --cutoff-type=<cutoff-type>  Which cutoff type to use if using full posterior draws [default: no-cutoff]
          --constraint-type=<constraint-type>  Aggregate or individual village welfare
          --welfare-function=<welfare-function>  Log or identity utility
          --data-input-path=<data-input-path>  Where village and PoT data is stored [default: {file.path('optim', 'data')}]
          --data-input-name=<data-input-name>  Filename of village and PoT data [default: full-experiment-4-extra-pots.rds]
          --demand-input-path=<demand-input-path>  Path where demand data is stored.
          --demand-input-filename=<demand-input-filename> Demand input filename. 
"),
  args = if (interactive()) "
                            --constraint-type=agg \
                            --welfare-function=identity \
                            --min-cost \
                            --output-path=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots \
                            --output-basename=summ-agg-identity \
                            --cutoff-type=cutoff
                            --data-input-name=full-many-pots-experiment.rds \
                            --posterior-median \
                            --demand-input-path=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots \
                            --demand-input-filename=pred-demand-dist-fit86-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv
                             " else commandArgs(trailingOnly = TRUE)
) 

library(tidyverse)
library(sf)
library(data.table)
library(latex2exp)


swf = eval(parse(text = script_options$welfare_function))

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
            social_welfare = sum(swf(demand))
        ) 
swf_summ_df %>%
    write_csv(
        file.path(
            script_options$output_path,
            paste0(script_options$output_basename,  "-experiment-target-constraint.csv")
        )
    )
