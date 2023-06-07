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
          --model=<model>
          --fit-version=<fit-version>
          --distance-constraint=<distance-constraint>  Distance constraint, in meters [default: 3500]
"),
  args = if (interactive()) "
                            --constraint-type=agg \
                            --welfare-function=identity \
                            --min-cost \
                            --output-path=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots \
                            --output-basename=target-rep-agg-identity-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median \
                            --cutoff-type=cutoff
                            --data-input-name=full-many-pots-experiment.rds \
                            --posterior-median \
                            --pdf-output-path=presentations/takeup-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-fig/
                            --demand-input-path=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots \
                            --demand-input-filename=pred-demand-dist-fit86-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP.csv
                             " else commandArgs(trailingOnly = TRUE)
) 

library(tidyverse)
library(sf)
library(data.table)
library(latex2exp)
library(tidybayes)



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


swf = eval(parse(text = script_options$welfare_function))

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


optimal_data %>%
    saveRDS(
        file.path(
            script_options$output_path, 
            str_glue(
                "{script_options$output_basename}-experimental-control-allocation-data.rds"
            )
        )
    )


swf_summ_df  = optimal_data %>%
    group_by(draw) %>%
        mutate(
            social_welfare = sum(swf(demand))
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
        # theme(legend.position = "bottom") +
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


#### Distance Walked Histogram ####


control_oa_files = str_glue(
    "optim/data/{script_options$model}/agg-full-many-pots/target-rep-distconstraint-{script_options$distance_constraint}-util-{script_options$welfare_function}-cutoff-b-control-mu-control-{script_options$model}-median-optimal-allocation.rds"
    )
bracelet_oa_files = str_glue(
    "optim/data/{script_options$model}/agg-full-many-pots/target-rep-distconstraint-{script_options$distance_constraint}-util-{script_options$welfare_function}-cutoff-b-control-mu-bracelet-{script_options$model}-median-optimal-allocation.rds"
    )

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
        str_glue("optim/data/{script_options$model}/agg-full-many-pots"),
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


max_x = comp_df %>%
    summarise(md = max(dist, na.rm = TRUE)) %>%
    pull()

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
    xlim(0, (max_x + 100)/1000)
}

c3 = comp_df %>%
    plot_fun() +
    annotate(
        "text", 
        x = 0.5, 
        y = 0.6, 
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

y_lim = layer_scales(c3)$y$range$range

c1 = comp_df %>%
    mutate(dist = if_else(type != "Experimental", Inf, dist)) %>%
    plot_fun() +
    ylim(y_lim) 
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
    ) +
    ylim(y_lim) 

imap(
    list(c1, c2, c3), 
    ~ggsave(
        plot = .x,
        filename = file.path(
            script_options$output_path,
            str_glue(
                "comp-dist-plot{.y}-fit{script_options$fit_version}-distconstraint-{script_options$distance_constraint}-util-{script_options$welfare_function}-{script_options$model}.pdf"
            )
        ),
        width = 8, 
        height = 6
        )
)
