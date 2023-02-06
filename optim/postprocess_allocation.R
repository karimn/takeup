#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        postprocess_allocation.R  [options] [--min-cost | --max-takeup]

        Options:
          --min-cost  Flag to minimise programme cost for a given takeup level.
          --max-takeup  Flag to maximise vaccine takeup for a given cost level.
          --optim-input-path=<optim-input-path>  Path where input data is stored.
          --demand-input-path=<demand-input-path>  Path where input data is stored.
          --optim-input-a-filename=<optim-input-a-filename>  Optim input a filename.
          --demand-input-a-filename=<demand-input-a-filename>  Estimated demand a for every village i, PoT j pair.  
          --output-path=<output-path>  Path where output should be saved.
          --output-basename=<output-basename>  Output basename.
          --map-plot  Plot map without adding lines for closest villages and active PoTs 
          --posterior-median  If posterior median over takeup demand is used.
          --cutoff-type=<cutoff-type>  Which cutoff type to use if using full posterior draws [default: no-cutoff]
          --constraint-type=<constraint-type>  Aggregate or individual village welfare
          --welfare-function=<welfare-function>  Log or identity utility
          --data-input-path=<data-input-path>  Where village and PoT data is stored [default: {file.path('optim', 'data')}]
          --data-input-name=<data-input-name>  Filename of village and PoT data [default: full-experiment-4-extra-pots.rds]
"),
  args = if (interactive()) "
                            --constraint-type=agg \
                            --welfare-function=log \
                            --min-cost \
                            --optim-input-path=optim/data/agg-log-full-many-pots \
                            --demand-input-path=optim/data/agg-log-full-many-pots \
                            --optim-input-a-filename=cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS-post-draws-optimal-allocation.rds \
                            --demand-input-a-filename=pred-demand-dist-fit71-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS.csv \
                            --output-path=optim/plots/agg-log-full-many-pots \
                            --output-basename=agg-log-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS-post-draws \
                            --cutoff-type=cutoff
                            --data-input-name=full-many-pots-experiment.rds
                             
                             " else commandArgs(trailingOnly = TRUE)
) 

                            # --target-constraint=0.33
                            #  --min-cost 
                            #  --target-constraint=0.33
                            #  --input-path=optim/data
                            #  --optim-input-a-filename=-no-rep-median-optimal-allocation.rds
                            #  --village-input-filename=village-df.csv
                            #  --pot-input-filename=pot-df.csv
                            #  --demand-input-a-filename=pred_demand_dist_fit71_no_rep.csv
                            #  --output-path=optim
                            #  --output-basename=structural-median-no-rep-test
                            #  --map-plot
                            #  --posterior-median

library(tidyverse)
library(sf)
library(data.table)



# script_options = script_options %>%
#   modify_at(numeric_options, as.numeric)

optim_type = if_else(script_options$min_cost, "min_cost", "max_takeup")
stat_type = if_else(script_options$posterior_median, "median", "post-draws")


## Input Paths
optim_input_a_filepath = file.path(script_options$optim_input_path, 
                                script_options$optim_input_a_filename)
demand_input_a_path = file.path(script_options$demand_input_path, 
                            script_options$demand_input_a_filename)

## Input Data
optimal_df = read_rds(optim_input_a_filepath)
demand_a_data = read_csv(demand_input_a_path) %>%
    as.data.table()

dist_data = read_rds(
  file.path(
    script_options$data_input_path,
    script_options$data_input_name
  )
)

village_data = dist_data$village_df
pot_data = dist_data$pot_df

n = nrow(village_data)
m = nrow(pot_data)

## Output Filepaths
map_plot_path = file.path(
    script_options$output_path,
    paste0(
        script_options$output_basename,
        "-map-plot.png"
    )
)


optimal_allocation_demand_comp_plot_path = file.path(
    script_options$output_path,
    str_glue(
        "{script_options$output_basename}-optimal-allocation-demand-comp-plot.png"
    )
)



wgs.84 = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
## Loading Data
data = list(
    n = n,
    m = m,
    pot_locations = st_as_sf(pot_data, coords = c("lon", "lat"), crs = wgs.84),
    village_locations = st_as_sf(village_data, coords = c("lon", "lat"), crs = wgs.84)
) 



library(scales)
#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(2)


if (script_options$map_plot){
    data$village_locations %>%
        ggplot() +
        geom_sf() +
        geom_sf(
            data = data$pot_locations,
            colour = hex[1],
            shape = 17,
            size = 4, 
            alpha = 0.7) +
        theme_bw() +
        labs(
            title = "Optimal PoT Allocation Problem",
            subtitle = "Black dots indicate villages. Triangles indicate potential clinic locations."
        )  +
        labs(x = "", y = "") + 
        # theme(
        # axis.title.x=element_blank(),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        # axis.title.y=element_blank(),
        # axis.text.y=element_blank(),
        # axis.ticks.y=element_blank()) +
        NULL

    ggsave(map_plot_path, width = 8, height = 6, dpi = 500)
}


if (stat_type == "median") {

    plot_optimal_allocation = function(village_data,
                                    pot_data,
                                    optimal_data){

        assigned_pots = unique(optimal_data$j)
        n_pots_used =  length(assigned_pots)
        pot_data$assigned_pot = pot_data$id %in% assigned_pots
        summ_optimal_data = optimal_data %>%
            mutate(target_optim = target_optim) %>%
            summarise(
                util = sum(log(demand)),
                mean_demand = mean(demand), 
                min_demand = min(demand), 
                n_pot = n_distinct(j), 
                mean_dist = mean(dist),
                target_optim = mean(target_optim)
            ) %>%
            mutate(
                target_optim = target_optim
            ) %>%
            mutate(
                overshoot = 100*(util/target_optim - 1)
            ) 

        takeup_hit = round(summ_optimal_data$mean_demand*100,1 )
        util_hit = round(summ_optimal_data$util)
        util_target = round(summ_optimal_data$target_optim)
        overshoot = round(abs(summ_optimal_data$overshoot), 3)
        mean_dist = round(summ_optimal_data$mean_dist,1)

        if (script_options$constraint_type == "agg") {
            sub_str = str_glue("Takeup target: Aggregate Welfare Under Control. Assigned PoTs: {n_pots_used}")
        }

        if (script_options$constraint_type == "indiv") {
            sub_str = str_glue("Takeup target: Pareto Improving Allocation, Village Takeup Fixed at Control. Assigned PoTs: {n_pots_used}")
        }

        caption_str = ifelse(
            script_options$welfare_function == "log", 
            str_glue("Log utility used in social welfare function. Utility target: {util_target}, utility achieved: {util_hit}, percentage over constraint: {overshoot}%, takeup achieved: {takeup_hit}%. Mean walking distance: {mean_dist}m"),
            "Identity utility function used, i.e. takeup enters SWF directly."
            )

        optimal_pot_plot = village_data %>%
            ggplot() +
            geom_sf(alpha = 1) +
            geom_sf(
                data = pot_data %>% filter(assigned_pot == FALSE),
                color = hex[1],
                shape = 17,
                size = 4, 
                alpha = 0.2) +
            geom_sf(
                data = pot_data %>% filter(assigned_pot == TRUE),
                color = hex[2],
                shape = 17,
                size = 4, 
                alpha = 1) +
            theme_bw() +
            geom_segment(
                data = optimal_data, 
                aes(x = village_lon, 
                    y = village_lat, 
                    xend = pot_lon, 
                    yend = pot_lat
                    ))  +
            labs(x = "", y = "", subtitle = sub_str, caption = caption_str) +
            scale_color_manual(
                values = hex,
                labels = c("PoT Used", "PoT Unused")
            ) 
        return(optimal_pot_plot)
    }


    # plot_optimal_allocation(
    #     village_data = data$village_locations,
    #     pot_data = data$pot_locations,
    #     optimal_data = optimal_df %>%
    #         filter(treatment == "bracelet" & str_detect(model, "RED")) %>%
    #         unnest(model_output)
    # )


    allocation_plots = map(
        optimal_df$model_output,
        ~ plot_optimal_allocation(
            village_data = data$village_locations,
            pot_data = data$pot_locations,
            optimal_data = .x
        )
    )
    iwalk(
        allocation_plots,
        ~ggsave(
            plot = .x,
            filename = file.path(
                script_options$output_path,
                str_glue(
                    "{script_options$output_basename}-optimal-allocation-plot.png"
                )
            ),
        width = 10,
        height = 10,
        dpi = 500
        )
    )
}

if (stat_type == "post-draws") {

    long_optimal_df = optimal_df %>%
        select(
            draw,
            model,
            private_benefit_z,
            visibility_z, 
            model_output
            )  %>%
        unnest(model_output)

    long_optimal_df = long_optimal_df %>%
        select(-c(`Level of Education`:cluster.id.y))

    long_optimal_df %>%
        saveRDS(
            file.path(
                script_options$optim_input_path, 
                str_replace(script_options$optim_input_a_filename, "\\.rds", "-subset-long-data.rds")
            )
        )
}
