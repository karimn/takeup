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
                            --posterior-median \
                            --optim-input-path=optim/data/agg-log-busia \
                            --demand-input-path=optim/data/agg-log-busia \
                            --optim-input-a-filename=no-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP-median-optimal-allocation.rds \
                            --demand-input-a-filename=pred-demand-dist-fit71-no-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP.csv \
                            --output-path=optim/plots/agg-log-busia \
                            --output-basename=agg-log-no-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_LOG_MU_REP-median \
                            --cutoff-type=cutoff
                            --data-input-name=BUSIA-experiment.rds

                             
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

# oa_files = fs::dir_ls(
#     "optim/data/agg-log-busia",
#     regexp = "rds"
#     )

# subset_oa_files = oa_files[str_detect(oa_files, "\\/no-cutoff") & str_detect(oa_files, "LOG")]

# oa_df = subset_oa_files %>%
#     map_dfr(read_rds)

# long_oa_df = oa_df %>%
#     unnest(model_output) 

# long_oa_df %>%
#     select(model) %>%
#     unique() %>%
#     pull()

# stop()


# long_oa_df %>%
#     unnest(demand) %>%
#     ggplot(aes(
#         x = dist,
#         y = demand, 
#         colour = visibility_z
#     )) +
#     geom_point() +
#     geom_vline(xintercept =  3500)


# long_oa_df %>%
#     unnest(demand_data) %>%
#     colnames()

# long_oa_df %>%
#     unnest(demand) %>%
#     ggplot(aes(
#         x = dist,
#         y = mu, 
#         colour = visibility_z
#     )) +
#     geom_point()




# walk_df = long_oa_df %>%
#     filter(visibility_z != "calendar") %>%
#     group_by(
#         i
#     ) %>%
#     mutate(
#         n_dists = n_distinct(dist) 
#     ) %>%
#     filter(n_dists > 1)  %>%
#     select(i, visibility_z, dist) %>%
#     spread(
#         visibility_z, dist
#     ) %>%
#     mutate( 
#         walk_increase = bracelet - control
#     )

# walk_df %>%
#     ungroup() %>%
#     summarise(
#         mean_increase = mean(walk_increase)
#     )

# long_oa_df



# long_oa_df %>%
#     ungroup() %>%
#     group_by(visibility_z) %>%
#     summarise(
#         pct_close = mean(dist < 1250),
#         pct_far = mean(dist > 1250)
#     )

# walk_df


# walk_df %>%
#     ungroup() %>%
#     summarise(
#         pct_close = mean(control < 1250), 
#         pct_far = mean(control > 1250), 
#     )


# walk_df     %>%
#     ggplot(aes(
#         x = walk_increase
#     )) +
#     geom_histogram() +
#     geom_vline(xintercept = 0, linetype = "longdash")


# long_oa_df %>%
#     group_by(visibility_z) %>%
#     summarise(
#         dist = mean(dist)
#     )

# long_oa_df %>%
#     filter(dist > 0) %>%
#     ggplot(aes(
#         x = dist,
#         fill = visibility_z
#     )) +
#     geom_histogram()


# long_oa_df %>%
#     ggplot(aes(
#         x = dist,
#         fill = visibility_z
#     )) +
#     geom_histogram() +
#     facet_wrap(~visibility_z, ncol = 1)



# optimal_allocation_demand_comp_plot_path = file.path(
#     script_options$output_path,
#     str_glue(
#         "{script_options$output_basename}-optimal-allocation-demand-comp-plot.png"
#     )
# )


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


        if (script_options$constraint_type == "agg") {
            sub_str = str_glue("Takeup target: Aggregate Welfare Under Control. Assigned PoTs: {n_pots_used}")
        }

        if (script_options$constraint_type == "indiv") {
            sub_str = str_glue("Takeup target: Pareto Improving Allocation, Village Takeup Fixed at Control. Assigned PoTs: {n_pots_used}")
        }

        caption_str = ifelse(
            script_options$welfare_function == "log", 
            "Log utility used in social welfare function",
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
    library(tidyverse)
    library(sf)

    cutoff_type = script_options$cutoff_type 

    draw_files = fs::dir_ls("optim/data/", regexp = "cutoff.*post-draws.*\\.rds$")


    if (cutoff_type == "cutoff") {
        draw_files = draw_files[str_detect(draw_files, "data\\/cutoff")]
    } else {
        draw_files = draw_files[str_detect(draw_files, "data\\/no-cutoff")]
    }



    draw_data = map(
        draw_files,
        read_rds
    ) %>%
        map(as_tibble) %>%
        imap(

            ~mutate(., file = draw_files[.y])
        ) 



    draw_data = map_dfr(
        draw_data,
        ~mutate(
            .x,
            cutoff_type = cutoff_type,
            rep_type = if_else(str_detect(file, "no-rep"), "no-rep", "rep")
        )
    ) %>%
    select(-file)

    draw_data = draw_data %>%
        filter(!str_detect(model, "REDUCED"))

    draw_data = draw_data %>%
        mutate(
            n_pots = map_dbl(
                model_output,
                ~n_distinct(.x$j)
            )
        )

    library(ggthemes)


    canva_palette_vibrant <- "Primary colors with a vibrant twist"



    plot_post_estimates = function(data, treatment_arm, cutoff_type) {
        
        plot_title = 
            str_glue("Posterior Histogram of Number of PoTs required - Visibility vs No Visibility")

        plot_subtitle = str_glue(
                "Using {treatment_arm} private incentive/visibility."
            )

        if (cutoff_type == "no-cutoff") {
            plot_subtitle = paste0(
                plot_subtitle,
                " Extrapolating past observed support."
            )
        } else {
            plot_subtitle = paste0(
                plot_subtitle,
                " Only assigning PoTs within observed support."

            )
        }

        plot_df = data %>%
            filter(n_pots != 0) %>%
            filter(treatment == treatment_arm) %>%
            mutate(rep_type = if_else(
                rep_type == "no-rep",
                "No Visibility", 
                "Visibility"
            ), 
            rep_type = factor(rep_type, levels = c("Visibility", "No Visibility"))
            )
            
            
            
        p = plot_df %>%
            ggplot(aes(
                x = n_pots, 
                fill = rep_type
            )) +
            geom_histogram(
                bins = 60,
                colour = "black"
            ) +
            theme_minimal() +
            theme(
                legend.position = "bottom"
            ) +
            labs(
                fill = element_blank(),
                x = "Number of PoTs", 
                y = "N Posterior Draws", 
                title = plot_title, 
                subtitle = plot_subtitle,
                caption = "Histogram shows 200 posterior draws. Decision maker targets takeup of 33% (equivalent to takeup in control arm)."
                ) +
            scale_fill_canva(
                palette = canva_palette_vibrant
            )
        return(p)
    }

    p = plot_post_estimates(
        draw_data,
        "bracelet",
        cutoff_type
    )
    p

    treatment_arms = c("bracelet", "control")

    p = map(
        treatment_arms,
        ~plot_post_estimates(
            data = draw_data,
            treatment_arm = .x,
            cutoff_type = cutoff_type
        )
    )
    iwalk(
        p,
        ~ggsave(
            plot = .x,
            filename = file.path(
                script_options$output_path,
                str_glue(
                    "posterior-pots-{cutoff_type}-{treatment_arms[.y]}.png"
                )
            ),
            width = 8,
            height = 6,
            dpi = 500
        )
    )
}
