#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        postprocess_allocation.R  [options] [--min-cost | --max-takeup]

        Options:
          --min-cost  Flag to minimise programme cost for a given takeup level.
          --max-takeup  Flag to maximise vaccine takeup for a given cost level.
          --target-constraint=<target-constraint>  Amount of takeup/budget constraint depending on [--min-cost|--max-takeup]  [default: 0.5]
          --input-path=<input-path>  Path where input data is stored.
          --optim-input-a-filename=<optim-input-a-filename>  Optim input a filename.
          --optim-input-b-filename=<optim-input-b-filename>   Optim input b filename
          --demand-input-a-filename=<demand-input-a-filename>  Estimated demand a for every village i, PoT j pair.  
          --demand-input-b-filename=<demand-input-b-filename>  Estimated demand b for every village i, PoT j pair.  
          --village-input-filename=<village-input-filename>  Index and location of each village - a csv path.
          --pot-input-filename=<pot-input-filename>  Index and location of each PoT - a csv path. 
          --output-path=<output-path>  Path where output should be saved.
          --output-basename=<output-basename>  Output basename.
          --comp-demand  Whether to compare demand under signaling vs naive
          --map-plot  Plot map without adding lines for closest villages and active PoTs 
          --posterior-median  If posterior median over takeup demand is used.
"),
  args = if (interactive()) "
                             --min-cost 
                             --target-constraint=0.33
                             --input-path=optim/data
                             --optim-input-a-filename=structural-post-draws-optimal-allocation.rds
                             --village-input-filename=village-df.csv
                             --pot-input-filename=pot-df.csv
                             --demand-input-a-filename=pred_demand_dist_fit66_STRUCTURAL_LINEAR_U_SHOCKS.csv
                             --output-path=optim
                             --output-basename=structural-post-draws-test
                             --map-plot
                             
                             " else commandArgs(trailingOnly = TRUE)
) 


                            #  --optim-input-b-filename=dry-run-subsidy-0-optimal-allocation.csv
                            #  --demand-input-b-filename=dry-run-subsidy-0.2-demand-data.csv
                            #  --comp-demand
library(tidyverse)
library(sf)
library(data.table)

numeric_options = c(
  "target_constraint"
)


script_options = script_options %>%
  modify_at(numeric_options, as.numeric)

optim_type = if_else(script_options$min_cost, "min_cost", "max_takeup")
stat_type = if_else(script_options$posterior_median, "median", "post-draws")

## Input Paths
optim_input_a_filepath = file.path(script_options$input_path, 
                                 script_options$optim_input_a_filename)
optim_input_b_filepath = file.path(script_options$input_path, 
                                 script_options$optim_input_b_filename)
demand_input_a_path = file.path(script_options$input_path, 
                            script_options$demand_input_a_filename)
demand_input_b_path = file.path(script_options$input_path, 
                            script_options$demand_input_b_filename)

village_input_path = file.path(script_options$input_path, 
                                script_options$village_input_filename)
pot_input_path = file.path(script_options$input_path, 
                            script_options$pot_input_filename)
## Input Data
optimal_df = read_rds(optim_input_a_filepath)
# optim_b_data = read_csv(optim_input_b_filepath)
demand_a_data = read_csv(demand_input_a_path) %>%
    as.data.table()
# demand_b_data = read_csv(demand_input_b_path) %>%
#     as.data.table()

village_data = read_csv(village_input_path)
pot_data = read_csv(pot_input_path)

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

optimal_df = optimal_df %>%
    mutate(allocation_infeasible = map(model_output, ~any(str_detect(colnames(.x), "fail")))) %>%
    filter(allocation_infeasible == FALSE)

if (stat_type == "median") {


    plot_optimal_allocation = function(village_data,
                                    pot_data,
                                    optimal_data){

        assigned_pots = unique(optimal_data$j)
        n_pots_used =  length(assigned_pots)
        pot_data$assigned_pot = pot_data$id %in% assigned_pots

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
            labs(x = "", y = "", subtitle = str_glue("
            Takeup target: {round(100*script_options$target_constraint)}%. Assigned PoTs: {n_pots_used}"))  +
            scale_color_manual(
                values = hex,
                labels = c("PoT Used", "PoT Unused")
            ) 
        return(optimal_pot_plot)
    }

    plot_optimal_allocation(
        village_data = data$village_locations,
        pot_data = data$pot_locations,
        optimal_data = optimal_df %>%
            filter(treatment == "bracelet" & str_detect(model, "RED")) %>%
            unnest(model_output)
    )


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
                    "{script_options$output_basename}-{optimal_df[.y, 'treatment']}-{optimal_df[.y, 'model']}-optimal-allocation-plot.png"
                )
            ),
        width = 10,
        height = 10,
        dpi = 500
        )
    )

}

if (stat_type == "post-draws") {
    hist_pots = optimal_df %>%
        mutate(n_pot_used = map_dbl(model_output, ~n_distinct(.x$j)))  %>%
        as_tibble() %>%
        ggplot(aes(
            x = n_pot_used, 
            fill = treatment
        )) +
        geom_histogram() +
        facet_wrap(~model) +
        theme_minimal()

    ggsave(
        file.path(
            script_options$output_path,
            str_glue(
                "{script_options$output_basename}-posterior-PoTs-used.png"
            )
        ),
        width = 10,
        height = 10,
        dpi = 500
    )
}
    



# summ_optimal_df = optimal_df %>%
#     group_by(
#         i,
#         j,
#         treatment
#     ) %>%
#     summarise( 
#         n_links = n(),
#         village_lon = unique(village_lon), 
#         village_lat = unique(village_lat), 
#         pot_lon = unique(pot_lon), 
#         pot_lat = unique(pot_lat)
#     ) %>%
#     mutate(
#         pct_links = n_links / 50
#         )

# pot_summ_optimal_df = optimal_df %>%
#     select(j, draw, treatment) %>%
#     unique() %>%
#     group_by(j, treatment) %>%
#     summarise(
#         n_used = n(),
#         across(contains("_lon"), unique), 
#         across(contains("_lat"), unique)
#         ) %>%
#     mutate(pct_used = n_used / 50)

# optimal_df %>%
#     select(solver_status)

# optimal_df  %>%
#     filter(is.na(solver_status)) %>%
#     group_by(draw, treatment) %>%
#     summarise(
#         n_pots = n_distinct(j)
#     ) %>%
#     ggplot(aes( 
#         x = n_pots, 
#         fill = treatment
#     )) + 
#     geom_histogram()



# data$pot_locations %>%
#     left_join(
#         pot_summ_optimal_df, 
#         by = c("id" = "j")
#     ) %>%
#     filter(!is.na(treatment)) %>%
#     ggplot() +
#     geom_sf(
#         aes(color = pct_used)
#     ) +
#     facet_wrap(~treatment)

#  data$pot_locations %>%
#     ggplot() +
#     geom_sf(
#         aes(color = pct_used)
#         shape = 17,
#         size = 4, 
#         alpha = 0.2) 

#  data$village_locations %>%
#     ggplot() +
#     geom_sf(alpha = 1) +
#     geom_sf(
#         data = data$pot_locations %>% filter(assigned_pot == FALSE),
#         color = hex[1],
#         shape = 17,
#         size = 4, 
#         alpha = 0.2) +
#     geom_sf(
#         data = data$pot_locations %>% filter(assigned_pot == TRUE),
#         color = hex[2],
#         shape = 17,
#         size = 4, 
#         alpha = 1) +
#     theme_bw() +
#     geom_segment(
#         data = summ_optimal_df, 
#         aes(x = village_lon, 
#             y = village_lat, 
#             xend = pot_lon, 
#             yend = pot_lat, 
#             alpha = pct_links))  +
#     labs(x = "", y = "", subtitle = str_glue("
#     Takeup target: {round(100*script_options$target_constraint)}%. Assigned PoTs: {n_pots_used}"))  +
#     scale_color_manual(
#         values = hex,
#         labels = c("PoT Used", "PoT Unused")
#     ) 

# optimal_pot_plot = data$village_locations %>%
#     ggplot() +
#     geom_sf(alpha = 1) +
#     geom_sf(
#         data = data$pot_locations %>% filter(assigned_pot == FALSE),
#         color = hex[1],
#         shape = 17,
#         size = 4, 
#         alpha = 0.2) +
#     geom_sf(
#         data = data$pot_locations %>% filter(assigned_pot == TRUE),
#         color = hex[2],
#         shape = 17,
#         size = 4, 
#         alpha = 1) +
#     theme_bw() +
#     geom_segment(
#         data = summ_optimal_df, 
#         aes(x = village_lon, 
#             y = village_lat, 
#             xend = pot_lon, 
#             yend = pot_lat, 
#             alpha = pct_links))  +
#     labs(x = "", y = "", subtitle = str_glue("
#     Takeup target: {round(100*script_options$target_constraint)}%. Assigned PoTs: {n_pots_used}"))  +
#     scale_color_manual(
#         values = hex,
#         labels = c("PoT Used", "PoT Unused")
#     ) 

# optimal_pot_plot

# ## If we provide comparison data

#     ms_pots = optim_b_data$j
#     dropped_pots = setdiff(ms_pots, assigned_pots)

#     optimal_pot_plot = optimal_pot_plot +
#         geom_point(
#             data = data$pot_locations %>% filter(id %in% dropped_pots),
#             color = "red",
#             shape = 17,
#             size = 4, 
#             alpha = 1)  +
#     labs(
#         title = str_glue("Optimal PoT Allocation Problem"),
#         caption = "Black dots indicate villages. Triangles indicate potential clinic locations. 
# Extraneous PoTs in red"
#     )  

# ggsave(
#     plot = optimal_pot_plot,
#     filename = optimal_allocation_plot_path,
#     width = 8,
#     height = 6,
#     dpi = 500
# )



# if (script_options$comp_demand) {

#     comp_demand_data = inner_join(
#         demand_a_data %>% rename(demand_a = demand),
#         demand_b_data %>% rename(demand_b = demand),
#         by = c("village_i", "pot_j")
#     )


#     comp_demand_plot = comp_demand_data %>%
#         ggplot(aes(
#             y = demand_b, 
#             x = demand_a
#         )) +
#         geom_point(alpha = 0.2) +
#         geom_abline(linetype = "longdash") +
#         theme_bw() +
#         labs(
#             y = "Takeup Probability Including Signal Benefit", 
#             x = "Takeup Probability Naive Social Planner"
#         )


#     total_overshoot = optim_a_data %>%
#         left_join(
#             demand_b_data %>% rename(demand_b = demand),
#             by = c("j" = "pot_j", "i" = "village_i")) %>%
#         mutate(
#             diff = abs(demand - demand_b)
#         ) %>%
#         summarise(total_overshoot = sum(diff)) %>%
#         pull()


#     print(str_glue("Total overshoot: {total_overshoot}"))

#     both_demand_df = data$optim_data %>%
#                 left_join(
#                     demand_b_data %>% rename(demand_b = demand),
#                     by = c("j" = "pot_j", "i" = "village_i"))  %>%
#                 mutate(
#                     diff = demand_b - demand, 
#                     pct_diff = 100*diff/demand
#                 )


#     data$village_locations %>%
#         ggplot(aes(
#             x = x,
#             y = y
#         )) +
#         geom_point(alpha = 0.5) +
#         geom_point(
#             data = data$pot_locations %>% filter(assigned_pot == FALSE),
#             color = hex[1],
#             shape = 17,
#             size = 4, 
#             alpha = 0.3) +
#         geom_point(
#             data = data$pot_locations %>% filter(assigned_pot == TRUE),
#             color = hex[2],
#             shape = 17,
#             size = 4, 
#             alpha = 1) +
#         theme_bw() +
#         labs(
#             title = "Optimal PoT Allocation Problem",
#             subtitle = "Black dots indicate villages. Triangles indicate potential clinic locations."
#         )  +
#         geom_segment(
#             data = both_demand_df,
#             aes(x = village_x, 
#                 y = village_y, 
#                 xend = pot_x, 
#                 yend = pot_y, 
#                 colour = pct_diff),
#             alpha = 0.5) +
#             scale_colour_viridis_c(option = "magma", direction = -1) +
#         geom_point(
#             data = both_demand_df,
#             aes(
#                 x = village_x, 
#                 y = village_y, 
#                 colour = pct_diff)
#         ) +
#         theme(legend.position = "bottom") +
#         labs(x = "", y = "") + 
#         theme(
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#         labs(
#             colour = "Percentage Takeup Underestimated",
#             caption = str_glue(
#                 "Ignoring reputational concerns leads to overshooting target takeup by {round(total_overshoot, 2)} percentage points"
#             )
#         )

#     ggsave(optimal_allocation_demand_comp_plot_path,
#     width = 8,
#     height = 6,
#     dpi = 500)


# }

