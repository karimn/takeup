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
"),
  args = if (interactive()) "
                             --min-cost 
                             --target-constraint=0.4
                             --input-path=optim/data
                             --optim-input-a-filename=init-approx-optimal-allocation.csv
                             --village-input-filename=school-df.csv
                             --pot-input-filename=school-df.csv
                             --demand-input-a-filename=approx-struct-demand.csv
                             --output-path=optim
                             --output-basename=init-approx
                             --map-plot
                             
                             " else commandArgs(trailingOnly = TRUE)
) 


                            #  --optim-input-b-filename=dry-run-subsidy-0-optimal-allocation.csv
                            #  --demand-input-b-filename=dry-run-subsidy-0.2-demand-data.csv
                            #  --comp-demand
library(tidyverse)
library(data.table)

numeric_options = c(
  "target_constraint"
)

script_options = script_options %>%
  modify_at(numeric_options, as.numeric)

optim_type = if_else(script_options$min_cost, "min_cost", "max_takeup")

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
optim_a_data = read_csv(optim_input_a_filepath)
optim_b_data = read_csv(optim_input_b_filepath)
demand_a_data = read_csv(demand_input_a_path) %>%
    as.data.table()
demand_b_data = read_csv(demand_input_b_path) %>%
    as.data.table()

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

optimal_allocation_plot_path = file.path(
    script_options$output_path,
    str_glue(
        "{script_options$output_basename}-optimal-allocation-plot.png"
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
    village_locations = st_as_sf(village_data, coords = c("lon", "lat"), crs = wgs.84), 
    optim_data = optim_a_data
) 

assigned_pots = unique(optim_a_data$j)

data$pot_locations$assigned_pot = data$pot_locations$id %in% assigned_pots


library(scales)
#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(2)


if (script_options$map_plot){

    data$village_locations %>%
        ggplot() +
        # geom_sf()
        geom_sf(
            data = data$pot_locations,
            colour = hex[1],
            shape = 17,
            size = 4, 
            alpha = 1) +
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





optimal_pot_plot = data$village_locations %>%
    ggplot() +
    # geom_sf(alpha = 0.5) +
    geom_sf(
        data = data$pot_locations %>% filter(assigned_pot == FALSE),
        color = hex[1],
        shape = 17,
        size = 4, 
        alpha = 0.3) +
    geom_sf(
        data = data$pot_locations %>% filter(assigned_pot == TRUE),
        color = hex[2],
        shape = 17,
        size = 4, 
        alpha = 1) +
    theme_bw() +
    labs(
        title = "Optimal PoT Allocation Problem"
    )  +
    geom_segment(
        data = data$optim_data, 
        aes(x = village_lon, 
            y = village_lat, 
            xend = pot_lon, 
            yend = pot_lat),
        alpha = 0.1)  +
    labs(x = "", y = "")  +
    scale_color_manual(
        values = hex,
        labels = c("PoT Used", "PoT Unused")
    )


## If we provide comparison data
if (!is.null(optim_b_data)) {

    ms_pots = optim_b_data$j
    dropped_pots = setdiff(ms_pots, assigned_pots)

    optimal_pot_plot = optimal_pot_plot +
        geom_point(
            data = data$pot_locations %>% filter(id %in% dropped_pots),
            color = "red",
            shape = 17,
            size = 4, 
            alpha = 1)  +
    labs(
        title = "Optimal PoT Allocation Problem",
        subtitle = "Black dots indicate villages. Triangles indicate potential clinic locations. 
Extraneous PoTs in red"
    )  

}

ggsave(
    plot = optimal_pot_plot,
    filename = optimal_allocation_plot_path,
    width = 8,
    height = 6,
    dpi = 500
)



if (script_options$comp_demand) {

    comp_demand_data = inner_join(
        demand_a_data %>% rename(demand_a = demand),
        demand_b_data %>% rename(demand_b = demand),
        by = c("village_i", "pot_j")
    )


    comp_demand_plot = comp_demand_data %>%
        ggplot(aes(
            y = demand_b, 
            x = demand_a
        )) +
        geom_point(alpha = 0.2) +
        geom_abline(linetype = "longdash") +
        theme_bw() +
        labs(
            y = "Takeup Probability Including Signal Benefit", 
            x = "Takeup Probability Naive Social Planner"
        )


    total_overshoot = optim_a_data %>%
        left_join(
            demand_b_data %>% rename(demand_b = demand),
            by = c("j" = "pot_j", "i" = "village_i")) %>%
        mutate(
            diff = abs(demand - demand_b)
        ) %>%
        summarise(total_overshoot = sum(diff)) %>%
        pull()


    print(str_glue("Total overshoot: {total_overshoot}"))

    both_demand_df = data$optim_data %>%
                left_join(
                    demand_b_data %>% rename(demand_b = demand),
                    by = c("j" = "pot_j", "i" = "village_i"))  %>%
                mutate(
                    diff = demand_b - demand, 
                    pct_diff = 100*diff/demand
                )


    data$village_locations %>%
        ggplot(aes(
            x = x,
            y = y
        )) +
        geom_point(alpha = 0.5) +
        geom_point(
            data = data$pot_locations %>% filter(assigned_pot == FALSE),
            color = hex[1],
            shape = 17,
            size = 4, 
            alpha = 0.3) +
        geom_point(
            data = data$pot_locations %>% filter(assigned_pot == TRUE),
            color = hex[2],
            shape = 17,
            size = 4, 
            alpha = 1) +
        theme_bw() +
        labs(
            title = "Optimal PoT Allocation Problem",
            subtitle = "Black dots indicate villages. Triangles indicate potential clinic locations."
        )  +
        geom_segment(
            data = both_demand_df,
            aes(x = village_x, 
                y = village_y, 
                xend = pot_x, 
                yend = pot_y, 
                colour = pct_diff),
            alpha = 0.5) +
            scale_colour_viridis_c(option = "magma", direction = -1) +
        geom_point(
            data = both_demand_df,
            aes(
                x = village_x, 
                y = village_y, 
                colour = pct_diff)
        ) +
        theme(legend.position = "bottom") +
        labs(x = "", y = "") + 
        theme(
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
        labs(
            colour = "Percentage Takeup Underestimated",
            caption = str_glue(
                "Ignoring reputational concerns leads to overshooting target takeup by {round(total_overshoot, 2)} percentage points"
            )
        )

    ggsave(optimal_allocation_demand_comp_plot_path,
    width = 8,
    height = 6,
    dpi = 500)


}

