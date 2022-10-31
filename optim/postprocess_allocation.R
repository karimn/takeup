#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        postprocess_allocation.R  [options] [--min-cost | --max-takeup]

        Options:
          --min-cost  Flag to minimise programme cost for a given takeup level.
          --max-takeup  Flag to maximise vaccine takeup for a given cost level.
          --target-constraint=<target-constraint>  Amount of takeup/budget constraint depending on [--min-cost|--max-takeup]  [default: 0.5]
          --input-path=<input-path>  Path where input data is stored.
          --optim-input-filename=<input-filename>  Input filename.
          --village-input-filename=<village-input-filename>  Index and location of each village - a csv path.
          --pot-input-filename=<pot-input-filename>  Index and location of each PoT - a csv path. 
          --demand-input-filename=<demand-input-filename>  Estimated demand for every village i, PoT j pair.  
          --true-demand-input-filename=<true-demand-input-filename>  Estimated demand for every village i, PoT j pair.  
          --output-path=<output-path>  Path where output should be saved.
          --output-filename=<output-filename>  Output filename.
          --misspecified-optim-input-filename=<misspecified-optim-input-filename>  Misspecified optim filename
          --comp-demand  Whether to compare demand under signaling vs naive
"),
  args = if (interactive()) "
                             --min-cost 
                             --target-constraint=0.85
                             --input-path=optim/data
                             --optim-input-filename=dry-run-subsidy-0.2-optimal-allocation.csv
                             --misspecified-optim-input-filename=dry-run-subsidy-0-optimal-allocation.csv
                             --village-input-filename=dry-run-subsidy-0.2-village-locations.csv
                             --pot-input-filename=dry-run-subsidy-0.2-pot-locations.csv
                             --demand-input-filename=dry-run-subsidy-0.2-demand-data.csv
                             --true-demand-input-filename=dry-run-subsidy-0.2-demand-data.csv
                             --output-path=optim
                             --output-filename=problem-c-optimal-pot-plot.png
                             
                             
                             " else commandArgs(trailingOnly = TRUE)
) 




library(tidyverse)
library(data.table)

numeric_options = c(
  "target_constraint"
)

script_options = script_options %>%
  modify_at(numeric_options, as.numeric)

optim_type = if_else(script_options$min_cost, "min_cost", "max_takeup")

## Input Paths
optim_input_filepath = file.path(script_options$input_path, 
                                 script_options$optim_input_filename)
village_input_path = file.path(script_options$input_path, 
                                script_options$village_input_filename)
pot_input_path = file.path(script_options$input_path, 
                            script_options$pot_input_filename)
demand_input_path = file.path(script_options$input_path, 
                            script_options$demand_input_filename)
true_demand_input_path = file.path(script_options$input_path, 
                            script_options$true_demand_input_filename)
ms_optim_input_filepath = file.path(script_options$input_path, 
                                 script_options$misspecified_optim_input_filename)
## Input Data
optim_data = read_csv(optim_input_filepath)
ms_optim_data = read_csv(ms_optim_input_filepath)
village_data = read_csv(village_input_path)
pot_data = read_csv(pot_input_path)
demand_data = read_csv(demand_input_path) %>%
    as.data.table()
true_demand_data = read_csv(true_demand_input_path) %>%
    as.data.table()
n = nrow(village_data)
m = nrow(pot_data)



data = list(
    n = n,
    m = m,
    pot_locations = pot_data,
    village_locations = village_data, 
    optim_data = optim_data
) 

assigned_pots = unique(optim_data$j)

data$pot_locations$assigned_pot = data$pot_locations$id %in% assigned_pots


library(scales)
#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(2)




optimal_pot_plot = data$village_locations %>%
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
        data = data$optim_data, 
        aes(x = village_x, 
            y = village_y, 
            xend = pot_x, 
            yend = pot_y),
        alpha = 0.1) 

if (!is.null(ms_optim_data)) {

    ms_pots = ms_optim_data$j
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

output_filepath = file.path(
    script_options$output_path,
    script_options$output_filename
)


ggsave(
    plot = optimal_pot_plot,
    filename = output_filepath,
    width = 8,
    height = 6,
    dpi = 500
)





comp_demand_data = inner_join(
    demand_data %>% rename(sp_obs_demand = demand),
    true_demand_data %>% rename(true_demand = demand),
    by = c("village_i", "pot_j")
)

comp_demand_plot = comp_demand_data %>%
    ggplot(aes(
        y = true_demand, 
        x = sp_obs_demand
    )) +
    geom_point(alpha = 0.2) +
    geom_abline(linetype = "longdash") +
    theme_bw() +
    labs(
        y = "Takeup Probability Including Signal Benefit", 
        x = "Takeup Probability Naive Social Planner"
    )




overshoot_df = optim_data %>%
    left_join(
        true_demand_data %>% rename(true_demand = demand),
        by = c("j" = "pot_j", "i" = "village_i")) %>%
    summarise(
        mean_takeup_naive_sp = mean(demand),
        mean_takeup_reality = mean(true_demand)
    ) %>%
    mutate(overshoot_pct = (mean_takeup_reality/mean_takeup_naive_sp - 1)*100)



if (script_options$comp_demand) {

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
            data = data$optim_data %>%
                left_join(
                    true_demand_data %>% rename(true_demand = demand),
                    by = c("j" = "pot_j", "i" = "village_i"))  %>%
                mutate(
                    diff = demand - true_demand, 
                    pct_diff = 100*diff/demand
                ),
            aes(x = village_x, 
                y = village_y, 
                xend = pot_x, 
                yend = pot_y, 
                colour = pct_diff),
            alpha = 0.5) +
            scale_colour_viridis_c(option = "magma") +
        geom_point(
            data = data$optim_data %>%
                left_join(
                    true_demand_data %>% rename(true_demand = demand),
                    by = c("j" = "pot_j", "i" = "village_i"))  %>%
                mutate(
                    diff = demand - true_demand, 
                    pct_diff = 100*diff/demand
                ),
            aes(
                x = village_x, 
                y = village_y, 
                colour = pct_diff)
        ) +
        theme(legend.position = "bottom") +
        labs(
            colour = "Percentage Takeup Misestimated"
        )
    ggsave("optim/problem-b-optimal-pot-plot.png",
    width = 8,
    height = 6,
    dpi = 500)


}


