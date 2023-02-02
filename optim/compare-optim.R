#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        compare-optim.R  [options] 

        Options:
        --input-path=<input-path>  The input path [default: {file.path('optim', 'data', 'agg-log-full')}]
        --output-path=<output-path>  The output path [default: {file.path('temp-data')}]
        --many-pots
"),
  args = if (interactive()) "
        --input-path=~/projects/takeup/optim/data/agg
                             
                             " else commandArgs(trailingOnly = TRUE)
) 

library(tidyverse)
library(sf)

input_path = script_options$input_path 
output_path = script_options$output_path 
dir.create(output_path)
many_pots = script_options$many_pots

oa_files = fs::dir_ls(
    input_path,
    regexp = "rds"
)

models_we_want = "STRUCTURAL_LINEAR_U_SHOCKS"

oa_df = map_dfr(
    oa_files, 
    read_rds
) %>% as_tibble() %>%
    mutate(
        file = oa_files
    ) %>%
    mutate(
        cutoff_type = if_else(str_detect(file, "no-cutoff"), "no_cutoff", "cutoff")
    ) 
treatments = c(
    "bracelet",
    "calendar",
    "ink", 
    "control"
)


oa_df = oa_df %>%
    mutate(
        private_benefit_z = factor(private_benefit_z, levels = treatments), 
        visibility_z = factor(visibility_z, levels = treatments), 
    )

subset_oa_df = oa_df %>%
    filter(model == "STRUCTURAL_LINEAR_U_SHOCKS")  %>%
    filter(
        model %in% models_we_want
    ) %>%
    filter(
        cutoff_type  == "cutoff"
    ) 


model_df = subset_oa_df %>%
    unnest(model_output) 

if (many_pots == TRUE) {
    target_optim = -174
} else {
    target_optim = -171.788

}

summ_optim_df = model_df %>%
    group_by(
        private_benefit_z,
        visibility_z, 
        model, 
        cutoff_type
    ) %>%
    summarise(
        util = sum(log(demand)),
        mean_demand = mean(demand), 
        min_demand = min(demand), 
        n_pot = n_distinct(j)
    ) %>%
    filter(
        model %in% models_we_want
    ) %>%
    filter(
        cutoff_type  == "cutoff"
    ) %>%
    mutate(
        target_optim = target_optim
    ) %>%
    mutate(
        overshoot = 100*(util/target_optim - 1)
    ) %>%
    select(-model, -cutoff_type) %>%
    ungroup() %>%
    select(-model)

summ_optim_df %>%
    filter(private_benefit_z == "control")

summ_optim_df %>%
    filter(visibility_z == "control")


summ_optim_df %>% 
    write_csv(
        file.path(
            output_path,
            "full-summ-optim.csv"
        ))


demand_df = subset_oa_df %>%
    unnest(demand_data)

demand_df %>%
    filter(dist < 3500) %>%
    filter(private_benefit_z == "control") %>%
    ggplot(aes(
        x =  dist,
        y = demand, 
        colour = visibility_z
    )) +
    geom_line()  +
    geom_point()  +
    theme_minimal() +
    theme(legend.position = "bottom") 

ggsave(
    file.path(
        output_path,
        "optim-demand-df.png"
    ),
    width = 8,
    height = 6,
    dpi = 500
)


demand_df %>%
    filter(
        visibility_z == "control"
    ) %>%
    filter(dist < 3500) %>%
    ggplot(aes(
         x = dist, 
         y = demand, 
         colour = private_benefit_z
    )) +
    geom_point() +
    geom_hline(yintercept =  0.33) +
    geom_vline(
        xintercept = 1000
    )


library(ggridges)
subset_oa_df %>%
    filter(private_benefit_z == "control") %>%
    unnest(model_output) %>%
    ggplot(aes(
        x = dist, 
        fill = visibility_z, 
        y = visibility_z
    )) +
    geom_density_ridges(alpha = 0.3) +
    theme_ridges() +
    theme(legend.position = "bottom")

subset_oa_df %>%
    filter(private_benefit_z == "control") %>%
    mutate(visibility_z = fct_relabel(visibility_z, str_to_title)) %>%
    unnest(model_output) %>%
    ggplot(aes(
        x = dist, 
        fill = visibility_z
    )) +
    geom_density(alpha = 0.3) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        fill = "",
        x = "Distance (m)", 
        y = "Density"
    ) +
    ggthemes::scale_fill_canva(palette = "Primary colors with a vibrant twist")

ggsave(
    file.path(
        output_path,
        "optim-walk-density-b-control.png"
    ),
    width = 8,
    height = 6,
    dpi = 500
)













