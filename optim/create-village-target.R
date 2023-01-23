#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        create-village-target.R <model-fit> [options]

        Options:
          --input-path=<input-path>  Path where input data is stored.
          --output-path=<output-path>  Path where output should be saved.
          --output-basename=<output-basename>  Output basename.
"),
  args = if (interactive()) "
                            pred-demand-dist-fit71-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS.csv \
                            --input-path=optim/data \
                            --output-path=optim/data \
                            --output-basename=target-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS \
                             
                             " else commandArgs(trailingOnly = TRUE)
) 


library(tidyverse)



demand_df = read_csv(
    file.path(
        script_options$input_path,
        script_options$model_fit
    )
)

target_demand_df = demand_df %>%
    group_by(village_i) %>%
    # adjustment for 3 villages which don't seem to have anyone close by
    mutate(all_far = all(dist > 2500 )) %>%
    ungroup() %>%
    filter(dist < 2500 | all_far == TRUE) %>%
    group_by(
        village_i
    ) %>%
    summarise(
        demand = mean(demand)
    ) %>%
    arrange(village_i)

target_demand_df %>%
    write_csv(
        file.path(
            script_options$output_path, 
            paste0(script_options$output_basename, ".csv")
        )
    )

