#!/usr/bin/Rscript

script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        optimal_allocation.R  [options] [--min-cost | --max-takeup]

        Options:
          --min-cost  Flag to minimise programme cost for a given takeup level.
          --max-takeup  Flag to maximise vaccine takeup for a given cost level.
          --target-constraint=<target-constraint>  Amount of takeup/budget constraint depending on [--min-cost|--max-takeup]  [default: 0.5]
"),
  args = if (interactive()) "--min-cost --target-constraint=0.5 " else commandArgs(trailingOnly = TRUE)
) 

library(tidyverse)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)

numeric_options = c(
  "target_constraint"
)

script_options = script_options %>%
  modify_at(numeric_options, as.numeric)


#' Create a MIP Model
#' 
#' @param sim_data list with n and m (n villages and PoTs). Also village_locations
#' df and treatment_locations df
#' @param demand_function function that takes as arguments: i, j, village_locations, pot_locations
#' i is village index, j is PoT index and village_locations, pot_locations are dfs indexed by i,j
#' @param optim_type Whether to solve takeup maximisation given budget or min cost given takeup target. 
#' Accepts: `min_cost` or `max_takeup`
#' @target_constraint What to put in the budget constraint 
define_baseline_MIPModel = function(sim_data, 
                                    demand_function, 
                                    optim_type, 
                                    target_constraint) {
    n = sim_data$n  # N villages
    m = sim_data$m  # M points of treatment
    village_locations = sim_data$village_locations # Village location df
    pot_locations = sim_data$pot_locations # PoT location df

    model = MIPModel() %>%
        # 1 iff village i gets assigned to PoT j
        add_variable(x[i, j], i = 1:n, j = 1:m, type = "binary") %>%
        # 1 iff PoT j is used
        add_variable(y[j], j = 1:m, type = "binary") %>%
        # every village needs to be assigned to a PoT
        add_constraint(sum_over(x[i, j], j = 1:m) == 1, i = 1:n) %>% 
        # if a village is assigned to a PoT, then this PoT must be online
        add_constraint(x[i,j] <= y[j], i = 1:n, j = 1:m) 

    if (optim_type == "min_cost") {
      model = model %>%
          # Takeup must be at least 
          add_constraint(sum_over(
            x[i, j]*demand_function(
              i,
              j,
              village_locations = village_locations,
              pot_locations = pot_locations
            ), i = 1:n, j = 1:m) >=  target_constraint*n)  %>%
          set_objective(
            sum_over(y[j], j = 1:m), "min"
          )
    }

    if (optim_type == "max_takeup") {
      model = model %>%
        # given budget of number PoTs
        add_constraint(
          sum_over(
            y[j],
            j = 1:m
          ) <= target_constraint
        ) %>%
        # maximize the takeup
        set_objective(
                sum_over(
                    x[i,j] * demand_function(
                        i,
                        j,
                        village_locations = village_locations,
                        pot_locations = pot_locations
                    ), i = 1:n, j = 1:m), "max")
    }

    return(model)
}









