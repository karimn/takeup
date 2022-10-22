#!/usr/bin/Rscript

script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        optimal_allocation.R  [options] [--min-cost | --max-takeup]

        Options:
          --dry-run  Run on simulated data
          --min-cost  Flag to minimise programme cost for a given takeup level.
          --max-takeup  Flag to maximise vaccine takeup for a given cost level.
          --target-constraint=<target-constraint>  Amount of takeup/budget constraint depending on [--min-cost|--max-takeup]  [default: 0.5]
          --input-path=<input-path>  Path where input data is stored.
          --output-path=<output-path>  Path where output should be saved.
          --output-filename=<output-filename>  Output filename.
          --village-input-filename=<village-input-filename>  Index and location of each village - a csv path.
          --pot-input-filename=<pot-input-filename>  Index and location of each PoT - a csv path. 
          --demand-input-filename=<demand-input-filename>  Estimated demand for every village i, PoT j pair.  
"),
  args = if (interactive()) "--min-cost 
                             --target-constraint=0.3
                             --dry-run " else commandArgs(trailingOnly = TRUE)
) 

library(tidyverse)
library(data.table)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)

numeric_options = c(
  "target_constraint"
)

script_options = script_options %>%
  modify_at(numeric_options, as.numeric)

optim_type = if_else(script_options$min_cost, "min_cost", "max_takeup")

#' Create a MIP Model
#' 
#' @param data list with n and m (n villages and PoTs). Also village_locations
#' df and treatment_locations df
#' @param demand_function function that takes as arguments: i, j, village_locations, pot_locations
#' i is village index, j is PoT index and village_locations, pot_locations are dfs indexed by i,j
#' @param optim_type Whether to solve takeup maximisation given budget or min cost given takeup target. 
#' Accepts: `min_cost` or `max_takeup`
#' @target_constraint What to put in the budget constraint 
define_baseline_MIPModel = function(data, 
                                    demand_data, 
                                    optim_type, 
                                    target_constraint) {
    n = data$n  # N villages
    m = data$m  # M points of treatment
    village_locations = data$village_locations # Village location df
    pot_locations = data$pot_locations # PoT location df

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
            x[i, j]*demand_data[village_i == i & pot_j == j, demand], i = 1:n, j = 1:m) >=  target_constraint*n)  %>%
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
                    x[i,j] * demand_data[village_i == i & pot_j == j, demand], 
                    i = 1:n, j = 1:m), "max")
    }

    return(model)
}





if (script_options$dry_run) {
  generate_data = function(..., 
                           seed, 
                           n_villages = 50, 
                           m_pots = 30){
      set.seed(seed)
      grid_size <- 10000
      n <- n_villages
      village_locations <- data.frame(
      id = 1:n,
      x = round(runif(n) * grid_size),
      y = round(runif(n) * grid_size)
      )

      m <- m_pots
      pot_locations <- data.frame(
      id = 1:m,
      x = round(runif(m) * grid_size),
      y = round(runif(m) * grid_size)
      )
      fixedcost <- abs(round(rnorm(m, mean = grid_size/20, sd = grid_size)))
      return(
          lst(
              grid_size,
              n,
              village_locations,
              m,
              pot_locations,
              fixedcost
          )
      )
  }
  data = generate_data(n_villages = 100, m_pots = 25, seed = 5)
  sim_demand_function = function(i, j, village_locations, pot_locations){
    subsidy = 0.1
    logit = function(x){
        y = exp(x)/(1 + exp(x))
        return(y)
    }
    village <- village_locations[i, ]
    pot <- pot_locations[j, ]
    distance_away = sqrt((village$x - pot$x)^2 + (village$y - pot$y)^2)
    distance_away = distance_away*(1 - subsidy)
    max_distance = sqrt(2)*10000 # gridsize hardcoded booo
    unit_interval_distance = distance_away/max_distance - 0.5
    demand = 1 - logit(5*unit_interval_distance)
    return(demand)
  }
  i_j_grid = expand.grid(i = 1:data$n, j = 1:data$m)
  demand_data = map2_dfr(
    i_j_grid$i,
    i_j_grid$j,
      ~tibble(
        demand = sim_demand_function(i = .x, 
                                    j = .y, 
                                    village_locations = data$village_locations,
                                    pot_locations = data$pot_locations),
        village_i = .x,
        pot_j = .y 
     )
  ) %>%
    as.data.table()
} else {
  village_input_path = file.path(script_options$input_path, 
                                 script_options$village_input_filename)
  pot_input_path = file.path(script_options$input_path, 
                             script_options$pot_input_filename)
  demand_input_path = file.path(script_options$input_path, 
                                script_options$demand_input_filename)
              
  village_data = read_csv(script_options$village_input_path)
  pot_data = read_csv(script_options$pot_input_path)
  demand_data = read_csv(script_options$demand_input_path) %>%
    as.data.table()
  
  n = nrow(village_data)
  m = nrow(pot_data)

  data = list(
    n = n,
    m = m,
    pot_locations = pot_data,
    village_locations = village_data
  ) 
}


model = define_baseline_MIPModel(
  data,
  demand_data,
  optim_type = optim_type,
  target_constraint = script_options$target_constraint
)



fit_model = solve_model(
  model,
  with_ROI(solver = "glpk", verbose = TRUE)
)




clean_output = function(model_fit, data, demand_data){
    matching = model_fit %>%
        get_solution(x[i,j]) %>%
        filter(value > .9) %>%  
        select(i, j) %>%
        as_tibble()
    tidy_output = matching %>%
        inner_join(data$village_locations %>% 
                        rename(village_x = x, village_y = y), by = c("i" = "id")) %>% 
        inner_join(data$pot_locations %>%
                        rename(pot_x = x, pot_y = y), by = c("j" = "id")) 
    tidy_output = left_join(
      tidy_output,
      demand_data,
      by = c("i" = "village_i", "j" = "pot_j")
    )
    return(tidy_output)
}


tidy_output = map(
    list(fit_model),
    clean_output,
    data = data,
    demand_data = demand_data
)

output_path = file.path(script_options$output_filepath, script_options$output_filename)


write_csv(
  tidy_output,
  output_path
)