#!/usr/bin/Rscript

script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        optimal_allocation.R  [options] [--min-cost | --max-takeup]

        Options:
          --dry-run  Run on simulated data
          --dry-run-subsidy=<dry-run-subsidy>  Amount of fake subsidy. 
          --min-cost  Flag to minimise programme cost for a given takeup level.
          --max-takeup  Flag to maximise vaccine takeup for a given cost level.
          --target-constraint=<target-constraint>  Amount of takeup/budget constraint depending on [--min-cost|--max-takeup]  [default: 0.5]
          --input-path=<input-path>  Path where input data is stored.
          --output-path=<output-path>  Path where output should be saved.
          --output-filename=<output-filename>  Output filename.
          --village-input-filename=<village-input-filename>  Index and location of each village - a csv path.
          --pot-input-filename=<pot-input-filename>  Index and location of each PoT - a csv path. 
          --demand-input-filename=<demand-input-filename>  Estimated demand for every village i, PoT j pair.  
          --num-cores=<num-cores>  Number of cores to use in parallel [default: 8]
          --time-limit=<time-limit>  GLPK solver time limit [default: 1000]
          --posterior-median  Use median demand rather than solve across all draws
"),
  args = if (interactive()) "
                             --posterior-median
                             --num-cores=12
                             --min-cost 
                             --target-constraint=0.32
                             --output-path=optim/data
                             --output-filename=structural-no-rep
                             --input-path=optim/data 
                             --village-input-filename=village-df.csv
                             --pot-input-filename=pot-df.csv
                             --demand-input-filename=pred_demand_dist_fit71_no_rep.csv
                             " else commandArgs(trailingOnly = TRUE)
) 
                            #  --dry-run 
                            #  --dry-run-subsidy=0.2

library(tidyverse)
library(data.table)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)
library(sf)

numeric_options = c(
  "target_constraint",
  "dry_run_subsidy",
  "num_cores",
  "time_limit"
)

script_options = script_options %>%
  modify_at(numeric_options, as.numeric)

optim_type = if_else(script_options$min_cost, "min_cost", "max_takeup")
stat_type = if_else(script_options$posterior_median, "median", "post-draws")

#' Create a MIP Model
#' 
#' @param data list with n and m (n villages and PoTs). Also village_locations
#' df and treatment_locations df
#' @param demand_function function that takes as arguments: i, j, village_locations, pot_locations
#' i is village index, j is PoT index and village_locations, pot_locations are dfs indexed by i,j
#' @param optim_type Whether to solve takeup maximisation given budget or min cost given takeup target. 
#' Accepts: `min_cost` or `max_takeup`
#' @target_constraint What to put in the budget constraint 
define_baseline_MIPModel = function(data) {
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
    return(model)
}

#' Add objectives to baseline MIPModel
#'
#' @param data list with n and m (n villages and PoTs). Also village_locations
#' df and treatment_locations df
#' @param demand_function function that takes as arguments: i, j, village_locations, pot_locations
#' i is village index, j is PoT index and village_locations, pot_locations are dfs indexed by i,j
#' @param optim_type Whether to solve takeup maximisation given budget or min cost given takeup target. 
#' Accepts: `min_cost` or `max_takeup`
#' @target_constraint What to put in the budget constraint 
add_MIPModel_objective = function(model, data, demand_data, optim_type, target_constraint) {
    n = data$n  # N villages
    m = data$m  # M points of treatment
    village_locations = data$village_locations # Village location df
    pot_locations = data$pot_locations # PoT location df

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


define_and_solve_model = function(baseline_model,
                                  data,
                                  demand_data,
                                  optim_type,
                                  target_constraint){
  model = baseline_model %>%
    add_MIPModel_objective(
      model =  .,
      data = data,
      demand_data = demand_data,
      optim_type = optim_type,
      target_constraint = target_constraint
    )                            
  fit_model = solve_model(
    model,
    with_ROI(solver = "glpk", verbose = TRUE, control = list(tm_limit = script_options$time_limit))
  )
  status = solver_status(fit_model)
  if (status == "success"){
    match_df = fit_model %>%
        get_solution(x[i,j]) %>%
        filter(value > .9) %>%  
        select(i, j) %>%
        as_tibble()
    tidy_output = clean_output(match_df, data, demand_data) 
    return(tidy_output)
  } else {
    return(tibble(fail = TRUE, solver_status = status))
  }
}

safe_define_and_solve_model = possibly(define_and_solve_model, otherwise = tibble(fail = TRUE))

clean_output = function(match_df, data, demand_data){
    tidy_output = match_df %>%
        inner_join(data$village_locations %>% 
                        rename(village_lon = lon, village_lat = lat), by = c("i" = "id")) %>% 
        inner_join(data$pot_locations %>%
                        rename(pot_lon = lon, pot_lat = lat), by = c("j" = "id")) 
    tidy_output = left_join(
      tidy_output,
      demand_data,
      by = c("i" = "village_i", "j" = "pot_j")
    )

  ## Reassign villages to closest PoT
  # SP doesn't care about distance cost atm so once he hits target takeup, can 
  # assign people to whatever PoT if they won't let him drop down to a smaller # of PoTs
  # therefore, we just reoptimise within allowed PoTs

  #### TEMPORARY FIX #####
  dd = copy(demand_data)
  assigned_pots = unique(tidy_output$j)
  new_tidy_output = dd[pot_j %in% assigned_pots, min_dist_avail := min(dist), village_i][
    dist == min_dist_avail
    ] %>%
    inner_join(data$village_locations %>% 
                    rename(village_lon = lon, village_lat = lat), by = c("village_i" = "id")) %>% 
    inner_join(data$pot_locations %>%
                    rename(pot_lon = lon, pot_lat = lat), by = c("pot_j" = "id"))  %>%
    rename(
      i = village_i, 
      j = pot_j
    ) %>%
    as_tibble() %>%
    select(-min_dist_avail)
  if (nrow(tidy_output) == nrow(new_tidy_output)) {
    tidy_output = new_tidy_output
  } else {
    warning("Temporary fix failed, using original allocation.")
  }
    return(tidy_output)
}


wgs.84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"



if (script_options$dry_run) {
  generate_data = function(..., 
                           seed, 
                           n_villages = 100, 
                           m_pots = 30){
      set.seed(seed)
      grid_size <- 10000
      n <- n_villages
      village_locations <- data.frame(
      id = 1:n,
      lon = round(runif(n) * grid_size),
      lat = round(runif(n) * grid_size)
      )

      m <- m_pots
      pot_locations <- data.frame(
      id = 1:m,
      lon = round(runif(m) * grid_size),
      lat = round(runif(m) * grid_size)
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

  data = generate_data(n_villages = 200, m_pots = 40, seed = 1)
  sim_demand_function = function(i, j, village_locations, pot_locations, subsidy){
    logit = function(x){
        y = exp(x)/(1 + exp(x))
        return(y)
    }
    village <- village_locations[i, ]
    pot <- pot_locations[j, ]
    distance_away = sqrt((village$lon - pot$lon)^2 + (village$lat - pot$lat)^2)
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
                                    pot_locations = data$pot_locations,
                                    subsidy = script_options$dry_run_subsidy),
        village_i = .x,
        pot_j = .y 
     )
  ) %>%
    as.data.table()
  data$village_locations %>%
    write_csv(str_glue("optim/data/dry-run-subsidy-{script_options$dry_run_subsidy}-village-locations.csv"))
    
  data$pot_locations %>%
    write_csv(str_glue("optim/data/dry-run-subsidy-{script_options$dry_run_subsidy}-pot-locations.csv"))
  demand_data %>%
    write_csv(str_glue("optim/data/dry-run-subsidy-{script_options$dry_run_subsidy}-demand-data.csv"))
} else {
  village_input_path = file.path(script_options$input_path, 
                                 script_options$village_input_filename)
  pot_input_path = file.path(script_options$input_path, 
                             script_options$pot_input_filename)
  demand_input_path = file.path(script_options$input_path, 
                                script_options$demand_input_filename)
              
  village_data = read_csv(village_input_path) %>%
    st_as_sf(
      coords = c("lon", "lat"), 
      crs = wgs.84) %>%
    mutate(
        lon = sf::st_coordinates(.)[, 1],
        lat = sf::st_coordinates(.)[, 2])
  pot_data = read_csv(pot_input_path) %>%
    st_as_sf(
      coords = c("lon", "lat"), 
      crs = wgs.84) %>%
    mutate(
        lon = sf::st_coordinates(.)[, 1],
        lat = sf::st_coordinates(.)[, 2])
  demand_data = read_csv(demand_input_path) %>%
    as.data.table() 

  n_max = Inf
  m_max = Inf
  demand_data = demand_data[village_i <= n_max & pot_j <= m_max]
  pot_data = pot_data %>%
    filter(id <= m_max)
  village_data = village_data %>%
    filter(id <= n_max)


  if (script_options$posterior_median) {
    demand_data = demand_data[
      , 
      .(demand = median(demand), dist = unique(dist)), 
      .(village_i,
        pot_j, 
        treatment, 
        model)]
  }
  n = nrow(village_data)
  m = nrow(pot_data)

  data = list(
    n = n,
    m = m,
    pot_locations = pot_data,
    village_locations = village_data
  ) 
}

library(furrr)
plan(multicore, workers = script_options$num_cores)
baseline_model = define_baseline_MIPModel(data)


demand_data = demand_data %>%
  nest(demand_data = -any_of(c("draw", "treatment", "model")))


test_demand_data = demand_data$demand_data[[3]]
full_model = baseline_model %>%
  add_MIPModel_objective(
    data = data,
    demand_data = test_demand_data,
    optim_type = optim_type,
    target_constraint = script_options$target_constraint 
  )


##
get_model_stuff = function(model) {
  modelConstraints = extract_constraints(model)
  variables = variable_types(model)
  variables = as.character(variables)
  variables[grepl("binary",variables)]="B"
  variables[grepl("continuous",variables)]="C"
  variables[grepl("integer",variables)]="I"
  modelConstraints$direction[grepl("==",modelConstraints$direction)]="="
  gurobiModel = list()
  gurobiModel$A = modelConstraints$matrix
  gurobiModel$obj = objective_function(model)$vector
  gurobiModel$sense = modelConstraints$direction
  gurobiModel$rhs = modelConstraints$rhs
  gurobiModel$vtype = variables
  gurobiModel$modelsense = "min"

  # solved = gurobi(gurobiModel)

  varNames = variable_keys(model)
  return(gurobiModel)
  # names(solved$x) = varNames
  # result = new_solution(model, solved$objval, tolower(solved$status), solved$x)
}
# ##



# gm$rhs

# str(objective_function(full_model))


# roi_model = as_ROI_model(full_model)

# str(roi_model)


# test_demand_data

# str(roi_model)

# roi_model$constraints$L

library(ROI)
library(Matrix)

create_availability_constraint = function(j, n, m) {
  y_matrix = sparseMatrix(
    i = 1:n, 
    j = rep(j, n), 
    x = -1, 
    dims = c(n, m))


  x_index = m*(1:n) - (m - j)


  x_matrix = sparseMatrix(
    i = 1:n,
    j = x_index,
    x = rep(1, n), 
    dims = c(n, m*n)
  )

  return(lst(y_matrix, x_matrix))
}



create_sum_constraint = function(i, n, m) {
  x_start_index = m*i  - (m - 1)
  x_end_index = m*i

  x_matrix = sparseMatrix(
    i = rep(1, m),
    j = x_start_index:x_end_index,
    x = rep(1, m),
    dims = c(1, m*n)
  )

  y_matrix = as(
    matrix(rep(0, m), nrow = 1),
    "sparseMatrix"
  )
  return(lst(
    y_matrix,
    x_matrix
  ))
}



create_takeup_constraint = function(takeup, n, m) {
  x_matrix = matrix(
    takeup,
    nrow = 1
  )
  y_matrix = as(
    matrix(rep(0, m), nrow = 1),
    "sparseMatrix"
  )
  return(lst(
    y_matrix,
    x_matrix
  ))
}


create_constraints = function(takeup, takeup_target, n, m) {
  availability_constraints = map(
    1:m,
    ~create_availability_constraint(
      j = .x,
      n = n,
      m = m
    )
  )
  sum_constraints = map(
    1:n,
    ~create_sum_constraint(
      i = .x,
      n = n,
      m = m
    )
  )

  takeup_constraint = create_takeup_constraint(
    takeup,
    n = n,
    m = m
  )


  x_constraint_matrix = rbind(
    do.call(rbind, map(sum_constraints, "x_matrix")),
    do.call(rbind, map(availability_constraints, "x_matrix")),
    takeup_constraint$x_matrix
  )



  y_constraint_matrix = rbind(
    do.call(rbind, map(sum_constraints, "y_matrix")),
    do.call(rbind, map(availability_constraints, "y_matrix")),
    takeup_constraint$y_matrix
  )

  constraint_matrix = cbind(
    y_constraint_matrix,
    x_constraint_matrix
  ) %>%
    as.matrix()


  dir = c(
    rep("==", n), 
    rep("<=", m*n),
    ">="
  ) 

  rhs = c(
    rep(1, n),
    rep(0, m*n),
    takeup_target*n
  )

  i_index = rep(1:n, each = m)
  j_index = rep(1:m, n)

  variable_names = c(
    paste0("y_", 1:m),
    paste0("x_", i_index, "_", j_index)
  )

  return(lst(
    constraint_matrix,
    dir,
    rhs,
    variable_names
    ))
}

create_objective = function(n, m){
  objective_matrix = c(rep(1, m), rep(0, n*m))
  i_index = rep(1:n, each = m)
  j_index = rep(1:m, n)

  variable_names = c(
    paste0("y_", 1:m),
    paste0("x_", i_index, "_", j_index)
  )
  return(lst(
    objective_matrix,
    variable_names
  ))
}



sorted_takeup = test_demand_data %>%
  arrange(pot_j, village_i) %>%
  pull(demand)

problem_constraints = create_constraints(
  takeup = sorted_takeup,
  takeup_target = script_options$target_constraint,
  n = data$n,
  m = data$m
)




problem_objective = create_objective(data$n, data$m)

gm = get_model_stuff(full_model) 






ed_lp = OP(
  objective = L_objective(
    problem_objective$objective_matrix,
    names = problem_objective$variable_names
  ), 
  constraints = L_constraint(
    L = problem_constraints$constraint_matrix,
    rhs = problem_constraints$rhs,
    dir = problem_constraints$dir,
    names = problem_constraints$variable_names
  ),
  types = rep(
    "B",
    length(problem_objective$objective_matrix)
  ),
  maximum = FALSE
)


ed_solution <- ROI_solve(
  ed_lp, 
  solver = "glpk",
  verbose = TRUE
)




full_soln = solve_model(
  full_model,
  with_ROI(
    solver = "glpk", 
    verbose = TRUE, 
    # control = list(tm_limit = script_options$time_limit)
    )
    )


tidy_full_soln = full_soln %>%
  get_solution(x[i,j]) %>%
  filter(value > .9) %>%  
  select(i, j) %>%
  as_tibble() %>%
  clean_output(
    match_df = .,
    data = data,
    demand_data = test_demand_data
  ) %>%
  mutate(type = "auto soln")

clean_ed_soln = ed_solution$solution %>%
  enframe() %>%
  mutate(
    variable = if_else(str_detect(name, "x"), "x", "y") 
  ) %>%
  mutate( 
    index_i = if_else(
      variable == "y", 
      NA_integer_, 
      str_extract(name, "(?<=_)\\d+(?=_)") %>% as.integer
    ),
    index_j = if_else(
      variable == "y", 
      str_extract(name, "\\d+") %>% as.integer,
      str_extract(name, "\\d+$") %>% as.integer
    )
  )


tidy_ed_soln = clean_ed_soln %>%
  filter(variable == "x") %>%
  filter(value == 1) %>%
  select( 
    i = index_i, 
    j = index_j
  ) %>%
  clean_output(
    match_df = ., 
    data = data, 
    demand_data = test_demand_data
  ) %>%
  mutate(
    type = "ed soln"
  )


tidy_comp_soln = bind_rows(
  tidy_ed_soln,
  tidy_full_soln
)



library(scales)
#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(2)


data$village_locations %>%
    ggplot() +
    geom_sf() +
    geom_sf(
        data = data$pot_locations,
        shape = 17,
        size = 4, 
        alpha = 0.7
        # colour = hex[1]
        ) +
    theme_bw() +
    geom_segment(
        data = tidy_comp_soln, 
        aes(x = village_lon, 
            y = village_lat, 
            xend = pot_lon, 
            yend = pot_lat,
            colour = type
            ))  +
    theme_bw()  +
    facet_wrap(~type)
ggsave("temp-plots/ed-manual.png", width = 8, height = 6, dpi = 500)






tictoc::tic()
tidy_output = demand_data %>%
  mutate(
    model_output = future_map(
      demand_data,
      ~safe_define_and_solve_model(
        baseline_model,
        data,
        .x,
        optim_type = optim_type,
        target_constraint = script_options$target_constraint
      ),
      .progress = TRUE,
      .options = furrr_options(seed = TRUE)
     )
  )
tictoc::toc()

if (script_options$dry_run) {
  output_path = file.path(
    script_options$output_path, 
    str_glue("{script_options$output_filename}-{stat_type}-subsidy-{script_options$dry_run_subsidy}-optimal-allocation.rds"))
} else {
  output_path = file.path(
    script_options$output_path, 
    str_glue("{script_options$output_filename}-{stat_type}-optimal-allocation.rds"))
}



saveRDS(
  tidy_output,
  output_path
)

