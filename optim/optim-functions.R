
################################################################################
# Functions ####################################################################
################################################################################

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
    if (optim_type == "min_cost_indiv") {
      model = model %>%
          # Takeup must be at least 
          add_constraint(
            sum_over(x[i, j]*demand_data[village_i == i & pot_j == j, demand], j = 1:m) >=  target_constraint, i = 1:n)  %>%
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
  # dd = copy(demand_data)
  # assigned_pots = unique(tidy_output$j)
  # new_tidy_output = dd[pot_j %in% assigned_pots, min_dist_avail := min(dist), village_i][
  #   dist == min_dist_avail
  #   ] %>%
  #   inner_join(data$village_locations %>% 
  #                   rename(village_lon = lon, village_lat = lat), by = c("village_i" = "id")) %>% 
  #   inner_join(data$pot_locations %>%
  #                   rename(pot_lon = lon, pot_lat = lat), by = c("pot_j" = "id"))  %>%
  #   rename(
  #     i = village_i, 
  #     j = pot_j
  #   ) %>%
  #   as_tibble() %>%
  #   select(-min_dist_avail)
  # if (nrow(tidy_output) == nrow(new_tidy_output)) {
  #   tidy_output = new_tidy_output
  # } else {
  #   warning("Temporary fix failed, using original allocation.")
  # }
    return(tidy_output)
}


#### Manual LP Creation Functions


find_x_index = function(n, m, i, j) {
  x_idx = m*i - (m - j)
  return(x_idx)
}

#' Create Availability Constraints
#' 
#' i.e. Village i can't use PoT j unless PoT j switched on.
#' 
#' x_{i,j} <= y_j for all i, j 
#'
#' y_j matrix is for n = 3, m = 2
#' 
#' -1 0
#' -1 0
#' -1 0
#'  0 -1
#'  0 -1
#'  0 -1
#'
#' x_matrix a bit tricky. Think of block matrix:
#' [ X_{vill 1} X_{vill 2} X_{vill 3} ]
#' For row j of constraint matrix we want the 
#' jth column from each village block matrix equal to 1 and 0 
#' otherwise
create_availability_constraint = function(n, m) {
  y_matrix = simple_triplet_matrix(
    i = 1:(n*m), 
    j = rep(1:m, each = n),
    v = rep(-1, n*m), 
    nrow = n*m,
    ncol = m)

  x_index = map(
    1:m,
    ~find_x_index(n = n, m = m, i = 1:n, j = .x)
  ) %>%
    unlist()

  i_mat = 1:(n*m)

  x_matrix = simple_triplet_matrix(
    i = i_mat, 
    j = x_index,
    v = rep(1, length(x_index)),
    nrow = n*m,
    ncol = n*m
  )

  return(lst(y_matrix, x_matrix))
}


#' Create Summation Constraint
#'
#' \sum^m x_{ij} = 1 \forall i
#'
#' i.e. every village is assigned to a single PoT
#' Format, n = 2 (village) m = 3 (PoT)
#'  x_11 x_12 x_13 x_21 x_22 x_23
#'   1    1     1   0    0    0
#'   0    0     0   1    1    1 
#' 
#' First index is village, second index is PoT
create_sum_constraint = function(n, m) {

  i_new = 1:n
  xsi = m*i_new - (m - 1) # x_start_index
  xei = m*i_new # x_end_index

  xi = map2(xsi, xei, ~.x:.y) %>%
    unlist()
  # index for constraint matrix ROW
  i_mat = rep(1:n, each = m)

  x_matrix = simple_triplet_matrix(
    i = i_mat,
    j = xi,
    v = rep(1, length(xi)),
    nrow = n,
    ncol = m*n
  )

  y_matrix = simple_triplet_zero_matrix(
    nrow = n,
    ncol = m
  )
  return(lst(
    y_matrix,
    x_matrix
  ))
}


#' Create takeup constraint
#'
#'
#' demand_data %>%
#'  unnest(demand_data) %>%
#'  arrange(pot_j, village_i) %>%
#'  select(pot_j, village_i)
#'
#' Takeup vector in order:
#' pot_1 vill_1
#' pot_1 vill_2
#' pot_1 vill_3
#' ...
#' pot_2 vill_1
#' pot_2 vill_2
#' pot_2 vill_3
#' ...
#' pot_m vill_i
#'
#' ## min-cost-agg ##
#' On average takeup must be greater than or equal to a target
#' 
#' Therefore A_takeup:
#' y_1 ... y_m x_11 ... x_ij ... x_nm
#' 0 ...   0   t_11 ... t_ij ... t_nm
#' i.e. summation of takeups (on RHS we will have target*n) to give average.
#' A_takeup  X = sum(takeup x_{ij, switched on}) 
#' 
#' ## min-cost-indiv ##
#' For each town, takeup must be greater than or equal to a target
#' 
#' Therefore A_takeup:
#' y_1 ... y_m x_11 x_12... x_21 x_22 ... x_nm
#' 0 ...   0   t_11 t_12 ...0    0    ... 0
#' 0 ...   0   0    0 ...   t_21    t_22 ... 0
#' (Basically a block diagonal matrix with each block a village's demand)
#' i.e. 0 everywhere but if row(A) == i, x_ji = t_ji
#' i.e. summation of takeups by row (on RHS we will have target_i for each 
#' town i).
#' A_takeup  X = vector(takeup x_{ij, switched on}) >= vector(target_i)
#' 
create_takeup_constraint = function(takeup, n, m, constraint_type, takeup_target) {
  rhs = takeup_target

  if (constraint_type == "agg") {
    x_matrix = matrix(
      takeup,
      nrow = 1
    )
    dir = ">="
    y_nrow = 1
  }
  if (constraint_type == "indiv") {

    xi = map(1:n, ~find_x_index(n = n, m = m, i = .x, j = 1:m)) %>%
      unlist()
    i_mat = rep(1:n, each = m)
    

    x_matrix = simple_triplet_matrix(
      i = i_mat, 
      j = xi, 
      v = takeup, 
      nrow = n, 
      ncol = m*n
    )

    # js = rep(1:n, each = m)
    # x_index = n*(1:m) - (n - js)

    # takeup_indices = map(1:n, ~x_index[js == .x]) %>%
    #   unlist()
    # x_matrix = simple_triplet_matrix(
    #   i = js,
    #   j = x_index,
    #   v = takeup[takeup_indices],
    #   nrow = n,
    #   ncol = m*n
    # )

    # x_matrix %>% as.matrix()

    dir = rep(">=", n)
  
    # Suppose we want to find demand in village 2 for all PoTs:
    # vill_we_want = 2
    # x_matrix[vill_we_want, find_x_index(n = n, m = m, i = vill_we_want, j = 1:m) ]$v

    # demand_data %>%
    #   unnest(demand_data) %>%
    #   group_by(village_i) %>%
    #   mutate(new_i_id = cur_group_id()) %>%
    #   arrange(pot_j, new_i_id) %>%
    #   select(pot_j, new_i_id, village_i, demand) %>%
    #   filter(new_i_id == vill_we_want) %>%
    #   filter(demand > 0)
    y_nrow = n
  }

  y_matrix = simple_triplet_zero_matrix(
    nrow = y_nrow,
    ncol = m
  )
  return(lst(
    y_matrix,
    x_matrix,
    rhs,
    dir
  ))
}


create_base_constraints = function(n, m) {

  availability_constraints = create_availability_constraint(
    n = n, 
    m = m
  )
  
  sum_constraints = create_sum_constraint(
    n = n, 
    m = m
  )
  



  x_constraint_matrix = rbind(
    sum_constraints$x_matrix,
    availability_constraints$x_matrix
  )



  y_constraint_matrix = rbind(
    sum_constraints$y_matrix,
    availability_constraints$y_matrix
  )

  constraint_matrix = cbind(
    y_constraint_matrix,
    x_constraint_matrix
  ) 


  dir = c(
    rep("==", n), 
    rep("<=", m*n)
  ) 

  rhs = c(
    rep(1, n),
    rep(0, m*n)
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

add_takeup_constraints = function(takeup, 
                                  takeup_target, 
                                  baseline_constraints, 
                                  n, 
                                  m, 
                                  constraint_type) {
  takeup_constraint = create_takeup_constraint(
    takeup,
    n = n,
    m = m,
    constraint_type =  constraint_type,
    takeup_target = takeup_target
  )

  baseline_constraints$constraint_matrix = rbind(
    baseline_constraints$constraint_matrix,
    cbind(takeup_constraint$y_matrix, takeup_constraint$x_matrix)
  ) 


  baseline_constraints$dir = c(
    baseline_constraints$dir,
    takeup_constraint$dir
  )

  baseline_constraints$rhs = c(
    baseline_constraints$rhs,
    takeup_constraint$rhs
  )

  return(baseline_constraints)
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



#' Define Model
#' 
#'
#' Solve in form Ax ><= B
#' A columns take form:
#'  y_1 y_2 ... y_m x_11 x_12 ... x_21 x_22 .. x_m1 ... x_mn
#' 
#' i.e. first m cols are PoT switched on indicators.
#' Next m*n variables are indicators for village i-PoT j pairing. 
#' First m sets are for village 1, next m sets for village 2 etc..
#' 
#' A rows take form: 
#'  summation constraints  (every village must have a PoT)
#'  availability constraints (every village can only use a switched on PoT)
#'  takeup constrains (we must hit a certain takeup target)
#'
define_model = function(takeup, data, target_takeup, baseline_constraints, constraint_type) {
  # sorted_takeup = takeup %>%
  #   arrange(pot_j, village_i) %>%
  #   pull(demand)
  # takeup
  # data
  sorted_takeup = takeup %>%
    arrange(village_i, pot_j) %>%
    pull(util)

  problem_constraints = add_takeup_constraints(
    takeup = sorted_takeup,
    takeup_target = target_takeup,
    n = data$n,
    m = data$m,
    baseline_constraints = baseline_constraints,
    constraint_type = constraint_type
  )

  problem_objective = create_objective(data$n, data$m)

  linear_programme = OP(
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
  return(linear_programme)
}








clean_solution = function(model_fit, data, takeup) {
  clean_soln = model_fit$solution %>%
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

  clean_model_output = clean_soln %>%
    filter(variable == "x") %>%
    filter(value == 1) %>%
    select( 
      i = index_i, 
      j = index_j
    ) %>%
    clean_output(
      match_df = ., 
      data = data, 
      demand_data = takeup
    ) 

  return(clean_model_output)
}
