library(tidyverse)
library(ggplot2)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)

set.seed(1234)

generate_data = function(...){
    grid_size <- 10000
    n <- 200
    customer_locations <- data.frame(
    id = 1:n,
    x = round(runif(n) * grid_size),
    y = round(runif(n) * grid_size)
    )

    m <- 50
    warehouse_locations <- data.frame(
    id = 1:m,
    x = round(runif(m) * grid_size),
    y = round(runif(m) * grid_size)
    )
    fixedcost <- abs(round(rnorm(m, mean = grid_size/20, sd = grid_size)))
    return(
        lst(
            grid_size,
            n,
            customer_locations,
            m,
            warehouse_locations,
            fixedcost
        )
    )
}


sim_data = generate_data()




transportcost <- function(i, j, customer_locations, warehouse_locations, subsidy = 0) {
  customer <- customer_locations[i, ]
  warehouse <- warehouse_locations[j, ]
  cost = round(sqrt((customer$x - warehouse$x)^2 + (customer$y - warehouse$y)^2))
  cost = pmax(cost*(1-subsidy), 0)
  return(cost)
}



p_customer = sim_data$customer_locations %>%
    ggplot(aes(
        x = x,
        y = y
    )) +
    geom_point() +
    geom_point(
        data = sim_data$warehouse_locations, 
        color = "red", 
        shape = 17,
        size = 4) +
    theme_bw() +
    labs(
        title = "Facility Location Problem",
        subtitle = "Black dots indicate villages. Red triangles indicate potential clinic locations."
    )

p_customer

ggsave(
    "optim/example-location-problem.png",
    width = 10,
    height = 10
)


define_MIPModel = function(sim_data, subsidy = 0) {
    n = sim_data$n 
    m = sim_data$m 
    customer_locations = sim_data$customer_locations
    warehouse_locations = sim_data$warehouse_locations
    fixedcost = sim_data$fixedcost

    model = MIPModel() %>%
        # 1 iff i gets assigned to warehouse j
        add_variable(x[i, j], i = 1:n, j = 1:m, type = "binary") %>%
        # 1 iff warehouse j is built
        add_variable(y[j], j = 1:m, type = "binary") %>%
        # maximize the preferences
        set_objective(
                sum_over(transportcost(i, 
                                    j, 
                                    customer_locations = customer_locations, 
                                    warehouse_locations = warehouse_locations,
                                    subsidy = subsidy) * x[i, j], i = 1:n, j = 1:m) + 
                        sum_over(fixedcost[j] * y[j], j = 1:m), "min") %>%
        # every customer needs to be assigned to a warehouse
        add_constraint(sum_over(x[i, j], j = 1:m) == 1, i = 1:n) %>% 
        # if a customer is assigned to a warehouse, then this warehouse must be built
        add_constraint(x[i,j] <= y[j], i = 1:n, j = 1:m)
    return(model)
}

model = define_MIPModel(sim_data, subsidy = 0)
sm_model = define_MIPModel(sim_data, subsidy = 0.5)

result = solve_model(model, with_ROI(solver = "glpk", verbose = TRUE))
sm_result = solve_model(sm_model, with_ROI(solver = "glpk", verbose = TRUE))


matching = result %>% 
  get_solution(x[i,j]) %>%
  filter(value > .9) %>%  
  select(i, j) %>%
  as_tibble()

sm_matching = sm_result %>%
    get_solution(x[i,j]) %>%
    filter(value > 0.9) %>%
    select(i, j) %>%
    as_tibble()

sm_tidy_output = sm_matching %>% 
  inner_join(sim_data$customer_locations %>% 
                rename(customer_x = x, customer_y = y), by = c("i" = "id")) %>% 
  inner_join(sim_data$warehouse_locations %>%
                rename(warehouse_x = x, warehouse_y = y), by = c("j" = "id")) %>%
    mutate( 
        fixed_cost = sim_data$fixedcost[j]
    )


tidy_output = matching %>% 
  inner_join(sim_data$customer_locations %>% 
                rename(customer_x = x, customer_y = y), by = c("i" = "id")) %>% 
  inner_join(sim_data$warehouse_locations %>%
                rename(warehouse_x = x, warehouse_y = y), by = c("j" = "id")) %>%
    mutate( 
        fixed_cost = sim_data$fixedcost[j]
    )



customer_count = matching %>% 
    group_by(j) %>% 
    summarise(n = n()) %>% 
    rename(id = j)
sm_customer_count = sm_matching %>% 
    group_by(j) %>% 
    summarise(n = n()) %>% 
    rename(id = j)


comp_output = bind_rows(
    tidy_output %>% mutate(type = "Original"),
    sm_tidy_output %>% mutate(type = "Subsidy")
)

comp_output

sim_data$customer_locations %>%
    ggplot(aes(
        x = x,
        y = y
    )) +
    geom_point() +
    geom_point(
        data = sim_data$warehouse_locations, 
        color = "red", 
        shape = 17,
        size = 4, 
        alpha = 0.2) +
    theme_bw() +
    labs(
        title = "Facility Location Problem",
        subtitle = "Black dots indicate villages. Red triangles indicate potential clinic locations."
    ) +
  geom_segment(
    data = comp_output, 
    aes(x = customer_x, 
        y = customer_y, 
        xend = warehouse_x, 
        yend = warehouse_y),
    alpha = 0.5) +
  facet_wrap(~type, ncol = 1)

removed_clinics = setdiff(tidy_output$j, sm_tidy_output$j)

sim_data$customer_locations %>%
    ggplot(aes(
        x = x,
        y = y
    )) +
    geom_point(alpha = 0.4) +
    geom_point(
        data = sim_data$warehouse_locations %>%
                mutate(clinic_removed = id %in% removed_clinics),
        aes(color = clinic_removed),
        shape = 17,
        size = 4, 
        alpha = 1) +
    theme_bw() +
    labs(
        title = "Facility Location Problem",
        subtitle = "Black dots indicate villages. Triangles indicate potential clinic locations."
    ) +
  geom_segment(
    data = comp_output, 
    aes(x = customer_x, 
        y = customer_y, 
        xend = warehouse_x, 
        yend = warehouse_y),
    alpha = 0.2) +
  facet_wrap(~type, ncol = 1) +
  labs(
    color = "Clinic Unnecessary After Social Planner Accounts For Subsidy") +
  theme(
    legend.position = "bottom"
  )


ggsave(
    "optim/example-location-solved-problem.png",
    width = 10,
    height = 10
)


comp_output %>%
    mutate(
        cost = sqrt((customer_x - warehouse_x)^2 + (customer_y - warehouse_y)^2),
        cost = pmax(cost*(1-0.5), 0)
    ) %>%
    group_by(
        type
    ) %>%
    summarise(
        total_cost = sum(cost) + sum(unique(fixed_cost)), 
        n_clinics = n_distinct(j)
    )

transportcost()


customer_count

sm_customer_count
