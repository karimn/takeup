#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        create-optim-paper-panel.R  [options] 

        Options:
        --model=<model>
        --output-path=<output-path>
        --fit-version=<fit-version>
        --welfare-function=<welfare-function> Which utility function to use [default: log]
        --input-path=<input-path>
        --distance-constraint=<distance-constraint>  Maximum distance DM can send someone in meters. [default: 3500]
"),
  args = if (interactive()) "
                            --input-path=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots \
                            --output-path=optim/plots/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots \
                            --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP \
                            --fit-version=86 \
                            --welfare-function=identity \
                            --distance-constraint=3500
                             " else commandArgs(trailingOnly = TRUE)
) 
library(tidyverse)
library(data.table)
library(sf)


dist_data = read_rds(
  file.path(
    "optim/data",
    "full-many-pots-experiment.rds" 
  )
)

village_data = dist_data$village_df
pot_data = dist_data$pot_df

n = nrow(village_data)
m = nrow(pot_data)


wgs.84 = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
## Loading Data
data = list(
    n = n,
    m = m,
    pot_locations = st_as_sf(pot_data, coords = c("lon", "lat"), crs = wgs.84),
    village_locations = st_as_sf(village_data, coords = c("lon", "lat"), crs = wgs.84)
) 

optim_input_path = script_options$input_path 
optim_files = c(
    str_glue("target-rep-distconstraint-{script_options$distance_constraint}-util-{script_options$welfare_function}-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-optimal-allocation.rds"),
    str_glue("target-rep-distconstraint-{script_options$distance_constraint}-util-{script_options$welfare_function}-cutoff-b-control-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-optimal-allocation.rds"),
    str_glue("target-rep-distconstraint-{script_options$distance_constraint}-util-{script_options$welfare_function}-static-cutoff-b-control-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-optimal-allocation.rds"),
    str_glue("target-rep-distconstraint-{script_options$distance_constraint}-util-{script_options$welfare_function}-suppress-rep-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-optimal-allocation.rds")
)
# experimental_file = str_glue("target-rep-util-{script_options$welfare_function}-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-experimental-control-allocation-data.rds")
experimental_file = str_glue(
    "target-rep-agg-{script_options$welfare_function}-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-experimental-control-allocation-data.rds"
)
experimental_demand = read_rds(
    file.path(
        optim_input_path,
        experimental_file
    )
)

experimental_allocation = experimental_demand %>%
    group_by(
        j = pot_id, 
        i = village_id
    ) %>%
    summarise(
        pot_lon = first(pot_lon), 
        pot_lat = first(pot_lat), 
        village_lon = first(village_lon), 
        village_lat = first(village_lat), 
        demand = mean(demand), 
        dist = mean(dist)
    ) %>%
    ungroup() %>%
    nest(data = everything()) %>%
    rename(
        model_output = data
    ) %>%
    mutate(
        private_benefit_z = "control", 
        visibility_z = "control", 
        model = "STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP",
        status_code = 0, 
        allocation_type = "experimental"
    ) %>%
    select( 
        private_benefit_z, 
        visibility_z, 
        model, 
        model_output, 
        status_code, 
        allocation_type
    )



optimisation_df = 
    file.path(
        optim_input_path, 
        optim_files
    ) %>%
    map_dfr(read_rds) %>%
    select(-demand_data, -optim_problem, -optim_fit)   

optimisation_df[, allocation_type := "policymaker"] 

optimisation_df = rbindlist(
    list(
        optimisation_df,
        experimental_allocation
    ) 
)


optimisation_df[, file := c(optim_files, experimental_file) ]

optimisation_df[, n_pots_used := lapply(model_output, function(x){length(unique(x$j))})]

optimisation_df[, suppress_rep := str_detect(file, "suppress-rep")]
optimisation_df[, static_vstar := str_detect(file, "static")]


optimisation_df[, pot_data := list(data$pot_locations)]
optimisation_df[suppress_rep == TRUE, visibility_z := "No Visibility"]
optimisation_df[
    static_vstar == TRUE & visibility_z == "bracelet", 
    visibility_z := "Signal value fixed at bracelet 0.5km"]


long_opt_df = bind_rows(
    optimisation_df %>%
        as_tibble()   %>%
        filter(status_code != 1 &  allocation_type == "policymaker") %>%
        unnest(cols = c(model_output)),
    optimisation_df %>%
        as_tibble()   %>%
        filter(status_code == 1 & allocation_type == "policymaker") %>%
        unnest(cols = c(model_output)),
    optimisation_df %>%
        as_tibble()   %>%
        filter(allocation_type == "experimental") %>%
        unnest(cols = c(model_output))
)



library(scales)
#extract hex color codes for a plot with three elements in ggplot2 
hex <- ggthemes::canva_pal("Primary colors with a vibrant twist")(2)



plot_comp_function = function(meta_df, long_df){
    # long_df = subset_long_opt_df

    # meta_df = subset_optimisation_df

    summ_df = long_df %>%
        group_by(
            private_benefit_z, 
            visibility_z,
            model,
            allocation_type
        ) %>%
        summarise(
            mean_demand = round(100*mean(demand, na.rm = TRUE), 2),
            mean_dist = mean(dist),
            n_pot = n_distinct(j)
        ) 

    optimal_comp_df = long_df

    meta_df[, pot_data := map2(
        model_output, 
        pot_data, 
        ~mutate(
            .y, 
            assigned_pot = id %in% unique(.x$j)
        )
        )]

    meta_df[, mu_z := factor(visibility_z)]     
    meta_df[, B_z := factor(private_benefit_z, levels = c("control", "ink", "calendar", "bracelet"))]
    n_row = nrow(meta_df)
    pot_data_list = map(1:n_row, ~meta_df[.x, pot_data][[1]])


    pot_comp_df = bind_rows(
        pot_data_list %>%
            imap(
                ~.x %>% mutate(
                    B_z = meta_df[.y, private_benefit_z], 
                    mu_z = meta_df[.y, visibility_z],
                    allocation_type = meta_df[.y, allocation_type]
                )
                )
            ) %>%
            mutate(
                type = paste0("pot_", assigned_pot)
            )



    village_comp_df = bind_rows(
        map(
            n_row, 
            ~ {
                village_data %>% mutate(type = "village") %>% select(-cluster.id) %>%
                    mutate(
                        B_z = meta_df[.x, private_benefit_z], 
                        mu_z = meta_df[.x, visibility_z], 
                        allocation_type = meta_df[.x, allocation_type]
                    )
            }
        )
    )
    
    


    optimal_comp_df = optimal_comp_df %>%
        mutate(B_z = private_benefit_z, mu_z = visibility_z) %>%
        mutate(
            mu_z = factor(mu_z), # levels = c("No Visibility", "control", "ink", "calendar", "bracelet")),
            cf_type = str_glue("B({B_z}, d),   mu({mu_z}, d) - {allocation_type} solution")
        ) 



    vil_pot_comp_df = bind_rows(
        pot_comp_df,
        village_comp_df
    ) %>%
    mutate(
        mu_z = factor(mu_z), # levels = c("No Visibility", "control", "ink", "calendar", "bracelet")),
        cf_type = str_glue("B({B_z}, d),   mu({mu_z}, d) - {allocation_type} solution")
    ) 

    summ_df = summ_df %>%
        rename(B_z = private_benefit_z, mu_z = visibility_z) %>%
        mutate(
            cf_type = str_glue("B({B_z}, d),   mu({mu_z}, d) - {allocation_type} solution")
        ) 


    meta_df = meta_df %>%
        as_tibble() %>%
        mutate(
            cf_type = str_glue("B({B_z}, d),   mu({mu_z}, d) - {allocation_type} solution")
        ) 

    meta_df = meta_df %>%
        mutate(
            levels = 1:n()
        )

    meta_df = meta_df %>%
        mutate(
            cf_type = fct_reorder(cf_type, levels)
        )

    relvl_cf_type = function(data) {
        data %>%
            mutate(
                cf_type = factor(cf_type, levels(meta_df$cf_type))
            )
    }
    vil_pot_comp_df = relvl_cf_type(vil_pot_comp_df)
    optimal_comp_df = relvl_cf_type(optimal_comp_df)
    summ_df = relvl_cf_type(summ_df)

    ed_label = function(string){string}

    plot_fun = function(vil_pot_df, opt_df, summ_df) {
        vil_pot_df %>%
            ggplot() +
            geom_sf(aes(
                color = type, 
                shape = type, 
                size = type, 
                alpha = type
            )) +
            geom_segment(
                data = opt_df, 
                aes(x = village_lon, 
                    y = village_lat, 
                    xend = pot_lon, 
                    yend = pot_lat, 
                    linetype = "Assigned Village-PoT Pair"
                    ))  +
            theme_bw()  +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            # facet_wrap(~cf_type, labeller = as_labeller(ed_label)) +
            geom_text(
                    data = summ_df,
                    x = 34.8, 
                    y = 0.75,
                aes(label = paste0("Assigned PoTs: ", n_pot))
            ) +
            # geom_text(
            #         data = summ_df,
            #         x = 34.77, 
            #         y = 0.0,
            #     aes(label = str_glue("Mean Distance: {round(mean_dist/1000, 2)}km"))
            # ) +
            labs(
                title = unique(opt_df$cf_type),
                subtitle = str_glue(
                    # "Assigned PoTs: {unique(summ_df$n_pot)}; Mean Distance: {round(unique(summ_df$mean_dist)/1000, 2)}km; Takeup: {round(unique(summ_df$mean_demand), 2)}"
                    "Mean Distance: {round(unique(summ_df$mean_dist)/1000, 2)}km; Takeup: {round(unique(summ_df$mean_demand), 2)}"
                ),
                x = "", y = "") +
            scale_size_manual(
                values = c("pot_TRUE" = 1.5, "pot_FALSE" = 1 , "village" = 1),
                breaks = c("pot_TRUE", "pot_FALSE", "village"),
                labels = c("PoT Used", "PoT Unused", "Village")
            ) +
            scale_alpha_manual(
                values = c("pot_TRUE" = 1, "pot_FALSE" = 0.3, "village" = 0.8),
                breaks = c("pot_TRUE", "pot_FALSE", "village"),
                labels = c("PoT Used", "PoT Unused", "Village")
            ) +
            scale_color_manual(
                values = c("pot_TRUE" = hex[2], "pot_FALSE" = hex[1], "village" = "black"),
                breaks = c("pot_TRUE", "pot_FALSE", "village"),
                labels = c("PoT Used", "PoT Unused", "Village")
            ) +
            scale_shape_manual(
                values = c("pot_TRUE" = 17, "pot_FALSE" = 17, "village" = 16),
                breaks = c("pot_TRUE", "pot_FALSE", "village"),
                labels = c("PoT Used", "PoT Unused", "Village")
            ) +
            scale_linetype_manual("Assigned Village-PoT Pair",values=c("Assigned Village-PoT Pair"=1)) +
            theme(legend.title = element_blank()) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            # theme(legend.position = "bottom") +
            theme(legend.position = "none") +
            guides(
                alpha = guide_legend(override.aes = list(alpha = 1))
            ) 
    }


    cf_levels = unique(vil_pot_comp_df$cf_type)     

    plots = map(
        cf_levels, 
        ~plot_fun(
            vil_pot_df = vil_pot_comp_df %>% filter(cf_type == .x | type == "village"),
            optimal_comp_df %>%
                filter(cf_type == .x), 
            summ_df %>%
                filter(cf_type == .x)
        )
    )


    library(cowplot)

    legend <- get_legend(
    # create some space to the left of the legend
    plots[[3]] + theme(
        legend.position = "bottom", 
        legend.box.margin = margin(0, 0, 0, 0, "cm")
        )
    )

    plots[[2]] = plots[[2]] + 
        theme(legend.position = "bottom") 
# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
    p_grid = plot_grid(
        plotlist = plots, 
        align = 'hv',
        axis = "tb",
        # hjust = -1,
        nrow = 1
    )
    p_grid

    return(p_grid)
}


subset_optimisation_df = optimisation_df %>%
    filter(
        file %in% c(
    str_glue("target-rep-distconstraint-{script_options$distance_constraint}-util-{script_options$welfare_function}-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-optimal-allocation.rds"),
    str_glue("target-rep-distconstraint-{script_options$distance_constraint}-util-{script_options$welfare_function}-cutoff-b-control-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-optimal-allocation.rds")
        ) | allocation_type == "experimental"
    ) %>%
    arrange(allocation_type)

subset_long_opt_df = long_opt_df %>%
    filter(
        file %in% c(
    str_glue("target-rep-distconstraint-{script_options$distance_constraint}-util-{script_options$welfare_function}-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-optimal-allocation.rds"),
    str_glue("target-rep-distconstraint-{script_options$distance_constraint}-util-{script_options$welfare_function}-cutoff-b-control-mu-bracelet-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-optimal-allocation.rds")
        ) | allocation_type == "experimental"
    ) %>%
    arrange(allocation_type)

comp_plot = plot_comp_function(subset_optimisation_df, subset_long_opt_df)


save_plot(
    plot = comp_plot,
    filename = file.path(
        str_glue("presentations/optim-takeup-{script_options$model}-fig"), 
        str_glue(
            "panel-scenarios-compare-optimal-allocation-plot-distconstraint-{script_options$distance_constraint}.pdf"
        )
    ), 
    base_width = 20, 
    base_height = 20
)
