#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        compare-optim.R  [options] 

        Options:
        --input-path=<input-path>  The input path [default: {file.path('optim', 'data', 'agg-full-many-pots')}]
        --output-path=<output-path>  The output path [default: {file.path('temp-data')}]
        --posterior-median  Whether to compare posterior medians or entire posterior draws
        --many-pots
        --welfare-function=<welfare-function>  The welfare function [default: log]
        --model=<model>
"),
  args = if (interactive()) "
        --input-path=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots
        --output-path=optim/data/STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP/agg-full-many-pots
        --many-pots 
        --model=STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP
        --welfare-function=identity
                             
                             " else commandArgs(trailingOnly = TRUE)
) 

library(tidyverse)
library(sf)

input_path = script_options$input_path 
output_path = script_options$output_path 
dir.create(output_path)
many_pots = script_options$many_pots


if (script_options$posterior_median) {
    oa_files = fs::dir_ls(
        input_path,
        regexp = "median-optimal-allocation\\.rds"
    )
} else {
    oa_files = fs::dir_ls(
        input_path,
        regexp = "optimal-allocation-subset-long-data\\.rds"
    )
}

models_we_want = script_options$model 

util_f = eval(parse(text = script_options$welfare_function))

experimental_file = str_glue("target-rep-agg-{script_options$welfare_function}-cutoff-b-control-mu-control-STRUCTURAL_LINEAR_U_SHOCKS_PHAT_MU_REP-median-experimental-control-allocation-data.rds")
oa_df = map_dfr(
    oa_files, 
    read_rds,
    .id = "file"
) %>% as_tibble()  %>%
    mutate(
        cutoff_type = if_else(str_detect(file, "no-cutoff"), "no_cutoff", "cutoff"), 
        rep_type = if_else(str_detect(file, "suppress-rep"), "suppress_rep", "rep"), 
        allocation_type = "optimal",
        distance_constraint = str_extract(file, "(?<=distconstraint-)\\d+") %>% as.numeric,
        fix_type = case_when(
            str_detect(file, "mu\\d+-") ~ "fix-mu",
            str_detect(file, "delta\\d+-") ~ "fix-delta",
            TRUE ~ "no-fix"
        ),
        fix_distance = if_else(fix_type != "no-fix", str_extract(file, "(?<=[mu|delta])\\d+") %>% as.numeric, NA_real_)
    ) 

experimental_oa_df = read_rds(
    file.path(
        input_path, 
        experimental_file
    )
) %>%
 as_tibble()  %>%
 mutate(file = experimental_file) %>%
    mutate(
        cutoff_type = if_else(str_detect(file, "no-cutoff"), "no_cutoff", "cutoff"), 
        rep_type = if_else(str_detect(file, "suppress-rep"), "suppress_rep", "rep"), 
        allocation_type = "experimental"
    ) %>%
    rename(j = pot_id, i = village_id)

oa_df = bind_rows(
    oa_df, 
    experimental_oa_df
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
        static_vstar = str_detect(file, "static")
    )


subset_oa_df = oa_df %>%
    filter(model == script_options$model)  %>%
    filter(
        cutoff_type  == "cutoff"
    ) %>%
    filter(!is.na(distance_constraint) | allocation_type == "experimental")


if (!script_options$posterior_median) { # if all draws

    demand_df = subset_oa_df %>%
        group_by(
            private_benefit_z,
            visibility_z, 
            static_vstar,
            allocation_type,
            distance_constraint,
            fix_type,
            fix_distance,
            model, 
            rep_type,
            cutoff_type,
            i,
            j
        ) %>%
        summarise(
            util = sum(util_f(demand)),
            mean_demand = mean(demand), 
            min_demand = min(demand), 
            median_demand = median(demand),
            mean_dist = mean(dist),
            n_pot = n_distinct(j, na.rm = TRUE),
            target_optim = mean(target_optim)
        )  %>%
        mutate(
            overshoot = 100*(util/target_optim - 1)
        )
treat_v = "bracelet"

    fix_oa_df = subset_oa_df %>%
        filter(private_benefit_z == treat_v  & visibility_z == treat_v) %>%
        filter(allocation_type != "experimental") %>%
        filter(static_vstar == FALSE) %>%
        filter(rep_type == "rep") %>%
        group_by(
            private_benefit_z,
            visibility_z, 
            static_vstar,
            allocation_type,
            distance_constraint,
            fix_type,
            fix_distance,
            model, 
            rep_type,
            cutoff_type,
            i,
            j
        )  %>%
        mutate(any_fix = fix_type != "no-fix")

    fix_demand_df = 
    demand_df %>%
        filter(private_benefit_z == treat_v  & visibility_z == treat_v) %>%
        filter(allocation_type != "experimental") %>%
        filter(static_vstar == FALSE) %>%
        filter(rep_type == "rep")

    draw_n = 200

    fix_oa_df %>%
        filter(!is.na(fix_distance)) %>%
        # filter(fix_distance == 0) %>%
        ggplot(aes(
            x = dist,
            y = demand,
            colour = fix_type,
            group =  interaction(fix_type, draw)
        )) +
        geom_line(alpha = 0.5)  +
        geom_line(
            inherit.aes = FALSE,
            aes(x = dist, y = demand, group = draw),
            data = fix_oa_df %>%
                filter(fix_type == "no-fix") %>%
                ungroup() %>%
                select(-fix_distance, -fix_type),               
                alpha = 0.3,
                colour = "black"
        ) +
        facet_grid(fix_distance~fix_type) +
        theme_bw() +
        theme(legend.position = "bottom") 

    draw_n = 20

    fix_oa_df %>%
        filter(!is.na(fix_distance)) %>%
        filter(draw == draw_n) %>%
        ggplot(aes(
            x = dist,
            y = demand,
            colour = fix_type,
            group =  interaction(fix_type, draw)
        )) +
        geom_line(alpha = 1)  +
        geom_line(
            data = fix_oa_df %>%
                filter(fix_type == "no-fix") %>%
                filter(draw == draw_n) %>%
                ungroup() %>%
                select(-fix_distance),               
                alpha = 1
        ) +
        facet_wrap(~fix_distance, ncol = 1) +
        theme_bw() +
        theme(legend.position = "bottom") +
        geom_vline(
            data = tibble(
                xint = c(0, 1250, 2500),
                fix_distance = c(0, 1250, 2500)
            ),
            aes(xintercept = xint),
            linetype = "longdash"
        )
    ggsave(
        "temp-plots/fixing-plot.pdf",
        width = 10,
        height = 10
    )
stop()



    demand_df %>%
        filter(!is.na(fix_distance)) %>%
        ggplot(aes(
            x = mean_dist,
            y = median_demand,
            colour = fix_type,
            group =  interaction(fix_type)
        )) +
        geom_line(alpha = 0.5)  +
        geom_line(
            data = demand_df %>%
                filter(fix_type == "no-fix") %>%
                ungroup() %>%
                select(-fix_distance),               
                alpha = 1
        ) +
        facet_wrap(~fix_distance, ncol = 1) +
        theme_bw() +
        theme(legend.position = "bottom") 


    summ_optim_df = subset_oa_df %>%
        group_by(
            draw,
            private_benefit_z,
            visibility_z, 
            static_vstar,
            allocation_type,
            distance_constraint,
            fix_type,
            fix_distance,
            model, 
            rep_type,
            cutoff_type
        ) %>%
        summarise(
            util = sum(util_f(demand)),
            mean_demand = mean(demand), 
            min_demand = min(demand), 
            mean_dist = mean(dist),
            n_pot = n_distinct(j, na.rm = TRUE),
            target_optim = mean(target_optim)
        )  %>%
        mutate(
            overshoot = 100*(util/target_optim - 1)
        )



    post_summ_optim_df = summ_optim_df %>%
        group_by(
            private_benefit_z,
            visibility_z, 
            static_vstar,
            allocation_type,
            distance_constraint,
            fix_type,
            fix_distance,
            model, 
            rep_type,
            cutoff_type
        ) %>%
        select(-draw) %>%
        filter(
            model %in% models_we_want
        ) %>%
        filter(
            cutoff_type  == "cutoff"
        )   %>%
        summarise(
            across(
                everything(),
                list(
                    estimate = ~mean(.x, na.rm = TRUE),
                    CI = ~paste0("(", signif(quantile(.x, 0.025, na.rm = TRUE), 3), ", ", signif(quantile(.x, 0.975, na.rm = TRUE), 3), ")")
                )
            )
        )   %>%
        ungroup() %>%
        select(-model, -cutoff_type) %>%
        arrange(
            private_benefit_z,
            visibility_z,
            allocation_type,
            rep_type
        ) %>%
    rename(
        B_z = private_benefit_z, 
        mu_z = visibility_z
    )

clean_summ_optim_df = summ_optim_df %>%
    group_by(
        private_benefit_z,
        visibility_z, 
        static_vstar,
        distance_constraint,
        fix_type,
        fix_distance,
        model, 
        rep_type,
        cutoff_type
    ) %>%
    select(-draw) %>%
    filter(
        model %in% models_we_want
    ) %>%
    filter(
        cutoff_type  == "cutoff"
    ) 





clean_summ_optim_df %>%
    filter(rep_type == "rep") %>%
    filter(!is.na(distance_constraint), distance_constraint != 2500) %>%
    ggplot(aes(
        x = n_pot, 
        fill = interaction(private_benefit_z, visibility_z, static_vstar, rep_type, allocation_type, fix_type, fix_distance)
    )) +
    geom_density(alpha = 0.6) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    facet_wrap(~distance_constraint)



clean_summ_optim_df %>%
    filter(rep_type == "rep") %>%
    ggplot(aes(
        x = n_pot, 
        fill = interaction(private_benefit_z, visibility_z, static_vstar, rep_type, allocation_type)
    )) +
    geom_histogram(
        colour = "black"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")




clean_summ_optim_df %>% 
    write_csv(
        file.path(
            output_path,
            "posterior-clean-summ-optim.csv"
        ))

post_summ_optim_df %>% 
    filter(rep_type == "suppress_rep") %>%
    select(-rep_type) %>%
    write_csv(
        file.path(
            output_path,
            "posterior-suppress-rep-summ-optim.csv"
        ))
}

if (script_options$posterior_median) {


model_df = subset_oa_df %>%
    unnest(model_output) 



summ_optim_df = model_df %>%
    group_by(
        private_benefit_z,
        visibility_z, 
        model, 
        rep_type,
        cutoff_type
    ) %>%
    summarise(
        util = sum(log(demand)),
        mean_demand = mean(demand), 
        min_demand = min(demand), 
        n_pot = n_distinct(j),
        target_optim = mean(target_optim)
    ) %>%
    filter(
        model %in% models_we_want
    ) %>%
    filter(
        cutoff_type  == "cutoff"
    ) %>%
    mutate(
        overshoot = 100*(util/target_optim - 1)
    ) %>%
    select(-model, -cutoff_type) %>%
    ungroup() %>%
    select(-model) %>%
    arrange(
        private_benefit_z,
        visibility_z,
        rep_type
    )




table_summ_optim_df = summ_optim_df %>%
    rename(
        B_z = private_benefit_z, 
        mu_z = visibility_z
    )

table_summ_optim_df %>% 
    filter(rep_type == "rep") %>%
    select(-rep_type) %>%
    write_csv(
        file.path(
            output_path,
            "rep-summ-optim.csv"
        ))

table_summ_optim_df %>% 
    filter(rep_type == "suppress_rep") %>%
    select(-rep_type) %>%
    write_csv(
        file.path(
            output_path,
            "suppress-rep-summ-optim.csv"
        ))

subset_oa_df %>%
    filter(private_benefit_z == "control") %>%
    mutate(visibility_z = fct_relabel(visibility_z, str_to_title)) %>%
    filter(rep_type == "rep") %>%
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


demand_df = subset_oa_df %>%
    unnest(demand_data)

demand_df %>%
    filter(dist < 3500) %>%
    filter(private_benefit_z == "control") %>%
    filter(rep_type == "rep") %>%
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
    filter(dist < 3500) %>%
    mutate(
        private_benefit_z = paste0("B: ", private_benefit_z),
        visibility_z = paste0("Mu: ", visibility_z)
    ) %>%
    ggplot(aes(
        x =  dist,
        y = demand, 
        colour = rep_type
    )) +
    geom_line()  +
    geom_point()  +
    theme(legend.position = "bottom")  +
    facet_grid(
        visibility_z ~ private_benefit_z
    ) +
    theme_bw() +
    theme(legend.position =  "bottom")


ggsave(
    file.path(
        output_path,
        "struct-naive-optim-demand-df.png"
    ),
    width = 8,
    height = 6,
    dpi = 500
)
}




