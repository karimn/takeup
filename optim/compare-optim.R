#!/usr/bin/Rscript
script_options <- docopt::docopt(
    stringr::str_glue("Usage:
        compare-optim.R  [options] 

        Options:
        --input-path=<input-path>  The input path [default: {file.path('optim', 'data', 'agg-log-full')}]
        --output-path=<output-path>  The output path [default: {file.path('temp-data')}]
        --posterior-median  Whether to compare posterior medians or entire posterior draws
        --many-pots
"),
  args = if (interactive()) "
        --input-path=~/projects/takeup/optim/data/agg-log-full-many-pots
        --output-path=~/projects/takeup/optim/data/agg-log-full-many-pots
        --many-pots
                             
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

models_we_want = "STRUCTURAL_LINEAR_U_SHOCKS"

oa_df = map_dfr(
    oa_files, 
    read_rds,
    .id = "file"
) %>% as_tibble()  %>%
    mutate(
        cutoff_type = if_else(str_detect(file, "no-cutoff"), "no_cutoff", "cutoff"), 
        rep_type = if_else(str_detect(file, "suppress-rep"), "suppress_rep", "rep")
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
        cutoff_type  == "cutoff"
    ) 


if (!script_options$posterior_median) { # if all draws
    summ_optim_df = subset_oa_df %>%
        group_by(
            draw,
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
        )  %>%
        mutate(
            overshoot = 100*(util/target_optim - 1)
        )

    post_summ_optim_df = summ_optim_df %>%
        group_by(
            private_benefit_z,
            visibility_z, 
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
        ) %>%
        summarise(
            across(
                everything(),
                list(
                    estimate = mean,
                    CI = ~paste0("(", signif(quantile(.x, 0.025), 3), ", ", signif(quantile(.x, 0.975), 3), ")")
                )
            )
        )   %>%
        ungroup() %>%
        select(-model, -cutoff_type) %>%
        arrange(
            private_benefit_z,
            visibility_z,
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




