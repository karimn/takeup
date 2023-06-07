library(tidyverse)

input_dir = "optim/plots/agg-log-full-many-pots"
stat_files = fs::dir_ls(input_dir, regexp = "*stat*")


stat_files = stat_files[str_detect(stat_files, "median")]


stat_df = map_dfr(stat_files, read_rds, .id = "file") 


stat_df %>%
    pull(file) %>%
    head(1)


colnames(stat_df)


stat_df %>%
    select(util_target) %>%
    pull()


stat_df %>%
    mutate(
        rep_type = if_else(str_detect(file, "suppress"), "suppress_rep", "rep"), 
        b_type = str_extract(file, "(?<=b-).*(?=-mu)"), 
        mu_type = str_extract(file, "(?<=mu-).*(?=-STRUCT)")
    ) %>%
    mutate(
        fig_cap = str_glue(
            "Optimal Allocation of Points of Treatment: $B_{{z = \\textrm{{{b_type}}}}}, \\mu_{{z = \\textrm{{{mu_type}}}}}$"
        ), 
        fig_note = str_glue(
            "Social welfare function uses log utility.
            Utility target: {util_target}, utility achieved: {util_hit}, percentage over 
            constraint: {overshoot}\\%. Take-up achieved: {takeup_hit}. Mean walking distance: {mean_dist}"
        )
    )  %>% 
    select(
        rep_type, 
        b_type, 
        mu_type, 
        fig_cap, 
        fig_note
    ) %>%
    kableExtra::kbl(
        format = "latex", 
        escape = FALSE, 
        col.names = c(
            "rep type", 
            "$B$ type", 
            "\\mu type", 
            "fig.cap", 
            "fig.note"
        )
    ) %>%
    kableExtra::save_kable("optim-notes.tex")


