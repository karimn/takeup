script_options = docopt::docopt(
    stringr::str_glue("Usage:
    analyse-reduced-form-sms.R  [options]

    Options:
        --load-fit
        --save-fit
        --fit-path=<fit-path>
        --fit-file=<fit-file>
        --include-reminder
        --output-path=<output-path>  Output path [default: temp-data]
    "),
    args = if (interactive()) "
        --load-fit
        --fit-path=data/stan_analysis_data
        --fit-file=SMS_BRMS_reminder_fit.rds
        --output-path=temp-data/sms-noreminder
        --include-reminder
    " else commandArgs(trailingOnly = TRUE)
    # args = if (interactive()) "takeup cv --models=REDUCED_FORM_NO_RESTRICT --cmdstanr --include-paths=stan_models --update --output-path=data/stan_analysis_data --outputname=test --folds=2 --sequential" else commandArgs(trailingOnly = TRUE)
) 

library(tidyverse)
library(posterior)
library(tidybayes)
library(brms)
library(marginaleffects)
library(ggthemes)

options(mc.cores = 4)

canva_palette_vibrant <- "Primary colors with a vibrant twist"


source("analysis_util.R")
source(file.path("multilvlr", "multilvlr_util.R"))
source("dist_structural_util.R")

# Data --------------------------------------------------------------------

load(file.path("data", "analysis.RData"))

standardize <- as_mapper(~ (.) / sd(.))
unstandardize <- function(standardized, original) standardized * sd(original)
# stick to monitored sms.treatment group
# remove sms.treatment.2
monitored_nosms_data <- analysis.data %>% 
  filter(mon_status == "monitored", sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()

monitored_sms_data <- analysis.data %>% 
  filter(mon_status == "monitored") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()



nosms_data <- analysis.data %>% 
  filter(sms.treatment.2 == "sms.control") %>% 
  left_join(village.centers %>% select(cluster.id, cluster.dist.to.pot = dist.to.pot),
            by = "cluster.id") %>% 
  mutate(standard_cluster.dist.to.pot = standardize(cluster.dist.to.pot)) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster_id = cur_group_id()) %>% 
  ungroup()


analysis_data <- monitored_sms_data %>% 
    mutate(
    assigned_treatment = assigned.treatment, 
    assigned_dist_group = dist.pot.group, 
    sms_treatment = sms.treatment.2, 
    phone_owner = if_else(phone_owner == TRUE, "phone", "nophone"), 
    sms_treatment = str_replace_all(sms_treatment, "\\.", "")) %>%
    # reminder.only only present in control condition
    filter(phone_owner == "phone") %>%
    mutate(sms_treatment = factor(sms_treatment))

if (!script_options$include_reminder) {
    analysis_data = analysis_data %>%
        filter(sms_treatment != "reminderonly")
}



default_priors = get_prior(
    data = analysis_data,
    dewormed ~ (assigned_treatment*assigned_dist_group | sms_treatment), 
    family = bernoulli(link = "probit")
)

manual_priors = default_priors

#TODO: gen this list of variable names
# coefs = c(prior_summary(rf_fit) %>%
#     as_tibble() %>%
#     filter(coef != "" & coef != "Intercept") %>%
#     select(coef) %>%
#     pull(coef))

# manual_priors[which(default_priors$coef %in% coefs), "prior"] = "normal(0, 0.25)"
# manual_priors[which(default_priors$coef == "Intercept"), "prior"] = "normal(0, 1)"


if (script_options$load_fit) {
    rf_fit = read_rds(
        file.path(
            script_options$fit_path,
            script_options$fit_file
        )
    )
} else {
    rf_fit = brm(
        data = analysis_data,
        dewormed ~ (assigned_treatment*assigned_dist_group | sms_treatment), 
        family = bernoulli(link = "probit"),
        prior = manual_priors,
	control = list(adapt_delta = 0.99)
    )
}

if (script_options$save_fit) {
    saveRDS(
        rf_fit,
        file.path(
            script_options$fit_path,
            script_options$fit_file
        )
    )
}


rf_priors = brm(
    data = analysis_data,
        dewormed ~ (assigned_treatment*assigned_dist_group | sms_treatment), 
        family = bernoulli(link = "probit"),
        sample_prior = "only", 
        prior = manual_priors
)





breakdown_sms_comp = comparisons(
    rf_fit, 
    variables = "assigned_treatment",
    newdata = datagrid(
        assigned_dist_group = unique(analysis_data$assigned_dist_group), 
        sms_treatment = unique(analysis_data$sms_treatment))
) %>%
    as_tibble()

create_comp_dfs = function(fit, interval) {
    close_far_comp_df = comparisons(
        fit,
        variable = "assigned_treatment",
        newdata = datagrid(
            assigned_dist_group = unique(analysis_data$assigned_dist_group),
            sms_treatment = unique(analysis_data$sms_treatment)), 
        conf.level = interval
    ) 

    new_low = paste0("conf.low_", (1 - interval)/2)
    new_high = paste0("conf.high_", 1 - (1 - interval)/2)
    comp_df = bind_rows(
        close_far_comp_df
    )  %>%
    mutate(interval = interval) %>%
    as_tibble()
    return(comp_df)

}


comp_df = comparisons(
        rf_fit,
        variable = list(sms_treatment = "all"),
        newdata = datagrid(
            assigned_dist_group = unique(analysis_data$assigned_dist_group),
            assigned_treatment = unique(analysis_data$assigned_treatment)
            ), 
        conf.level = 0.95
    ) 


sms_comp_df = comparisons(
        rf_fit,
        variable = list(assigned_treatment = "all"),
        newdata = datagrid(
            assigned_dist_group = unique(analysis_data$assigned_dist_group),
            sms_treatment = c("smscontrol", "socialinfo")
            ), 
        conf.level = 0.95
    ) 

make_quantiles = function(y, data) {
    q_data = data %>%
        summarise(
            mean_est = mean({{ y }}), 
            per_0.025 = quantile({{ y }}, 0.025),
            per_0.05 = quantile({{ y }}, 0.05), 
            per_0.1 = quantile({{ y }}, 0.1), 
            per_0.25 = quantile({{ y }}, 0.25), 
            per_0.5 = quantile({{ y }}, 0.5), 
            per_0.75 = quantile({{ y }}, 0.75), 
            per_0.9 = quantile({{ y }}, 0.9),
            per_0.95 = quantile({{ y }}, 0.95),
            per_0.975 = quantile({{ y }}, 0.975)
        )
    return(q_data)
}


sms_comp_draws = posteriordraws(sms_comp_df) %>%
    as_tibble()

comp_draws = posteriordraws(comp_df) %>%
    as_tibble()


sms_summ_close_far = sms_comp_draws %>%
    group_by(
        contrast, 
        sms_treatment, 
        assigned_dist_group
    ) %>%
    make_quantiles(y = draw)


sms_summ_combined = sms_comp_draws %>%
    group_by(
        contrast, 
        sms_treatment
    ) %>%
    make_quantiles(y = draw) %>%
    mutate(assigned_dist_group = "combined")

comp_summ_close_far = comp_draws %>%
    group_by( 
        contrast,
        assigned_treatment, 
        assigned_dist_group
    ) %>%
    make_quantiles(y = draw) 

comp_summ_combined = comp_draws %>%
    group_by( 
        contrast,
        assigned_treatment
    ) %>%
    make_quantiles(y = draw) %>%
    mutate(assigned_dist_group = "combined")

comp_summ_df = bind_rows(
    comp_summ_close_far,
    comp_summ_combined
)  %>%
    mutate(
        contrast_left = str_extract(contrast, "^\\w+"),
        contrast_right = str_extract(contrast, "\\w+$")
    ) %>%
    ungroup() 

sms_comp_summ_df = bind_rows(
    sms_summ_close_far,
    sms_summ_combined
) %>%
    mutate(
        contrast_left = str_extract(contrast, "^\\w+"),
        contrast_right = str_extract(contrast, "\\w+$")
    )  %>%
    ungroup() 



 
clean_sms_comp_draws = sms_comp_draws %>%
    bind_rows(sms_comp_draws %>% mutate(assigned_dist_group = "combined")) %>%
    mutate(
        contrast_left = str_extract(contrast, "^\\w+"),
        contrast_right = str_extract(contrast, "\\w+$")
    )  %>%
    ungroup()  %>%
    mutate(
        assigned_dist_group = factor(assigned_dist_group, levels = c("combined", "close", "far")), 
        contrast_left = factor(contrast_left, levels = c("ink", "calendar", "bracelet"))
    ) %>%
    filter(
        contrast_right == "control", 
        contrast_left != "control"
        ) %>%
    mutate(assigned_treatment = contrast_left) %>%
    mutate(assigned_treatment = fct_relabel(assigned_treatment, str_to_title)) %>%
    mutate(assigned_dist_group = factor(assigned_dist_group, levels = c("combined", "close", "far" ))) 

sms_comp_te_diff_df = clean_sms_comp_draws %>%
    select(
         draw,  
        assigned_treatment, assigned_dist_group, sms_treatment
    ) %>%
    pivot_wider(
        id_cols = c(assigned_treatment, assigned_dist_group),
        names_from = sms_treatment, 
        values_from = draw, 
        values_fn = list
    )  %>%
    unnest(c(socialinfo, smscontrol)) %>%
    mutate(
        te_diff = (socialinfo - smscontrol)
    )  %>%
    group_by(
        assigned_treatment, 
        assigned_dist_group
    ) %>%
    make_quantiles(y = te_diff)

pooled_sms_comp_te_diff_df = clean_sms_comp_draws %>%
    select(
         draw,  
        assigned_treatment, assigned_dist_group, sms_treatment
    ) %>%
    pivot_wider(
        id_cols = c(assigned_treatment, assigned_dist_group),
        names_from = sms_treatment, 
        values_from = draw, 
        values_fn = list
    )  %>%
    unnest(c(socialinfo, smscontrol)) %>%
    mutate(
        te_diff = (socialinfo - smscontrol)
    )  %>%
    make_quantiles(y = te_diff)

plot_single_sms_est = function(sms_df, 
                                color_var,
                                top_title = NULL, 
                                width = 0.3, 
                                crossbar_width = 0.2, 
                                vline = TRUE) {
  pos_dodge <- position_dodge(width = width)
  sms_plot = 
      sms_df %>% 
        ggplot(aes(
            y = assigned_treatment, 
            group = interaction(assigned_dist_group, sms_treatment))) +
        geom_linerange(aes(xmin = per_0.05, xmax = per_0.95, color = {{color_var}}), position = pos_dodge, size = 0.3) +
        geom_crossbar(aes(x = per_0.5, xmin = per_0.1, xmax = per_0.9, color = {{ color_var }}), position = pos_dodge, fatten = 2, size = 0.4, width = crossbar_width) +
        geom_linerange(aes(xmin = per_0.25, xmax = per_0.75, color = {{ color_var }}), position = pos_dodge, alpha = 0.4, size = 2.25) +
        geom_point(aes(x = per_0.5, color = {{ color_var }}), position = pos_dodge, size = 1.8) +
        geom_point(aes(x = per_0.5), position = pos_dodge, color = "white", size = 0.6) +
        scale_y_discrete(drop = FALSE) +
        scale_color_canva("", labels = str_to_title, palette = canva_palette_vibrant) + 
        labs(
          subtitle = "",
          x = "", y = "") +
        theme_minimal() +
        theme(
          legend.position = "bottom"
        ) + 
        NULL

  if (vline == TRUE) {
    sms_plot = sms_plot +
         geom_vline(xintercept = 0, linetype = "dotted") 
  }
        
  return(sms_plot)
} 

sms_comp_summ_df = sms_comp_summ_df %>%
    mutate(
        assigned_dist_group = factor(assigned_dist_group, levels = c("combined", "close", "far")), 
        contrast_left = factor(contrast_left, levels = c("ink", "calendar", "bracelet")), 
    )


sms_comp_summ_df  %>%
    filter(
        contrast_right == "control", 
        contrast_left != "control"
        ) %>%
    mutate(assigned_treatment = contrast_left) %>%
    mutate(assigned_treatment = fct_relabel(assigned_treatment, str_to_title)) %>%
    mutate(assigned_dist_group = factor(assigned_dist_group, levels = c("combined", "close", "far" ))) %>%
    mutate(sms_treatment = case_when(
        sms_treatment == "socialinfo" ~ "Social info", 
        sms_treatment == "smscontrol" ~ "SMS control"
    )) %>%
    plot_single_sms_est(
        color_var = sms_treatment
        ) +
    facet_wrap(~assigned_dist_group, ncol = 1, labeller = as_labeller(str_to_title)) 

ggsave(
    file.path(
        script_options$output_path,
        "TE-by-sms-incentive.pdf"
    ),
    width = 8,
    height = 6,
    dpi = 500
)


sms_comp_te_diff_df %>%
    mutate(assigned_dist_group = factor(assigned_dist_group, levels = c("combined", "close", "far" ))) %>%
    mutate(sms_treatment = NA) %>%
    plot_single_sms_est(
        color_var = NULL
        ) +
    facet_wrap(~assigned_dist_group, ncol = 1, labeller = as_labeller(str_to_title)) 




ggsave(
    file.path(
        script_options$output_path,
        "TE-diff-by-sms-incentive.pdf"
    ),
    width = 8,
    height = 6,
    dpi = 500
)


#### LOO ####
# rf_no_sms = brm(
#     data = analysis_data,
#         dewormed ~ assigned_treatment*assigned_dist_group, 
#         family = bernoulli(link = "probit")
# )

# loo_fit = loo(rf_no_sms, rf_fit)



#### Freq ####
library(fixest)
freq_sms = feglm(
    data = analysis_data %>% filter(sms_treatment != "reminderonly") %>% mutate(sms_treatment = fct_drop(sms_treatment)), 
        dewormed ~ 0 + assigned_treatment:assigned_dist_group:sms_treatment, 
        family = binomial(link = "probit"), 
        cluster = ~county
)

tidy_freq_sms = tidy(freq_sms)

create_hyp_matrix = function(tidy_fit, treat_z, treat_d) {
    treat_df = tidy_fit %>%
        select(term) %>%
        mutate(
            treatment = str_extract(term, "(?<=treatment).*(?=:assigned_dist)"), 
            dist = str_extract(term, "(?<=assigned_dist_group).*(?=:sms_treatment)"), 
            sms = str_extract(term, "(?<=sms_treatment).*$")
        ) %>%
        select(
            -term
        ) 

    hyp_m = treat_df %>%
        mutate(
            hyp_m = 0
        ) %>%
        # control diff
        mutate(
            hyp_m = if_else(dist == treat_d & treatment == "control", -1, hyp_m)
        ) %>%
        mutate(
            hyp_m = if_else(dist == treat_d & treatment == treat_z, 1, hyp_m)
        ) %>%
        pull(hyp_m)
    return(hyp_m)
}

hyp_grid = expand_grid(
    treat_z = c("ink", "calendar", "bracelet"), 
    treat_d = c("close", "far")
)

hyp_ms = map2(
    hyp_grid$treat_z,
    hyp_grid$treat_d,
    ~create_hyp_matrix(tidy_freq_sms, .x, .y)
) 


hyp_matrix = hyp_ms %>%
    unlist() %>%
    matrix(byrow = TRUE, ncol = 16)

hyp_matrix

# computationally singular warning :(
# freq_sms_diff_test = car::linearHypothesis(freq_sms, hyp_matrix)

hyp_pvals = map(
    hyp_ms, 
    ~car::linearHypothesis(freq_sms, .x)
    ) %>%
    map("Pr(>Chisq)") %>%
    map_dbl(~.x[2])  
hyp_grid$p_value = hyp_pvals 

hyp_grid = hyp_grid %>%
    mutate(
        p_adj = p.adjust(p_value),
        signif_5pct = p_value < 0.05, 
        signif_1pct = p_value < 0.01
    )

hyp_grid %>%
    write_csv(
    file.path(
        script_options$output_path,
        "freq-probit-sms-hypothesis-tests.csv"
    )
)



#### Old Stuff ####
# comp_summ_df %>%
#     filter(
#         contrast_right == "smscontrol",
#         contrast_left != "smscontrol",
#         (contrast_left != "reminderonly" | assigned_treatment == "control")
#     ) %>%
#     mutate(
#         sms_treatment = if_else(str_detect(contrast, "reminder"), "Reminder Only", "Social Info")
#     ) %>%
#     mutate(assigned_treatment = fct_relabel(assigned_treatment, str_to_title)) %>%
#     mutate(assigned_dist_group = factor(assigned_dist_group, levels = c("combined", "close", "far" )) %>% fct_rev) %>%
#     plot_single_sms_est(
#         color_var = assigned_dist_group
#         ) +
#     labs(
#         title = "SMS Treatment Effect By Incentive and Distance Condition"
#     ) +
#     facet_wrap(~sms_treatment) 

# comp_summ_df %>%
#     mutate(
#         contrast_left = str_extract(contrast, "^\\w+"),
#         contrast_right = str_extract(contrast, "\\w+$")
#     ) %>%
#     filter(
#         contrast_right == "smscontrol",
#         contrast_left != "smscontrol",
#         (contrast_left != "reminderonly" | assigned_treatment == "control")
#     ) %>%
#     mutate(
#         sms_treatment = if_else(str_detect(contrast, "reminder"), "Reminder Only", "Social Info")
#     ) %>%
#     mutate(assigned_treatment = fct_relabel(assigned_treatment, str_to_title)) %>%
#     mutate(assigned_dist_group = factor(assigned_dist_group, levels = c("combined", "close", "far" )) %>% fct_rev) %>%
#     plot_single_sms_est(
#         color_var = assigned_dist_group
#         ) +
#     labs(
#         title = "SMS Treatment Effect By Incentive and Distance Condition"
#     ) +
#     facet_wrap(~sms_treatment) 

# ggsave(
#     file.path(
#         script_options$output_path,
#         "sms-TE-by-dist-incentive.png"
#     ),
#     width = 7.5,
#     height = 5,
#     dpi = 500
# )

# create_pred_dfs = function(fit, interval){
#     combined_pred_df = predictions(
#         fit,
#         newdata = datagrid(
#             assigned_treatment = unique(analysis_data$assigned_treatment), 
#             assigned_dist_group = unique(analysis_data$assigned_dist_group),
#             sms_treatment = unique(analysis_data$sms_treatment)), 
#         by = c("assigned_treatment", "sms_treatment"), 
#         conf.level = interval
#     ) %>% 
#         mutate(assigned_dist_group = "combined")
#     close_far_pred_df = predictions(
#         fit,
#         newdata = datagrid(
#             assigned_treatment = unique(analysis_data$assigned_treatment), 
#             assigned_dist_group = unique(analysis_data$assigned_dist_group),
#             sms_treatment = unique(analysis_data$sms_treatment)), 
#         conf.level = interval
#     )

#     new_low = paste0("conf.low_", (1 - interval)/2)
#     new_high = paste0("conf.high_", 1 - (1 - interval)/2)
#     pred_df = bind_rows(
#         combined_pred_df, 
#         close_far_pred_df
#     )  %>%
#     mutate(interval = interval) %>%
#     as_tibble()
#     return(pred_df)
# }


# prior_df = map_dfr(
#     c(0.9, 0.8, 0.5), 
#     create_pred_dfs, fit = rf_priors) %>%
#     mutate(fit_type = "prior-predict")
# fit_pred_df = map_dfr(
#     c(0.9,
#       0.8,
#       0.5), 
#     create_pred_dfs, fit =  rf_fit) %>%
#     mutate(fit_type = "fit")
# pred_df = bind_rows(
#     fit_pred_df,
#     # prior_df
# ) %>%
#     mutate(fit_type = factor(fit_type, levels = c("prior-predict", "fit")))


# pred_df = pred_df %>%
#     mutate(
#         assigned_dist_group = factor(assigned_dist_group, levels = c("combined", "close", "far")), 
#         assigned_dist_group = fct_relabel(assigned_dist_group, str_to_title), 
#         assigned_treatment = fct_relabel(assigned_treatment, str_to_title)
#     )

# plot_brm_estimands = function(fit_data, 
#                               pos_height = 0.8, 
#                               color_var,
#                               center_bar_size = 3, 
#                               top_levels = c("Bracelet", "Combined")) {

#     plot_pos <- position_dodge2(pos_height, reverse = FALSE)
#     # plot_pos = ggstance::position_dodgev(height = 0.8)
#     # spacing/height-width ratio remains the same for other plots
#     fit_data = fit_data %>%
#         mutate(across(c(estimate, conf.low, conf.high),
#                       ~if_else(
#                        fct_match(fit_type, "prior-predict") & !(fct_match(assigned_treatment, top_levels[1]) & fct_match(assigned_dist_group, top_levels[2])),
#                         NA_real_,
#                         .x
#                         ))) %>%
#         mutate( 
#             fit_alpha = if_else(
#                 fct_match(fit_type, "prior-predict"), 
#                 0.1,
#                 1
#             )
#         ) %>%
#         mutate(model = "ed")

#     p = fit_data %>%
#         ggplot(aes(
#             x = estimate, 
#             xmin = conf.low, 
#             xmax = conf.high, 
#             y = assigned_treatment, 
#             group = model)) +
#       NULL

#     p = p  +
#         geom_linerange(
#             data = . %>% filter(interval == 0.5), 
#             aes(color = {{ color_var }}),
#             alpha = 0.4, 
#             size = center_bar_size,
#             position = plot_pos
#         ) +
#         geom_crossbar(
#             data = . %>% filter(interval == 0.8), 
#             aes(color = {{ color_var }}),
#             fatten = 0, 
#             size = 0.4, 
#             width = 0.8,
#             position = plot_pos
#         )  +
#         geom_linerange(
#             data = . %>% filter(interval == 0.9), 
#             aes(color = {{ color_var }}),
#             size = 0.4, 
#             position = plot_pos
#         ) +
#         geom_point(
#             data = . %>% filter(interval == 0.9),
#             aes(color = {{ color_var }}),
#             position = plot_pos
#         ) +
#         geom_point(
#             data = . %>% filter(interval == 0.9),
#             color = "white", size = 0.75, position = plot_pos) +
#         theme(legend.position = "top", legend.direction = "vertical") +
#         guides(color = guide_legend(ncol = 3)) +
#         scale_y_discrete("") +
#         labs(
#             x = "",
#         caption = #"Dotted line range: 98% credible interval. 
#                     "Line range: 90% credible interval. 
#                     Outer box: 80% credible interval. Inner box: 50% credible interval. 
#                     Thick vertical line: median. Point: mean."
        
#         )  
#     return(p)
# }









# theme_set(theme_minimal() +
#             theme(legend.position = "bottom"))

# p_social_info_levels = pred_df %>%
#     filter(sms_treatment == "socialinfo") %>%
#     plot_brm_estimands(color_var = fit_type) +
#     geom_vline(
#         xintercept = 0, linetype = "dotted"
#     ) +
#     guides(colour = "none", alpha = "none")  +
#     ggforce::facet_col(vars(assigned_dist_group), 
#                     space = "free",
#                     scales = "free_y") +
#     scale_color_manual("", 
#                         values = c("darkgrey", "black"), 
#                     aesthetics = c("color", "fill"))  
# ggsave(
#     plot = p_social_info_levels,
#     file.path(
#         script_options$output_path,
#         "sms-socialinfo-levels.png" 
#     ),
#   width = 7.5,
#   height = 5,
#   dpi = 500
# )


# p_reminder_levels = pred_df %>%
#     filter(sms_treatment == "smscontrol") %>%
#     plot_brm_estimands(color_var = fit_type) +
#     geom_vline(
#         xintercept = 0, linetype = "dotted"
#     ) +
#     guides(colour = "none", alpha = "none")  +
#     ggforce::facet_col(vars(assigned_dist_group), 
#                     space = "free",
#                     scales = "free_y") +
#     scale_color_manual("", 
#                         values = c("darkgrey", "black"), 
#                     aesthetics = c("color", "fill"))  


# ggsave(
#     plot = p_reminder_levels,
#     filename = file.path(
#         script_options$output_path,
#         "sms-control-levels.png"
#         ),
#     width = 7.5,
#     height = 5,
#     dpi = 500
# )


# p_comp_levels = pred_df %>%
#     filter(fct_match(fit_type, "fit")) %>%
#     filter(
#         sms_treatment != "reminderonly" | (sms_treatment == "reminderonly" & assigned_treatment == "Control")
#     ) %>%
#     mutate(sms_treatment = case_when(
#         sms_treatment == "socialinfo" ~ "Social Info",
#         sms_treatment == "smscontrol" ~ "SMS Control", 
#         sms_treatment == "reminderonly" ~ "Reminder Only"
#     )) %>%
#     plot_brm_estimands(
#         color_var = sms_treatment
#     ) +
#     ggforce::facet_col(vars(assigned_dist_group), 
#                     space = "free",
#                     scales = "free_y") +
#     scale_color_canva("", labels = str_to_title, palette = canva_palette_vibrant) +
#     theme(legend.position =  "bottom")  +
#     geom_vline(xintercept = 0, linetype = "longdash") +
#     scale_x_continuous("", breaks = seq(-1, 1, 0.1), limits = c(0, 0.9))  
# p_comp_levels
# ggsave(
#     plot = p_comp_levels,
#     filename = file.path(
#         script_options$output_path,
#         "sms-comp-levels-axis.png"),
#     width = 7.5,
#     height = 5,
#     dpi = 500
# )


# # ------------------------
