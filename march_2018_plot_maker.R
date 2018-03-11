
# Combined Take-up ----------------------------------------------------------------

plot_takeup(
  takeup_summ_data = est_deworming_takeup_all,
  data_preparer = . %>%
    filter(
      !name_matched,
      sms.treatment.2 == "control",
      incentive_treatment != "Bracelet Social"
      ),
  incentive_treatment_col = "incentive_treatment"
    ) +
  labs(y = "") 

ggsave(filename = "static_level_smsctrl.png", path = "march_2018_plots", 
       width = 6.5, height = 3.25)

plot_takeup(
  takeup_summ_data = est_deworming_takeup_all,
  data_preparer = . %>%
    filter(
      !name_matched,
      incentive_treatment != "Bracelet Social"
      ),
  incentive_treatment_col = "incentive_treatment",
  include_sms_treatment = TRUE
    ) +
  labs(y = "") 

ggsave(filename = "static_level_all.png", path = "march_2018_plots", 
       width = 7.5, height = 4.25)


# Combined ATE ------------------------------------------------------------

est_deworming_takeup_ate_all %>% 
  filter(incentive_treatment_left != "bracelet social",
         !(incentive_treatment_left == "bracelet" & incentive_treatment_right == "calendar")) %>% 
  # mutate(incentive_treatment_left = if_else(incentive_treatment_left == "bracelet" & incentive_treatment_right == "calendar", 
  #                                           "bracelet social",
  #                                           as.character(incentive_treatment_left)) %>% 
  #          incentive_treatment_factor()) %>% 
  plot_ate(
    ate_summ_data = .,
      data_preparer = . %>% 
      filter(incentive_treatment_right %in% c("Control", "Calendar"),
             incentive_treatment_left != incentive_treatment_right,
             sms.treatment.2_left == "No SMS",
             sms.treatment.2_left == sms.treatment.2_right),
    incentive_treatment_left_col = "incentive_treatment_left" 
  ) + 
  labs(y = "") 

ggsave(filename = "static_ate_smsctrl.png", path = "march_2018_plots", 
       width = 7.5, height = 3.25)

# est_deworming_takeup_ate_all %>%
est_deworming_takeup_ate_phone %>%
# est_deworming_takeup_ate_phone_name_matched %>%
  mutate(same_incentive = incentive_treatment_left == incentive_treatment_right,
         same_sms = sms.treatment.2_left == sms.treatment.2_right) %>% 
  filter(phone_owner_left, #!name_matched_left, 
         incentive_treatment_left != "bracelet social",
         !same_incentive & same_sms,
         sms.treatment.2_right != "reminder.only",
         !(incentive_treatment_left == "bracelet" & incentive_treatment_right == "calendar")) %>% 
  plot_ate(
    ate_summ_data = .,
    incentive_treatment_left_col = "incentive_treatment_left",
    include_sms_treatment = TRUE
  ) + 
  labs(y = "") 

ggsave(filename = "static_ate_all.png", path = "march_2018_plots", 
       width = 7.5, height = 4.25)

# Combined Signal Diff by Distance Diff -----------------------------------

est_deworming_takeup_ate %>% 
  filter(incentive_treatment_left %in% c("ink", "bracelet"),
         incentive_treatment_right == "control",
         dist.pot.group_left == dist.pot.group_right,
         sms.treatment.2_left == sms.treatment.2_right,
         sms.treatment.2_left == "control") %>% { 
   inner_join(filter(., dist.pot.group_left == "close"), 
              filter(., dist.pot.group_left == "far"),
              c("incentive_treatment_left", "name_matched_left", "iter_id", "phone_owner_left", "treatment_size"), suffix = c("_close", "_far")) 
 } %>% 
 mutate(iter_ate = iter_ate_close - iter_ate_far) %>% 
 prepare_est_deworming_ate() %>% 
 plot_ate(
   ate_summ_data = .,
   incentive_treatment_left_col = "incentive_treatment_left",
   include_sms_treatment = FALSE
 ) + 
 labs(y = "", subtitle = "Estimating Difference in Social Effect Between Far and Close") 

ggsave(filename = "static_ate_signal_effect_dist_diff.png", path = "march_2018_plots", 
       width = 7.5, height = 3.25)

# Combined Signal Effect by SMS -------------------------------------------

est_deworming_takeup_ate %>% 
  filter(incentive_treatment_left %in% c("ink", "bracelet"),
         incentive_treatment_right == "control",
         phone_owner_left, phone_owner_right, !name_matched_left,
         dist.pot.group_left == dist.pot.group_right,
         sms.treatment.2_left == sms.treatment.2_right) %>% {
   inner_join(filter(., sms.treatment.2_left == "social.info"),
                filter(., sms.treatment.2_left == "control"),
                c("incentive_treatment_left", "dist.pot.group_left", "phone_owner_left", "treatment_size", "iter_id"),
                suffix = c("_sms", "_nosms"))  
 } %>% 
 mutate(iter_ate = iter_ate_sms - iter_ate_nosms) %>% 
 prepare_est_deworming_ate()  %>% 
 plot_ate(
   ate_summ_data = .,
   incentive_treatment_left_col = "incentive_treatment_left",
   include_sms_treatment = FALSE
 ) + 
 labs(y = "", subtitle = "Estimating Difference in Social Effect Between Social Info SMS Treatment and No SMS Treatment") 

ggsave(filename = "static_ate_signal_effect_sms_diff.png", path = "march_2018_plots", 
       width = 7.5, height = 3.25)

# By Distance Take-up -------------------------------------------------------------

plot_takeup(
  takeup_summ_data = est_deworming_takeup_dist,
  data_preparer = . %>%
    filter(
      !name_matched,
      sms.treatment.2 == "control",
      incentive_treatment != "Bracelet Social"
      ) %>% 
    mutate_at(vars(starts_with("dist.pot.group")), funs(fct_relabel(., str_to_title))), 
    incentive_treatment_col = "incentive_treatment"
    ) +
  labs(y = "", subtitle = "Split By Distance") +
  facet_wrap(~ dist.pot.group, ncol = 1)

ggsave(filename = "static_level_smsctrl_dist.png", path = "march_2018_plots", 
       width = 7.5, height = 5.5)

plot_takeup(
  takeup_summ_data = est_deworming_takeup_dist,
  data_preparer = . %>%
    filter(
      !name_matched,
      incentive_treatment != "Bracelet Social"
      ) %>% 
    mutate_at(vars(starts_with("dist.pot.group")), funs(fct_relabel(., str_to_title))), 
    incentive_treatment_col = "incentive_treatment",
    include_sms_treatment = TRUE
  ) +
  labs(y = "", subtitle = "Split By Distance") +
  facet_wrap(~ dist.pot.group, ncol = 1) 

ggsave(filename = "static_level_all_dist.png", path = "march_2018_plots", 
       width = 7.5, height = 7)


# By Distance ATE ---------------------------------------------------------

est_deworming_takeup_ate_dist %>% 
  filter(incentive_treatment_left != "bracelet social",
         !(incentive_treatment_left == "bracelet" & incentive_treatment_right == "calendar")) %>% 
  # mutate(incentive_treatment_left = if_else(incentive_treatment_left == "bracelet" & incentive_treatment_right == "calendar", 
  #                                           "bracelet social",
  #                                           as.character(incentive_treatment_left)) %>% 
  #          incentive_treatment_factor()) %>% 
  plot_ate(
    ate_summ_data = .,
    data_preparer = . %>% 
    filter(incentive_treatment_right %in% c("Control", "Calendar"),
           incentive_treatment_left != incentive_treatment_right,
           dist.pot.group_left == dist.pot.group_right,
           sms.treatment.2_left == "No SMS",
           sms.treatment.2_left == sms.treatment.2_right) %>% 
      mutate_at(vars(starts_with("dist.pot.group")), funs(fct_relabel(factor(.), str_to_title))), 
    incentive_treatment_left_col = "incentive_treatment_left" 
  ) + 
  labs(y = "") +
  facet_wrap(~ dist.pot.group_left, ncol = 1)

ggsave(filename = "static_ate_smsctrl_dist.png", path = "march_2018_plots", 
       width = 7.5, height = 4.75)

est_deworming_takeup_ate_phone_dist %>% 
  mutate(same_incentive = incentive_treatment_left == incentive_treatment_right,
         same_sms = sms.treatment.2_left == sms.treatment.2_right) %>% 
  filter(phone_owner_left,
         incentive_treatment_left != "bracelet social",
         !same_incentive & same_sms,
         dist.pot.group_left == dist.pot.group_right,
         !(incentive_treatment_left == "bracelet" & incentive_treatment_right == "calendar")) %>% 
  # mutate(incentive_treatment_left = if_else(incentive_treatment_left == "bracelet" & incentive_treatment_right == "calendar", 
                                            # "bracelet social",
                                            # as.character(incentive_treatment_left)) %>% 
           # incentive_treatment_factor()) %>% 
  plot_ate(
    ate_summ_data = .,
    incentive_treatment_left_col = "incentive_treatment_left",
    data_preparer = . %>% 
      mutate_at(vars(starts_with("dist.pot.group")), funs(fct_relabel(factor(.), str_to_title))), 
    include_sms_treatment = TRUE,
    sms_treatment_right = FALSE
  ) + 
  labs(y = "") +
  facet_wrap(~ dist.pot.group_left, ncol = 1)

ggsave(filename = "static_ate_all_dist.png", path = "march_2018_plots", 
       width = 9, height = 7.5)

# Dynamic -----------------------------------------------------------------

est_deworming_day_all %>% 
  mutate(full_observation = incentive_treatment_static == incentive_treatment_dynamic) %>% 
  filter(full_observation,
         incentive_treatment_static != "bracelet social") %>% 
  plot_dyn_takeup_daily(control_observation = F) +
  coord_cartesian(ylim = c(0, 0.09)) +
  theme(legend.position = "bottom")

ggsave(filename = "dyn_daily_all.png", path = "march_2018_plots", 
       width = 7.5, height = 6.5)

est_deworming_day_dist %>% 
  mutate(full_observation = incentive_treatment_static == incentive_treatment_dynamic) %>% 
  filter(full_observation,
         incentive_treatment_static != "bracelet social",
         incentive_treatment_static == incentive_treatment_dynamic,
         dist.pot.group_static == dist.pot.group_dynamic) %>% 
  mutate_at(vars(starts_with("dist.pot.group_static")), funs(fct_relabel(., str_to_title))) %>% 
  plot_dyn_takeup_daily() +
  facet_wrap(~ dist.pot.group_static, nrow = 1) +
  coord_cartesian(ylim = c(0, 0.125)) +
  theme(legend.position = "bottom")

ggsave(filename = "dyn_daily_all_dist.png", path = "march_2018_plots", 
       width = 8.5, height = 6)
