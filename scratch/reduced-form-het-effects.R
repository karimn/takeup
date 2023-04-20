library(tidyverse)
library(fixest)
library(broom)
library(marginaleffects)

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


analysis_data <- monitored_nosms_data %>% 
    mutate(
    assigned_treatment = assigned.treatment, 
    assigned_dist_group = dist.pot.group)

social_perception_data = baseline.data %>% 
  select(matches("^(praise|stigma)_[^_]+$")) %>% 
  gather(key = key, value = response) %>% 
  separate(key, c("praise.stigma", "topic"), "_") %>% 
  separate(topic, c("topic", "question.group"), -2) %>% 
  filter(!is.na(response))  %>%
  count(praise.stigma, topic, response) %>% 
  group_by(praise.stigma, topic) %>% 
  mutate(n = n/sum(n)) %>%
  mutate_at(vars(praise.stigma, response), ~ fct_relabel(factor(.), str_to_title)) %>% 
  mutate(topic = fct_recode(factor(topic), 
                            "Wearing/not wearing nice clothes to church" = "clothe",
                            "Use Latrine/open defecation" = "defecat",
                            "Deworming/not deworming during MDA" = "dewor",
                            "Immunize/not immunize children" = "immuniz")) 


externality_tab_df = baseline.data %>%
    select(worms_affect, neighbours_worms_affect) %>%
    group_by(worms_affect, neighbours_worms_affect) %>%
    summarise(
        n = n()
    )  %>%
    filter(!is.na(worms_affect)) %>%
    spread(
        worms_affect, 
        n
    ) %>%
    rename(
        "(neighbour)/worms_affect" = neighbours_worms_affect
    )  


baseline.data = baseline.data %>%
    mutate(
        any_externality_knowledge = worms_affect == "yes" | neighbours_worms_affect == "yes"
    )


baseline.data %>%
    group_by(
        cluster.id
    ) %>%
    summarise(
        frac_externality_knowledge = mean(any_externality_knowledge, na.rm = TRUE), 
        treatment = unique(assigned.treatment)
    ) %>%
    ggplot(aes(
        x = frac_externality_knowledge
    )) +
    geom_histogram(
        colour = "black"
    ) +
    theme_minimal() + 
    labs(
        title = "Fraction Externality Knowledge Per Community", 
        x = "Fraction Aware of Worm Externality"
    )
ggsave(
    "temp-plots/hist-frac-aware-externalities.pdf", 
    width = 8, 
    height = 6
)


#### Externalities ####

cluster_clean_externality_df  = baseline.data %>%
    group_by(
        cluster.id
    ) %>%
    summarise(
        frac_externality_knowledge = mean(any_externality_knowledge, na.rm = TRUE)
    ) %>%
    inner_join(
        analysis_data %>%
            group_by(cluster.id) %>%
            summarise(
                takeup = mean(dewormed), 
                assigned_treatment = unique(assigned_treatment), 
                assigned_dist_group = unique(assigned_dist_group)
            ) %>%
            select(
                cluster.id, 
                assigned_treatment, 
                assigned_dist_group, 
                takeup
            ), 
    by = "cluster.id"
    )


clean_externality_df = analysis_data %>%
    left_join(
        baseline.data %>%
            group_by(
                cluster.id
            ) %>%
            summarise(
                frac_externality_knowledge = mean(any_externality_knowledge, na.rm = TRUE)
            ) %>% 
            ungroup() %>%
            mutate(
                above_med_externality_knowledge = frac_externality_knowledge > median(frac_externality_knowledge, na.rm = TRUE)
            ),
    by = "cluster.id"
    )




#### Perceptions ####
clean_perception_data = baseline.data %>% 
  select(cluster.id, matches("^(praise|stigma)_[^_]+$")) %>% 
  gather(key = key, value = response, -cluster.id) %>% 
  separate(key, c("praise.stigma", "topic"), "_") %>% 
  separate(topic, c("topic", "question.group"), -2)  %>%
  filter(!is.na(response))  




long_topic_judgement_score_df = clean_perception_data %>%
        count(cluster.id, praise.stigma, topic, response)  %>%
        group_by(cluster.id, praise.stigma, topic) %>% 
        mutate(n = n/sum(n))   %>%
        group_by(topic, cluster.id) %>%
        summarise(
            judge_score = n[praise.stigma == "praise" & response == "yes"] * n[praise.stigma == "stigma" & response == "yes"],
            stigma_score = n[praise.stigma == "stigma" & response == "yes"],
            praise_score = n[praise.stigma == "praise" & response == "yes"]
        )  %>%
        ungroup() %>%
        group_by(topic) %>%
        mutate(
            median_judge_score = median(judge_score, na.rm = TRUE),
            median_stigma_score = median(stigma_score, na.rm = TRUE)
        )



overall_judgement_score_df = clean_perception_data %>%
  count(cluster.id, praise.stigma, topic, response)  %>%
  group_by(cluster.id, praise.stigma, topic) %>% 
  mutate(n = n/sum(n))   %>%
  group_by(cluster.id) %>%
  filter(response == "yes") %>%
  summarise(
    judge_score = prod(n), 
  ) %>% 
  ungroup() %>%
  mutate(
    median_judge_score = median(judge_score, na.rm = TRUE)
  ) %>%
  mutate(topic = "overall")

judgement_score_df = bind_rows(
    long_topic_judgement_score_df, 
    overall_judgement_score_df

)

judgement_score_df %>%
    filter(topic == "dewor") %>%
    ggplot(aes(
        x = judge_score, 
    )) +
    geom_histogram(colour = "black") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
        title = "Deworming 'Judgemental' Score Distribution"
    )

ggsave(
    "temp-plots/deworm-judgemental-score-histogram.pdf", 
    width = 10, 
    height = 10
)



wide_judgement_score_df = judgement_score_df %>%
    pivot_wider(
        names_from = topic, 
        values_from = c(judge_score, stigma_score, praise_score, median_judge_score, median_stigma_score)
    )


judge_analysis_data = left_join(
    analysis_data, 
    wide_judgement_score_df,
    by = "cluster.id"
) 

judge_analysis_data = judge_analysis_data %>%
    mutate(
        judgemental_dewor = judge_score_dewor >= median_judge_score_dewor,
        judgemental_immuniz = judge_score_immuniz >= median_judge_score_immuniz,
        judgemental_defecat = judge_score_defecat >= median_judge_score_defecat,
        judgmental_overall = judge_score_overall >= median_judge_score_overall,
    )

#### Community Centrality ####

analysis_data %<>% 
  nest_join(
    endline.know.table.data %>% 
      filter(fct_match(know.table.type, "table.A")),
    by = "KEY.individ", 
    name = "knowledge_data"
  ) %>% 
  mutate(
    map_dfr(knowledge_data, ~ {
      tibble(
        obs_know_person = sum(.x$num.recognized),
        obs_know_person_prop = mean(.x$num.recognized),
        knows_other_dewormed = sum(fct_match(.x$dewormed, c("yes", "no")), na.rm = TRUE),
        knows_other_dewormed_yes = sum(fct_match(.x$dewormed, "yes"), na.rm = TRUE),
        knows_other_dewormed_no = sum(fct_match(.x$dewormed, "no"), na.rm = TRUE),
        thinks_other_knows = sum(fct_match(.x$second.order, c("yes", "no")), na.rm = TRUE),
        thinks_other_knows_yes = sum(fct_match(.x$second.order, "yes"), na.rm = TRUE),
        thinks_other_knows_no = sum(fct_match(.x$second.order, "no"), na.rm = TRUE),
      )
    }
  ))


comp_belief_data = analysis_data %>%
    filter(obs_know_person > 0)  %>%
    select(KEY.individ, contains("know"), assigned.treatment, dist.pot.group, assigned_dist_group) %>%
    mutate(
        doesnt_know_other_dewormed = obs_know_person - knows_other_dewormed, 
        doesnt_think_other_knows = obs_know_person - thinks_other_knows
    ) %>% 
    select(KEY.individ, 
           assigned.treatment,
           assigned_dist_group,
           obs_know_person,
           knows_other_dewormed_yes,
           knows_other_dewormed_no,
           doesnt_know_other_dewormed, 
           thinks_other_knows_yes, 
           thinks_other_knows_no, 
           doesnt_think_other_knows
           ) %>%
    gather(variable, value, 
        knows_other_dewormed_yes:doesnt_think_other_knows)   %>%
    mutate(knowledge_type = case_when(
        str_detect(variable, "_yes") ~ "yes",
        str_detect(variable, "_no") ~ "no",
        str_detect(variable, "doesnt") ~ "doesn't know"
    )) %>%
    mutate(belief_type = if_else(str_detect(variable, "think"), "2ord", "1ord")) %>%
    mutate(prop = value/obs_know_person) 


comp_belief_data  %>%
    select(KEY.individ, obs_know_person) %>%
    ggplot(aes(
        x = obs_know_person
    )) +
    geom_bar() +
    theme_minimal() +
    labs(
        x = "Number of People Recognised", 
        y = "Count"
    ) +
    scale_x_continuous(breaks = 1:10)

ggsave(
    "temp-plots/num-recognised-hist.pdf", 
    width = 8,
    height = 6
)





#### Regressions ####
probit_fit = analysis_data %>%
    feglm(
        dewormed ~ 0 + assigned_treatment:assigned_dist_group, 
        data = ., 
        family = binomial(link = "probit"), 
        cluster = ~county
    ) 



externality_fit = clean_externality_df %>%
    feglm(
        dewormed ~  frac_externality_knowledge + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit")
    )

judgement_geo_fit = judge_analysis_data %>%
    feglm(
        dewormed ~  judge_score_dewor + assigned_treatment:assigned_dist_group,
        family = binomial(link = "probit")
    )

etable(
    judgement_geo_fit,
    externality_fit, 
    order = "!assigned", 
    keep = c("!assigned"),
    drop = "Constant",
    cluster = ~county
)

