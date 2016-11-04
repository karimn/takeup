# Setup -------------------------------------------------------------------

library(plyr)
library(dplyr)
library(magrittr)
library(purrr)
library(readr)
library(lubridate)
library(httr)
library(stringr)
library(jsonlite)

config <- yaml::yaml.load_file("local_config.yaml")

# arthur.phone.number <- "+254728401235"
# carol.phone.number <- "+254707012589"

pot.desciption.wave1 <- read_csv("pot_description_wave1.csv") %>% 
  rename(desciption.english = English,
         description.kiswahili = Kiswahili)

pot.desciption.wave2 <- read_csv("pot_description_wave2.csv") %>% 
  select(cluster.id, starts_with("description"))

# SMS Script --------------------------------------------------------------

reminder.message <- "Matibabu ya minyoo ya bure iko %s. Jumbe zote ni bure."

amount.dict <- list("no" = "hakuna",
                    "few" = "wachache",
                    "almost half" = "karibu nusu",
                    "half" = "nusu",
                    "more than half" = "zaidi ya nusu",
                    "almost all" = "karibu wote",
                    "all" = "wote")

amount.dict <- tribble(~ en, ~ sw, ~ lower, ~ upper,
                       "no", "hakuna", 0, 0,
                       "few", "wachache", 0, 2,
                       "almost half", "karibu nusu", 2, 4,
                       "half", "nusu", 4, 6,
                       "more than half", "zaidi ya nusu", 6, 8,
                       "almost all",  "karibu wote", 8, 10,
                       "all", "wote", 10, 10)

# social.info.message <- "Matibabu ya minyoo ya bure iko %s, hakuna/wachache/karibu nusu/Nusu/zaidi ya nusu/Karibu wote/wote wa kijiji chako walikuja, hao ni 10 kwa 10 watu wazima.Jibu kwa nambari 111 kusimamisha ujumbe. jumbe zote ni bure"

social.info.message <- "Matibabu ya minyoo ya bure iko %s, %s wa kijiji chako walikuja, hao ni karibu %d kwa 10 watu wazima. Jumbe zote ni bure."

# Send SMS ----------------------------------------------------------------

send.message <- function(phone.numbers, msg, api.url = "http://api.africastalking.com/version1/messaging") {
  POST(api.url,
       accept_json(),
       add_headers(Apikey = config$africastalking_api_key),
       body = list(username = "karimn",
                   from = "EvidenceAct",
                   to = paste(phone.numbers, collapse = ","),
                   message = msg),
       encode = "form",
       verbose())
}

cluster.reminder.messenger <- . %>%
  group_by(cluster.id) %>%
  do(response = send.message(phone.numbers = .$phone,
                             msg = first(.$cluster.msg))) # first(.$alt_name))
                             # api.url = "https://nth-fort-117516.appspot.com/sms/test_sending"))

mid.intervention.messaging <- function(.data) { 
  bind_rows(filter(.data, sms.treatment == "reminder.only") %>%
              group_by(sms.treatment, cluster.id) %>%
              mutate(cluster.msg = sprintf(reminder.message, description.kiswahili)) %>%
              do(response = send.message(phone.numbers = .$phone, 
                                         # api.url = "https://nth-fort-117516.appspot.com/sms/test_sending",
                                         msg = first(.$cluster.msg))) %>% 
              ungroup,
            filter(.data, sms.treatment == "social.info") %>%
              group_by(sms.treatment, cluster.id, village) %>%
              mutate(cluster.msg = sprintf(social.info.message, description.kiswahili, takeup.prop.grp, takeup.prop.10)) %>%
              do(response = send.message(phone.numbers = .$phone, 
                                         # api.url = "https://nth-fort-117516.appspot.com/sms/test_sending",
                                         msg = first(.$cluster.msg))) %>%  
              ungroup)
}
# response <- POST("http://api.africastalking.com/version1/messaging",
#                  accept_json(),
#                  add_headers(Apikey = config$africastalking_api_key),
#                  body = list(username = "karimn",
#                              from = "EvidenceAct",
#                              to = carol.phone.number,
#                              message = "Hi! I might send you a couple of test messages. I'll try to keep it to a minimum."),
#                  encode = "form",
#                  verbose())

# Reward pilot sample ------------------------------------------------------------

reward.pilot <- sms.recruit.data %>% 
  filter(cluster.id %in% c(277, 491, 492)) %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13,
         pot.cluster.id = ifelse(cluster.id == 277, 503, cluster.id)) %>% 
  left_join(select(wide.sms.status.data, to, total, Success), c("phone" = "to")) %>% 
  rename(sms.status.success = Success,
         sms.status.total = total) %>% 
  filter(sms.status.success >= 5) %>% 
  sample_frac(1)


# Sending out reward ------------------------------------------------

reward.msg.sender <- . %>% 
  do(response = send.message(.$phone,
                             msg = "Asante kukubali ujumbe mfupi kutoka Evidence Action, umeshinda muda wa maongezi 50Ksh. Tuma ujumbe 123 kwa 20880 ili upokee zawadi hii sasa hivi. Kutuma ujumbe ni bure."))
                             # msg = "Asante kwa kukubali kupokea jumbe zetu fupi kutoka Evidence Action umejishindia kadi ya simu ya 50 KSh tuma ujumbe 123 kwa nambari 20880 ili kupokea zawadi yako."))

reward.reminder.msg.sender <- . %>% 
  do(response = send.message(.$phone,
                             msg = "Asante kukubali ujumbe mfupi kutoka Evidence Action. Usisahau kuitisha muda wako maongezi wa 50Ksh. Tuma ujumbe 123 kwa 20880 ili upokee zawadi hii sasa hivi. Kutuma ujumbe ni bure."))

# reward.response.1 <- reward.pilot %>% 
#   filter(cluster.id == 277) %>% 
#   slice(1:5) %>% 
#   reward.pilot.msg.sender

# reward.response.2 <- reward.pilot %>% 
#   filter(cluster.id == 491) %>% 
#   slice(1:5) %>% 
#   reward.pilot.msg.sender

# reward.response.3 <- reward.pilot %>% 
#   filter(cluster.id == 492) %>% 
#   slice(2:5) %>% 
#   reward.pilot.msg.sender

# Send Airtime ------------------------------------------------------------

pilot.rewards <- c("+254710398102", "+254712936051", "+254702761933", "+254795038104")

# pilot.airtime.response.1 <- POST("http://api.africastalking.com/version1/airtime/send",
#                  accept_json(),
#                  add_headers(Apikey = config$africastalking_api_key),
#                  body = list(username = "karimn",
#                              recipients = toJSON(list(list(phoneNumber = "+254795038104",
#                                                            amount = "KES 50")),
#                                                  auto_unbox = TRUE),
#                              from = "EvidenceAct"),
#                  encode = "form",
#                  verbose())

# First Reminder for All (Wave 1) -----------------------------------------

wave1.sms.subjects <- sms.recruit.data %>% 
  filter(county != "Kakamega", cluster.id != 1) %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13,
         pot.cluster.id = ifelse(cluster.id == 277, 503, cluster.id)) %>% 
  # left_join(select(known.pot.locations, cluster.id, alt_name), "cluster.id") %>% 
  left_join(pot.desciption.wave1, c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster.msg = sprintf(reminder.message, description.kiswahili)) %>% 
  ungroup


# first.reminder.response.wave1.batch1 <- wave1.sms.subjects %>% 
#   filter(cluster.id == 277) %>% 
#   cluster.reminder.messenger
#   
# first.reminder.response.wave1.batch2 <- wave1.sms.subjects %>% 
#   filter(cluster.id %in% c(491, 492)) %>% 
#   cluster.reminder.messenger
# 
# first.reminder.response.wave1.batch3 <- wave1.sms.subjects %>% 
#   filter(cluster.id == 147) %>% 
#   cluster.reminder.messenger
# 
# first.reminder.response.wave1.batch4 <- wave1.sms.subjects %>%
#   filter(!cluster.id %in% c(27, 491, 492, 147)) %>%
#   cluster.reminder.messenger
#  
# first.reminder.response.wave1 <- bind_rows(first.reminder.response.wave1.batch1, 
#                                            first.reminder.response.wave1.batch2,
#                                            first.reminder.response.wave1.batch3,
#                                            first.reminder.response.wave1.batch4) 
# 
# first.reminder.response.wave1.data <- first.reminder.response.wave1$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)


# Update take-up proportions ----------------------------------------------

takeup.prop.data.10 <- census.data %>% 
  count(wave, cluster.id, village) %>% 
  ungroup %>% 
  rename(village.pop = n) %>% 
  left_join(takeup.data %>% 
              count(village) %>% 
              rename(village.takers = n),
            "village") %>% 
  mutate(village.takers = coalesce(village.takers, as.integer(0)),
         takeup.prop = pmin(1, village.takers / village.pop),
         takeup.prop.10 = ceiling(takeup.prop * 10) %>% as.integer) %>% 
  mutate(takeup.prop.grp = case_when(is.na(.$takeup.prop) ~ "",
                                .$takeup.prop == 0 ~ "hakuna",
                                .$takeup.prop <= 0.2 ~ "wachache",
                                .$takeup.prop <= 0.4 ~ "karibu nusu",
                                .$takeup.prop <= 0.6 ~ "nusu",
                                .$takeup.prop <= 0.8 ~ "zaidi ya nusu",
                                .$takeup.prop < 1 ~ "karibu wote",
                                .$takeup.prop == 1 ~ "wote") %>% 
           na_if(""))

# takeup.prop.data.3.cluster <- census.data %>% 
#   count(cluster.id) %>% 
#   rename(cluster.pop = n) %>% 
#   left_join(takeup.data %>% 
#               count(cluster.id) %>% 
#               rename(cluster.takers = n),
#             "cluster.id") %>% 
#   mutate(cluster.takers = coalesce(cluster.takers, as.integer(0)),
#          takeup.prop = pmin(1, cluster.takers / cluster.pop),
#          takeup.prop.10 = ceiling(takeup.prop * 10) %>% as.integer) %>% 
#   mutate(takeup.prop.grp = case_when(is.na(.$takeup.prop) ~ "",
#                                 .$takeup.prop == 0 ~ "hakuna",
#                                 .$takeup.prop <= 0.2 ~ "wachache",
#                                 .$takeup.prop <= 0.4 ~ "karibu nusu",
#                                 .$takeup.prop <= 0.6 ~ "nusu",
#                                 .$takeup.prop <= 0.8 ~ "zaidi ya nusu",
#                                 .$takeup.prop < 1 ~ "karibu wote",
#                                 .$takeup.prop == 1 ~ "wote") %>% 
#            na_if(""))

# Second Message (Wave 1) -------------------------------------------------

wave1.sms.subjects <- sms.recruit.data %>% 
  filter(county != "Kakamega", cluster.id != 1) %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13,
         pot.cluster.id = ifelse(cluster.id == 277, 503, cluster.id)) %>% 
  left_join(pot.desciption.wave1, c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  left_join(select(old.census.villages, KEY.individ, village), "KEY.individ") %>% {
    wrong.key.type.data <- filter(., !is.na(KEY.individ), !str_detect(KEY.individ, "\\d\\]$"))
    
    anti_join(., wrong.key.type.data, "KEY") %>% 
      bind_rows(wrong.key.type.data %>% 
                  select(-village) %>% 
                  left_join(select(hh.census.data, KEY, village), c("KEY.individ" = "KEY")))
  } %>% 
  left_join(takeup.prop.data, "village") %>% 
  mutate(sms.treatment = ifelse(cluster.id %in% clusters.to.drop, "reminder.only", as.character(sms.treatment.recruit)))


# batch1 <- wave1.sms.subjects %>% 
#   filter(cluster.id == 277) %>% 
#   mid.intervention.messaging
#   
# batch2 <- wave1.sms.subjects %>% 
#   filter(cluster.id %in% c(491, 492)) %>% 
#   mid.intervention.messaging
#     
# batch3 <- wave1.sms.subjects %>% 
#   filter(cluster.id == 147) %>% 
#   mid.intervention.messaging
# 
# batch4 <- wave1.sms.subjects %>% 
#   filter(!cluster.id %in% c(277, 491, 492, 147)) %>% 
#   mid.intervention.messaging
# 
# second.sms.response.wave1 <- bind_rows(batch1, batch2, batch3, batch4)
# 
# second.sms.response.wave1.data <- second.sms.response.wave1$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)
# 
# save(wave1.sms.subjects, takeup.prop.data, second.sms.response.wave1, file = "sms_2_wave1.RData")

# Third Message (Wave 1) -------------------------------------------------

wave1.sms.subjects <- sms.recruit.data %>% 
  filter(county != "Kakamega", cluster.id != 1) %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13,
         pot.cluster.id = ifelse(cluster.id == 277, 503, cluster.id)) %>% 
  left_join(pot.desciption.wave1, c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  left_join(select(old.census.villages, KEY.individ, village), "KEY.individ") %>% {
    wrong.key.type.data <- filter(., !is.na(KEY.individ), !str_detect(KEY.individ, "\\d\\]$"))
    
    anti_join(., wrong.key.type.data, "KEY") %>% 
      bind_rows(wrong.key.type.data %>% 
                  select(-village) %>% 
                  left_join(select(hh.census.data, KEY, village), c("KEY.individ" = "KEY")))
  } %>% 
  left_join(takeup.prop.data.2, "village") %>% 
  mutate(sms.treatment = ifelse(cluster.id %in% clusters.to.drop, "reminder.only", as.character(sms.treatment.recruit)))

# batch.2.1 <- wave1.sms.subjects %>%
#   filter(cluster.id == 277) %>%
#   mid.intervention.messaging
# 
# batch.2.2 <- wave1.sms.subjects %>%
#   filter(cluster.id %in% c(491, 492)) %>%
#   mid.intervention.messaging
# 
# batch.2.3 <- wave1.sms.subjects %>%
#   filter(cluster.id == 147) %>%
#   mid.intervention.messaging
# 
# batch.2.4 <- wave1.sms.subjects %>%
#   filter(!cluster.id %in% c(277, 491, 492, 147)) %>%
#   mid.intervention.messaging

# sms.2.response.wave1 <- bind_rows(batch.2.1, batch.2.2, batch.2.3, batch.2.4)
# 
# sms.2.response.wave1.data <- sms.2.response.wave1$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)
# 
# save(wave1.sms.subjects, takeup.prop.data.2, sms.2.response.wave1, file = "sms_3_wave1.RData")


# Fourth Message (Wave 1) -------------------------------------------------

wave1.sms.subjects <- sms.recruit.data %>% 
  filter(county != "Kakamega", cluster.id != 1) %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13,
         pot.cluster.id = ifelse(cluster.id == 277, 503, cluster.id)) %>% 
  left_join(pot.desciption.wave1, c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  left_join(select(old.census.villages, KEY.individ, village), "KEY.individ") %>% {
    wrong.key.type.data <- filter(., !is.na(KEY.individ), !str_detect(KEY.individ, "\\d\\]$"))
    
    anti_join(., wrong.key.type.data, "KEY") %>% 
      bind_rows(wrong.key.type.data %>% 
                  select(-village) %>% 
                  left_join(select(hh.census.data, KEY, village), c("KEY.individ" = "KEY")))
  } %>% 
  left_join(takeup.prop.data.3, c("cluster.id", "village")) %>% 
  mutate(sms.treatment = ifelse(cluster.id %in% clusters.to.drop, "reminder.only", as.character(sms.treatment.recruit)))

# batch.3.1 <- wave1.sms.subjects %>%
#   filter(cluster.id == 277) %>%
#   mid.intervention.messaging
# 
# batch.3.2 <- wave1.sms.subjects %>%
#   filter(cluster.id %in% c(491, 492)) %>%
#   mid.intervention.messaging
# 
# batch.3.3 <- wave1.sms.subjects %>%
#   filter(cluster.id == 147) %>%
#   mid.intervention.messaging
# 
# batch.3.4 <- wave1.sms.subjects %>%
#   filter(!cluster.id %in% c(277, 491, 492, 147)) %>%
#   mid.intervention.messaging
# 
# sms.3.response.wave1 <- bind_rows(batch.3.1, batch.3.2, batch.3.3, batch.3.4)
# 
# sms.3.response.wave1.data <- sms.3.response.wave1$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)
# 
# save(wave1.sms.subjects, takeup.prop.data.3, sms.3.response.wave1, file = "sms_4_wave1.RData")

# Fifth Message (Wave 1) -------------------------------------------------

wave1.sms.subjects <- sms.recruit.data %>% 
  filter(county != "Kakamega", cluster.id != 1) %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13,
         pot.cluster.id = ifelse(cluster.id == 277, 503, cluster.id)) %>% 
  left_join(pot.desciption.wave1, c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  left_join(select(old.census.villages, KEY.individ, village), "KEY.individ") %>% {
    wrong.key.type.data <- filter(., !is.na(KEY.individ), !str_detect(KEY.individ, "\\d\\]$"))
    
    anti_join(., wrong.key.type.data, "KEY") %>% 
      bind_rows(wrong.key.type.data %>% 
                  select(-village) %>% 
                  left_join(select(hh.census.data, KEY, village), c("KEY.individ" = "KEY")))
  } %>% 
  left_join(takeup.prop.data.4, c("cluster.id", "village")) %>% 
  mutate(sms.treatment = ifelse(cluster.id %in% clusters.to.drop, "reminder.only", as.character(sms.treatment.recruit)))

# batch.4.1 <- wave1.sms.subjects %>%
#   filter(cluster.id == 277) %>%
#   mid.intervention.messaging
# 
# batch.4.2 <- wave1.sms.subjects %>%
#   filter(cluster.id %in% c(491, 492)) %>%
#   mid.intervention.messaging
# 
# batch.4.3 <- wave1.sms.subjects %>%
#   filter(cluster.id == 147) %>%
#   mid.intervention.messaging
# 
# batch.4.4 <- wave1.sms.subjects %>%
#   filter(!cluster.id %in% c(277, 491, 492, 147)) %>%
#   mid.intervention.messaging

# sms.4.response.wave1 <- bind_rows(batch.4.1, batch.4.2, batch.4.3, batch.4.4)
# 
# sms.4.response.wave1.data <- sms.4.response.wave1$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)
# 
# save(wave1.sms.subjects, takeup.prop.data.4, sms.4.response.wave1, file = "sms_5_wave1.RData")

# Sixth Message (Wave 1) -------------------------------------------------

wave1.sms.subjects <- sms.recruit.data %>% 
  filter(county != "Kakamega", cluster.id != 1) %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13,
         pot.cluster.id = ifelse(cluster.id == 277, 503, cluster.id)) %>% 
  left_join(pot.desciption.wave1, c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  left_join(select(old.census.villages, KEY.individ, village), "KEY.individ") %>% {
    wrong.key.type.data <- filter(., !is.na(KEY.individ), !str_detect(KEY.individ, "\\d\\]$"))
    
    anti_join(., wrong.key.type.data, "KEY") %>% 
      bind_rows(wrong.key.type.data %>% 
                  select(-village) %>% 
                  left_join(select(hh.census.data, KEY, village), c("KEY.individ" = "KEY")))
  } %>% 
  left_join(takeup.prop.data.5, c("cluster.id", "village")) %>% 
  mutate(sms.treatment = ifelse(cluster.id %in% clusters.to.drop, "reminder.only", as.character(sms.treatment.recruit)))

# batch.5.1 <- wave1.sms.subjects %>%
#   filter(cluster.id == 277) %>%
#   mid.intervention.messaging
# 
# batch.5.2 <- wave1.sms.subjects %>%
#   filter(cluster.id %in% c(491, 492)) %>%
#   mid.intervention.messaging
# 
# batch.5.3 <- wave1.sms.subjects %>%
#   filter(cluster.id == 147) %>%
#   mid.intervention.messaging
# 
# batch.5.4 <- wave1.sms.subjects %>%
#   filter(!cluster.id %in% c(277, 491, 492, 147)) %>%
#   mid.intervention.messaging
# 
# sms.5.response.wave1 <- bind_rows(batch.5.1, batch.5.2, batch.5.3, batch.5.4)
# 
# sms.5.response.wave1.data <- sms.5.response.wave1$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)
# 
# save(wave1.sms.subjects, takeup.prop.data.5, sms.5.response.wave1, file = "sms_6_wave1.RData")

# First Reminder for All (Wave 2) -----------------------------------------

wave2.sms.subjects <- sms.recruit.data %>% 
  filter(county == "Kakamega") %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13) %>% 
  left_join(pot.desciption.wave2, "cluster.id") %>% # c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  group_by(cluster.id) %>% 
  mutate(cluster.msg = sprintf(reminder.message, description.kiswahili)) %>% 
  ungroup

# first.reminder.response.wave2.batch1 <- wave2.sms.subjects %>%
#   filter(cluster.id == 517) %>%
#   cluster.reminder.messenger

# first.reminder.response.wave2.batch2 <- wave2.sms.subjects %>%
#   filter(cluster.id %in% c(551, 555)) %>%
#   cluster.reminder.messenger

# first.reminder.response.wave2.batch3 <- wave2.sms.subjects %>%
#   filter(!cluster.id %in% c(517, 551, 555)) %>%
#   cluster.reminder.messenger

# first.reminder.response.wave1 <- bind_rows(first.reminder.response.wave1.batch1,
#                                            first.reminder.response.wave1.batch2,
#                                            first.reminder.response.wave1.batch3,
#                                            first.reminder.response.wave1.batch4)
# 
# first.reminder.response.wave1.data <- first.reminder.response.wave1$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)

# Second Message (Wave 2) -------------------------------------------------

wave2.sms.subjects <- sms.recruit.data %>% 
  filter(county == "Kakamega") %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13) %>% 
  left_join(select(census.data, KEY.individ, airtime.sample, airtime.order), "KEY.individ") %>% 
  group_by(KEY.individ, airtime.sample) %>% 
  mutate(num.phone.for.same.person = n()) %>% { # Recruited multiple times and with different numbers, take first one 
    if (first(.$airtime.sample)) {
      return(.) # Keep them all!
    } else {
      filter(., min_rank(starttime) == 1) 
    }
  } %>% 
  group_by(phone) %>% # We have a number of people with the same phone number. We'll only message one time
  mutate(num.phone.repeats = n()) %>% 
  slice(1) %>% 
  ungroup %>% 
  left_join(pot.desciption.wave2, "cluster.id") %>% # c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  left_join(select(old.census.villages, KEY.individ, village), "KEY.individ") %>% {
    wrong.key.type.data <- filter(., !is.na(KEY.individ), !str_detect(KEY.individ, "\\d\\]$"))
    
    anti_join(., wrong.key.type.data, "KEY") %>% 
      bind_rows(wrong.key.type.data %>% 
                  select(-village) %>% 
                  left_join(select(hh.census.data, KEY, village), c("KEY.individ" = "KEY")))
  } %>% 
  left_join(takeup.prop.data.6, c("cluster.id", "village")) %>% 
  mutate(sms.treatment = as.character(sms.treatment.recruit)) %>% 
  filter(!is.na(takeup.prop.10))

# wave2.batch1.1 <- wave2.sms.subjects %>%
#   filter(cluster.id == 517) %>%
#   mid.intervention.messaging

# wave2.batch1.2 <- wave2.sms.subjects %>%
#   filter(cluster.id %in% c(551, 555)) %>%
#   mid.intervention.messaging

# wave2.batch1.3 <- wave2.sms.subjects %>%
#   filter(!cluster.id %in% c(517, 551, 555)) %>%
#   mid.intervention.messaging

# wave2.sms.subjects %>%
#   filter(airtime.sample, sms.treatment == "reminder.only" | airtime.order == 1) %>% 
#   reward.pilot.msg.sender

# wave2.sms.subjects %>% 
#   filter(airtime.sample, sms.treatment != "reminder.only" & airtime.order == 2) %>% dim

# sms.1.response.wave2 <- bind_rows(wave2.batch1.1, wave2.batch1.2, wave2.batch1.3)
# 
# sms.1.response.wave2.data <- sms.1.response.wave2$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)
# 
# save(wave2.sms.subjects, takeup.prop.data.6, sms.1.response.wave2, file = "sms_1_wave2.RData")

# Third Message (Wave 2) -------------------------------------------------

wave2.sms.subjects <- sms.recruit.data %>% 
  filter(county == "Kakamega") %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13) %>% 
  left_join(select(census.data, KEY.individ, airtime.sample, airtime.order), "KEY.individ") %>% 
  group_by(KEY.individ, airtime.sample) %>% 
  mutate(num.phone.for.same.person = n()) %>% { # Recruited multiple times and with different numbers, take first one 
    if (first(.$airtime.sample)) {
      return(.) # Keep them all!
    } else {
      filter(., min_rank(starttime) == 1) 
    }
  } %>% 
  group_by(phone) %>% # We have a number of people with the same phone number. We'll only message one time
  mutate(num.phone.repeats = n()) %>% 
  slice(1) %>% 
  ungroup %>% 
  left_join(pot.desciption.wave2, "cluster.id") %>% # c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  left_join(select(old.census.villages, KEY.individ, village), "KEY.individ") %>% {
    wrong.key.type.data <- filter(., !is.na(KEY.individ), !str_detect(KEY.individ, "\\d\\]$"))
    
    anti_join(., wrong.key.type.data, "KEY") %>% 
      bind_rows(wrong.key.type.data %>% 
                  select(-village) %>% 
                  left_join(select(hh.census.data, KEY, village), c("KEY.individ" = "KEY")))
  } %>% 
  left_join(takeup.prop.data.7, c("cluster.id", "village")) %>% 
  mutate(sms.treatment = as.character(sms.treatment.recruit)) %>% 
  filter(!is.na(takeup.prop.10))

# wave2.batch2.1 <- wave2.sms.subjects %>%
#   filter(cluster.id == 517) %>%
#   mid.intervention.messaging

# wave2.batch2.2 <- wave2.sms.subjects %>%
#   filter(cluster.id %in% c(551, 555)) %>%
#   mid.intervention.messaging

# wave2.batch2.3 <- wave2.sms.subjects %>%
#   filter(!cluster.id %in% c(517, 551, 555)) %>%
#   mid.intervention.messaging

# wave2.sms.subjects %>%
#   filter(airtime.sample, sms.treatment == "reminder.only" | airtime.order == 1) %>% 
#   reward.pilot.msg.sender

# sms.1.response.wave2 <- bind_rows(wave2.batch1.1, wave2.batch1.2, wave2.batch1.3)
# 
# sms.1.response.wave2.data <- sms.1.response.wave2$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)
# 
# save(wave2.sms.subjects, takeup.prop.data.6, sms.1.response.wave2, file = "sms_1_wave2.RData")

# Fourth Message (Wave 2) -------------------------------------------------

wave2.sms.subjects <- sms.recruit.data %>% 
  filter(county == "Kakamega") %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13) %>% 
  left_join(select(census.data, KEY.individ, airtime.sample, airtime.order), "KEY.individ") %>% 
  group_by(KEY.individ, airtime.sample) %>% 
  mutate(num.phone.for.same.person = n()) %>% { # Recruited multiple times and with different numbers, take first one 
    if (first(.$airtime.sample)) {
      return(.) # Keep them all!
    } else {
      filter(., min_rank(starttime) == 1) 
    }
  } %>% 
  group_by(phone) %>% # We have a number of people with the same phone number. We'll only message one time
  mutate(num.phone.repeats = n()) %>% 
  slice(1) %>% 
  ungroup %>% 
  left_join(pot.desciption.wave2, "cluster.id") %>% # c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  left_join(select(old.census.villages, KEY.individ, village), "KEY.individ") %>% {
    wrong.key.type.data <- filter(., !is.na(KEY.individ), !str_detect(KEY.individ, "\\d\\]$"))
    
    anti_join(., wrong.key.type.data, "KEY") %>% 
      bind_rows(wrong.key.type.data %>% 
                  select(-village) %>% 
                  left_join(select(hh.census.data, KEY, village), c("KEY.individ" = "KEY")))
  } %>% 
  left_join(takeup.prop.data.8, c("cluster.id", "village")) %>% 
  mutate(sms.treatment = as.character(sms.treatment.recruit)) %>% 
  filter(!is.na(takeup.prop.10))

# wave2.batch3.1 <- wave2.sms.subjects %>%
#   filter(cluster.id == 517) %>%
#   mid.intervention.messaging

# wave2.batch3.2 <- wave2.sms.subjects %>%
#   filter(cluster.id %in% c(551, 555)) %>%
#   mid.intervention.messaging

# wave2.batch3.3 <- wave2.sms.subjects %>%
#   filter(!cluster.id %in% c(517, 551, 555)) %>%
#   mid.intervention.messaging

# wave2.sms.subjects %>%
#   filter(airtime.sample, sms.treatment != "reminder.only" & airtime.order == 2) %>%
#   reward.pilot.msg.sender

# sms.1.response.wave2 <- bind_rows(wave2.batch1.1, wave2.batch1.2, wave2.batch1.3)
# 
# sms.1.response.wave2.data <- sms.1.response.wave2$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)
# 
# save(wave2.sms.subjects, takeup.prop.data.6, sms.1.response.wave2, file = "sms_1_wave2.RData")

# Fifth Message (Wave 2) -------------------------------------------------

wave2.sms.subjects <- sms.recruit.data %>% 
  filter(county == "Kakamega") %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13) %>% 
  left_join(select(census.data, KEY.individ, airtime.sample, airtime.order), "KEY.individ") %>% 
  group_by(KEY.individ, airtime.sample) %>% 
  mutate(num.phone.for.same.person = n()) %>% { # Recruited multiple times and with different numbers, take first one 
    if (first(.$airtime.sample)) {
      return(.) # Keep them all!
    } else {
      filter(., min_rank(starttime) == 1) 
    }
  } %>% 
  group_by(phone) %>% # We have a number of people with the same phone number. We'll only message one time
  mutate(num.phone.repeats = n()) %>% 
  slice(1) %>% 
  ungroup %>% 
  left_join(pot.desciption.wave2, "cluster.id") %>% # c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  left_join(select(old.census.villages, KEY.individ, village), "KEY.individ") %>% {
    wrong.key.type.data <- filter(., !is.na(KEY.individ), !str_detect(KEY.individ, "\\d\\]$"))
    
    anti_join(., wrong.key.type.data, "KEY") %>% 
      bind_rows(wrong.key.type.data %>% 
                  select(-village) %>% 
                  left_join(select(hh.census.data, KEY, village), c("KEY.individ" = "KEY")))
  } %>% 
  left_join(takeup.prop.data.9, c("cluster.id", "village")) %>% 
  mutate(sms.treatment = as.character(sms.treatment.recruit)) %>% 
  filter(!is.na(takeup.prop.10))

# wave2.batch4.1 <- wave2.sms.subjects %>%
#   filter(cluster.id == 517) %>%
#   mid.intervention.messaging

# wave2.batch4.2 <- wave2.sms.subjects %>%
#   filter(cluster.id %in% c(551, 555)) %>%
#   mid.intervention.messaging

# wave2.batch4.3 <- wave2.sms.subjects %>%
#   filter(!cluster.id %in% c(517, 551, 555)) %>%
#   mid.intervention.messaging

# sms.1.response.wave2 <- bind_rows(wave2.batch1.1, wave2.batch1.2, wave2.batch1.3)
# 
# sms.1.response.wave2.data <- sms.1.response.wave2$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)
# 
# save(wave2.sms.subjects, takeup.prop.data.6, sms.1.response.wave2, file = "sms_1_wave2.RData")

# Sixth Message (Wave 2) -------------------------------------------------

wave2.sms.subjects <- sms.recruit.data %>% 
  filter(county == "Kakamega") %>% 
  mutate(phone = ifelse(nchar(phone) == 10 & str_sub(phone, 1, 1) == "0", str_sub(phone, 2), phone),
         phone = paste0("+254", phone),
         valid.phone = nchar(phone) == 13) %>% 
  left_join(select(census.data, KEY.individ, airtime.sample, airtime.order), "KEY.individ") %>% 
  group_by(KEY.individ, airtime.sample) %>% 
  mutate(num.phone.for.same.person = n()) %>% { # Recruited multiple times and with different numbers, take first one 
    if (first(.$airtime.sample)) {
      return(.) # Keep them all!
    } else {
      filter(., min_rank(starttime) == 1) 
    }
  } %>% 
  group_by(phone) %>% # We have a number of people with the same phone number. We'll only message one time
  mutate(num.phone.repeats = n()) %>% 
  slice(1) %>% 
  ungroup %>% 
  left_join(pot.desciption.wave2, "cluster.id") %>% # c("pot.cluster.id" = "cluster.id")) %>% 
  filter(valid.phone) %>% 
  left_join(select(old.census.villages, KEY.individ, village), "KEY.individ") %>% {
    wrong.key.type.data <- filter(., !is.na(KEY.individ), !str_detect(KEY.individ, "\\d\\]$"))
    
    anti_join(., wrong.key.type.data, "KEY") %>% 
      bind_rows(wrong.key.type.data %>% 
                  select(-village) %>% 
                  left_join(select(hh.census.data, KEY, village), c("KEY.individ" = "KEY")))
  } %>% 
  left_join(takeup.prop.data.10, c("cluster.id", "village")) %>% 
  mutate(sms.treatment = as.character(sms.treatment.recruit)) %>% 
  filter(!is.na(takeup.prop.10))

# wave2.batch5.1 <- wave2.sms.subjects %>%
#   filter(cluster.id == 517) %>%
#   mid.intervention.messaging

# wave2.batch5.2 <- wave2.sms.subjects %>%
#   filter(cluster.id %in% c(551, 555)) %>%
#   mid.intervention.messaging

wave2.batch5.3 <- wave2.sms.subjects %>%
  filter(!cluster.id %in% c(517, 551, 555)) %>%
  group_by(cluster.id) %>% 
  do(response = mid.intervention.messaging(.))

# sms.1.response.wave2 <- bind_rows(wave2.batch1.1, wave2.batch1.2, wave2.batch1.3)
# 
# sms.1.response.wave2.data <- sms.1.response.wave2$response %>% map_df(~ content(.)$SMSMessageData$Recipients %>% bind_rows)
# 
# save(wave2.sms.subjects, takeup.prop.data.6, sms.1.response.wave2, file = "sms_1_wave2.RData")

# Find Airtime Reward Recipients --------------------------------------------------

airtime.outbox <- read_csv("takeup_airtime_outbox.csv",
                           skip = 1,
                           col_names = c("airtime.datetime", "transact.id", "recipient", "amount", "discount", "status"),
                           col_types = list(airtime.datetime = col_datetime("%I:%M %p %B %d, %Y"))) %>% 
  mutate(recipient = paste0("+", recipient))

all.airtime.inbox <- read_csv("takeup_airtime_inbox.csv", 
         skip = 1,
         col_names = c("message.datetime", "from", "to", "text"),
         col_types = list(message.datetime = col_datetime("%I:%M %p %B %d, %Y"),
                          from = col_character())) %>% 
  filter(date(message.datetime) >= "2016-10-25") %>% 
  distinct(from, .keep_all = TRUE) %>% 
  left_join(wave2.sms.subjects, c("from" = "phone")) %>% 
  filter(airtime.sample) 

airtime.inbox <- all.airtime.inbox %>% 
  anti_join(airtime.outbox, c("from" = "recipient")) %>% 
  arrange(desc(message.datetime)) %>% 
  select(message.datetime, from, KEY.individ, airtime.sample, cluster.id, sms.treatment) 

# Reward Reminder ---------------------------------------------------------

# wave2.sms.subjects %>%
#   filter(airtime.sample, airtime.order == 2 & sms.treatment != "reminder.only") %>%
#   anti_join(all.airtime.inbox, c("phone" = "from")) %>%
#   reward.reminder.msg.sender

# Send Airtime ------------------------------------------------------------

pilot.airtime.response.1 <- POST("http://api.africastalking.com/version1/airtime/send",
                 accept_json(),
                 add_headers(Apikey = config$africastalking_api_key),
                 body = list(username = "karimn",
                             recipients = map(airtime.inbox$from, ~ list(phoneNumber = ., amount = "KES 50")) %>% toJSON(auto_unbox = TRUE),
                             from = "EvidenceAct"),
                 encode = "form",
                 verbose())

