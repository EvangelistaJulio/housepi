library(tidyverse)

url <- "https://gitlab.pasteur.fr/mmmi-pasteur/hh_israel_sarscov2_delta2/-/raw/main/codes/db_filtered.txt"

date_fmt <- function(x) {
  as.numeric(as.Date(x, format = "%Y-%m-%d") - as.Date("2021-01-01"))
}

#study_period <- 21

sarscov2_2021_israel <- read.table(url) %>%
  mutate(
    isolation = as.integer(isolation != "no quarantine"),
    status = as.integer(!is.na(first_pos)),
    test_date = date_fmt(first_pos),
    isolation = as.logical(isolation),
    asymptomatic = as.logical(asymptomatic)
  ) %>%
  ungroup() %>%
  select(indid, hhid, status, test_date, vaccination_cat,
         age, isolation, sex, asymptomatic) %>%
  rename(all_of(c(id = "indid", house_id = "hhid", vacc_status = "vaccination_cat")))

head(sarscov2_2021_israel)

usethis::use_data(sarscov2_2021_israel, overwrite = TRUE)
