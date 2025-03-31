
url <- "https://gitlab.pasteur.fr/mmmi-pasteur/hh_israel_sarscov2_delta2/-/raw/main/codes/db_filtered.txt"
dataset <- read.table(url)

date_fmt <- function(x) {
  as.numeric(as.Date(x, format = "%Y-%m-%d")  - as.Date("2021-01-01"))
}

study_period <- 21

sarscov2_2021_israel <- dataset %>%
  mutate(
    isolations = ifelse(isolation == "no quarantine", 0, 1),
    status = ifelse(is.na(first_pos), 0, 1),
    dds = ifelse(is.na(first_pos), NA, date_fmt(first_pos)),
    severity = asymptomatic
  ) %>%
  ungroup() %>%
  select(indid, hhid, dds, status, vaccination_cat, age, isolations, sex, severity)


head(sarscov2_2021_israel)

# usethis::use_data(sarscov2_2021_israel)
