
url <- "https://gitlab.pasteur.fr/mmmi-pasteur/hh_israel_sarscov2_delta2/-/raw/main/codes/db_filtered.txt"
dataset <- read.table(url)

date_fmt <- function(x) {
  as.numeric(as.Date(x, format = "%Y-%m-%d") - as.Date("2021-01-01"))
}

study_period <- 21

sarscov2_2021_israel <- dataset %>%
  mutate(
    isolation = as.integer(isolation != "no quarantine"),
    status = as.integer(!is.na(first_pos)),
    onset = date_fmt(first_pos)
  ) %>%
  ungroup() %>%
  select(indid, hhid, status, onset, vaccination_cat,
         age, isolation, sex, asymptomatic)

head(sarscov2_2021_israel)

# usethis::use_data(sarscov2_2021_israel)
