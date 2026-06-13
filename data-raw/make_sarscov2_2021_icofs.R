# path to public pre-processed data
url <- "https://gitlab.pasteur.fr/mmmi-pasteur/hh_israel_sarscov2_delta2/-/raw/main/codes/db_filtered.txt"

# format dates as number of days from beginning of 2020
date_fmt <- function(x) as.numeric(as.Date(x) - as.Date("2020-01-01"))

sarscov2_2021_icofs <- utils::read.table(url) |>
  dplyr::group_by(hhid) |>
  dplyr::mutate(
    house_id = sprintf("hh%03d", dplyr::cur_group_id()),
    id = sprintf("%sm%02d", house_id, dplyr::row_number()),
    # follow-up starts at the household index's first positive PCR
    start_followup = date_fmt(dplyr::first(index_first_pos_pcr)),
    # 21-day follow-up from the household's index case (paper's studyPeriod)
    end_followup = start_followup + 21,
    # using original vaccination category (7-day effectiveness rule)
    vacc_status = dplyr::case_when(
      vaccination_cat == "unvaccinated" ~ "unvaccinated",
      vaccination_cat == "1dose" ~ "1dose",
      vaccination_cat == "2doses" & vaccine_delay <= 90 ~ "2doses_leq90",
      vaccination_cat == "2doses" & vaccine_delay > 90 ~ "2doses_gtr90",
      vaccination_cat == "boosted" ~ "3doses"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    status = as.integer(!is.na(first_pos)),
    test_date = date_fmt(first_pos),

    isolation = isolation != "no quarantine",
    recovered = as.logical(recovered),
    # raw `asymptomatic == 1` actually means symptomatic
    # compare to ../raw/main/codes/data_filtering.R:188
    asymptomatic = !as.logical(asymptomatic)
  ) |>
  # zero-out positives whose test falls after the 21-day window
  dplyr::mutate(
    status = ifelse(status > 0 & test_date > end_followup, 0L, status),
    test_date = ifelse(status == 0, NA, test_date)
  ) |>
  dplyr::group_by(house_id) |>
  dplyr::filter(sum(status > 0) > 0) |>
  dplyr::ungroup() |>
  dplyr::select(
    id, house_id, status, start_followup, end_followup, test_date, age, sex,
    vacc_status, recovered, isolation, asymptomatic
  ) |>
  as.data.frame()

usethis::use_data(sarscov2_2021_icofs, overwrite = TRUE)
