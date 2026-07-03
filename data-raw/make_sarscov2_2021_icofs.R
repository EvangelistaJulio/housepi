# path to public pre-processed data
url <- "https://gitlab.pasteur.fr/mmmi-pasteur/hh_israel_sarscov2_delta2/-/raw/main/codes/db_filtered.txt"

sarscov2_2021_icofs <- utils::read.table(url) |>
  # rename household and respective members IDs
  dplyr::group_by(hhid) |>
  dplyr::mutate(
    house_id = sprintf("hh%03d", dplyr::cur_group_id()),
    id = sprintf("%sm%02d", house_id, dplyr::row_number()),
    start_date = as.Date(min(index_first_pos_pcr, na.rm = TRUE))
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    age = dplyr::case_when(age > 12 ~ "adult", .default = "child"),
    # raw `asymptomatic == 1` actually means symptomatic
    # compare to ../raw/main/codes/data_filtering.R:188
    asymptomatic = !as.logical(asymptomatic),
    status = dplyr::case_when(
      is.na(first_pos) ~ 0,
      asymptomatic ~ 2,
      .default = 1
    ),
    # 3-weeks follow-up from the household's index case
    study_period = 21L,
    end_date = start_date + study_period,
    test_date = as.Date(first_pos),
    # drop late cases
    asymptomatic = dplyr::if_else(
      status != 0 & test_date > end_date, TRUE, asymptomatic
    ),
    test_date = dplyr::if_else(
      status != 0 & test_date > end_date, NA, test_date
    ),
    status = dplyr::if_else(is.na(test_date), 0L, status),
    # Cortier et al. '25 vaccination categories (7-day effectiveness rule)
    vacc_status = dplyr::case_when(
      vaccination_cat == "unvaccinated" ~ "unvaccinated",
      vaccination_cat == "1dose" ~ "1dose",
      vaccination_cat == "2doses" & vaccine_delay <= 90 ~ "2doses_leq90",
      vaccination_cat == "2doses" & vaccine_delay > 90 ~ "2doses_gtr90",
      vaccination_cat == "boosted" ~ "3doses"
    ),
    isolated = isolation != "no quarantine",
    recovered = as.logical(recovered)
  ) |>
  dplyr::select(
    id, house_id, status, study_period, test_date, age, sex,
    vacc_status, recovered, isolated, asymptomatic
  ) |>
  as.data.frame()

usethis::use_data(sarscov2_2021_icofs, overwrite = TRUE)
