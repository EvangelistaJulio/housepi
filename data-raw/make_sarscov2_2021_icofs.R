# path to public pre-processed data
url <- "https://gitlab.pasteur.fr/mmmi-pasteur/hh_israel_sarscov2_delta2/-/raw/main/codes/db_filtered.txt"

sarscov2_2021_icofs <- utils::read.table(url) |>
  # rename household and respective members IDs
  dplyr::group_by(hhid) |>
  dplyr::mutate(
    house_id = sprintf("hh%03d", dplyr::cur_group_id()),
    id = sprintf("%sm%02d", house_id, dplyr::row_number())
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    # raw `asymptomatic == 1` actually means symptomatic
    # compare to ../raw/main/codes/data_filtering.R:188
    asymptomatic = !as.logical(asymptomatic),
    status = dplyr::case_when(
      is.na(first_pos)  ~ 0L,
      asymptomatic ~ 1L,
      .default = 2L
    ),
    # 3-weeks follow-up from the household's index case
    study_period = 21L,
    start_date = as.Date(index_first_pos_pcr),
    end_date = start_date + study_period,
    test_date = as.Date(first_pos),
    # Cortier et al. '25 vaccination categories (7-day effectiveness rule)
    vacc_status = dplyr::case_when(
      vaccination_cat == "unvaccinated" ~ "unvaccinated",
      vaccination_cat == "1dose" ~ "1dose",
      vaccination_cat == "2doses" & vaccine_delay <= 90 ~ "2doses_leq90",
      vaccination_cat == "2doses" & vaccine_delay > 90 ~ "2doses_gtr90",
      vaccination_cat == "boosted" ~ "3doses"
    ),
    isolation = isolation != "no quarantine",
    recovered = as.logical(recovered)
  ) |>
  dplyr::select(
    id, house_id, status, start_date, study_period, test_date, age, sex,
    vacc_status, recovered, isolation, asymptomatic
  ) |>
  as.data.frame()

usethis::use_data(sarscov2_2021_icofs, overwrite = TRUE)
