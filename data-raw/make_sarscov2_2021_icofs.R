# path to public pre-processed data
url <- paste0(
  "https://gitlab.pasteur.fr/mmmi-pasteur/hh_israel_sarscov2_delta2/",
  "-/raw/main/codes/db_filtered.txt"
)

# starting period of Delta wave cluster of interest in ICoFS
icofs_origin <- as.Date("2021-07-01")

# follow-up from the first positive test, in days
follow_up <- 21

sarscov2_2021_icofs <- utils::read.table(url) |>
  # NOTE: dates anchored on earliest positive PCR leaked prior-wave infections
  dplyr::mutate(
    age = dplyr::case_when(age > 12 ~ "adult", .default = "child"),
    # raw `asymptomatic == 1` actually means symptomatic
    # compare to ../raw/main/codes/data_filtering.R:188
    asymptomatic = !as.logical(asymptomatic),
    first_pos = as.Date(first_pos),
    # pre-wave positive is prior infection, over the questionnaire
    recovered = as.logical(recovered) |
      (!is.na(first_pos) & first_pos < icofs_origin),
    # in-window positive is the observed case date; NA otherwise
    test_date = dplyr::if_else(
      !is.na(first_pos) & first_pos >= icofs_origin, first_pos, as.Date(NA)
    )
  ) |>
  # rename household and respective members IDs
  dplyr::group_by(raw_hh_id = hhid) |>
  dplyr::mutate(
    hh_id = sprintf("hh%03d", dplyr::cur_group_id()),
    id = sprintf("%sm%02d", hh_id, dplyr::row_number()),
    # first positive test
    start_date = min(test_date, na.rm = TRUE),
    # fixed follow-up, as in the original analysis
    end_date = start_date + follow_up
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    # positive after the follow-up are assumed escaped
    test_date = dplyr::if_else(test_date > end_date, as.Date(NA), test_date),
    # escaped: 0, infected symptomatic: 1, asymptomatic: 2
    status = dplyr::case_when(
      is.na(test_date) ~ 0L, asymptomatic ~ 2L, .default = 1L
    ),
    # observed index case: positive on the household's first test day
    index_case = as.integer(!is.na(test_date) & test_date == start_date),
    # re-derive vaccination categories with fixed icofs_origin
    # 7-day activation, >90-day split on dose
    delay_1 = as.numeric(start_date - as.Date(vaccine_1)),
    delay_2 = as.numeric(start_date - as.Date(vaccine_2)),
    delay_3 = as.numeric(start_date - as.Date(vaccine_3)),
    n_doses = (!is.na(delay_1) & delay_1 > 7) +
      (!is.na(delay_2) & delay_2 > 7) +
      (!is.na(delay_3) & delay_3 > 7),
    vacc_status = dplyr::case_when(
      n_doses == 0 ~ "unvaccinated",
      n_doses == 1 ~ "1dose",
      n_doses == 2 & delay_2 > 90 ~ "2doses_gtr90",
      n_doses == 2 ~ "2doses_leq90",
      n_doses >= 3 ~ "3doses"
    ),
    isolated = isolation != "no quarantine"
  ) |>
  # case-ascertained model assumes a single introduction per household
  dplyr::filter(sum(index_case) == 1, .by = hh_id) |>
  dplyr::select(
    id, hh_id, status, age, sex,
    start_date, end_date, test_date,
    vacc_status, recovered, isolated, asymptomatic, index_case
  ) |>
  as.data.frame()

usethis::use_data(sarscov2_2021_icofs, overwrite = TRUE)
