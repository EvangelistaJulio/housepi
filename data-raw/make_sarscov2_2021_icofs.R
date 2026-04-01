
# path to public pre-processed data
url <- "https://gitlab.pasteur.fr/mmmi-pasteur/hh_israel_sarscov2_delta2/-/raw/main/codes/db_filtered.txt"

date_fmt <- function(x) as.numeric(as.Date(x) - as.Date("2020-01-01"))

sarscov2_2021_icofs <- utils::read.table(url) |> 
  dplyr::group_by(hhid) |>
  # define followup period
  dplyr::mutate(
    house_id = sprintf("hh%03d", dplyr::cur_group_id()),
    id = sprintf("%sm%02d", house_id, dplyr::row_number()),
    end_followup = if(!all(is.na(last_neg))) {
      max(date_fmt(last_neg), na.rm = TRUE)
    } else {
      date_fmt(dplyr::first(index_first_pos_pcr)) + 21
    }
  ) |>
  dplyr::ungroup() |> 
  # format variables
  dplyr::mutate(
    status = as.integer(!is.na(first_pos)),
    test_date = date_fmt(first_pos),
    last_vacc = pmax(vaccine_1, vaccine_2, vaccine_3, na.rm = TRUE),
    vacc_status = dplyr::case_when(
      is.na(last_vacc) ~ "unvaccinated",
      last_vacc == vaccine_3 ~ "3doses",
      last_vacc == vaccine_2 ~ "2doses",
      last_vacc == vaccine_1 ~ "1dose"
    ),
    # whether last vaccine is more than 3 months old
    vacc_status = dplyr::case_when(
      vacc_status == "2doses" ~ ifelse(
        date_fmt(index_first_pos_pcr) - date_fmt(last_vacc) <= 90,
        "2doses_leq90", "2doses_gtr90"
      ),
      .default = vacc_status
    ),
    isolation = isolation != "no quarantine",
    asymptomatic = as.logical(asymptomatic)
  ) |>
  # drop positive status if test after end followup
  dplyr::mutate(
    status = ifelse(status > 0 & test_date > end_followup, 0, status),
    test_date = ifelse(status == 0, NA, test_date),
  ) |>
  # drop resulting households with no cases
  dplyr::group_by(house_id) |>
  dplyr::filter(sum(status > 0) > 0) |>
  dplyr::ungroup() |>
  dplyr::select(id, house_id, status, end_followup, test_date, age, sex,
                vacc_status, isolation, asymptomatic) |>
  as.data.frame()

usethis::use_data(sarscov2_2021_icofs, overwrite = TRUE)
