
format_flu_1918 <- function(location) {
  # path to hand-copied file
  hand_copied <- sprintf("data-raw/flu_1918_%s.csv", location)

  read.csv(hand_copied) |>
    tidyr::pivot_longer(cols = tidyr::starts_with("cases_"),
                        names_to = "case_type",
                        values_to = "n",
                        values_drop_na = TRUE) |>
    dplyr::mutate(inf = readr::parse_number(case_type)) |>
    dplyr::select(-case_type) |>
    tidyr::uncount(n) |>
    dplyr::mutate(
      house_id = seq_len(dplyr::n()),
      status = lapply(seq_len(dplyr::n()), function(i) {
      c(rep(1L, inf[i]), rep(0L, size[i] - inf[i]))
    })) |>
    tidyr::unnest(status) |>
    dplyr::select(-size, -inf) |>
    as.data.frame()
}

flu_1918_frederick <- format_flu_1918("frederick")
usethis::use_data(flu_1918_frederick, overwrite = TRUE)

flu_1918_baltimora <- format_flu_1918("baltimora")
usethis::use_data(flu_1918_baltimora, overwrite = TRUE)
