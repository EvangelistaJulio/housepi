
# path to public pre-processed data
url <- "https://raw.githubusercontent.com/akira-endo/HHstudy_FluMatsumoto2014-15/master/HHdata_FluMatsumoto2014-15.csv"

flu_2014_matsumoto <- read.csv(url) |>
  tidyr::uncount(Frequency) |>
  dplyr::mutate(house_id = seq_len(dplyr::n())) |>
  tidyr::pivot_longer(cols = matches("\\.(num|inf)$"),
                      names_to = c("member", ".value"),
                      names_sep = "\\.",
                      values_drop_na = FALSE) |>
  dplyr::filter(num > 0, num > inf) |>
  dplyr::mutate(
    status = lapply(seq_len(dplyr::n()), function(i) {
      c(rep(1L, inf[i]), rep(0L, num[i] - inf[i]))
    }),
    member = stringr::str_to_lower(member)
  ) |>
  tidyr::unnest(status) |>
  dplyr::select(house_id, status, member) |>
  as.data.frame()

usethis::use_data(flu_2014_matsumoto, overwrite = TRUE)
