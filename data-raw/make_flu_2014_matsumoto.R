library(tidyverse)

url <- "https://raw.githubusercontent.com/akira-endo/HHstudy_FluMatsumoto2014-15/master/HHdata_FluMatsumoto2014-15.csv"

read_csv(url) |>
  uncount(Frequency) |>
  mutate(house_id = seq_len(n())) |>
  pivot_longer(cols = matches("\\.(num|inf)$"),
               names_to = c("member", ".value"),
               names_sep = "\\.",
               values_drop_na = FALSE) |>
  filter(num > 0) |>
  filter(num > inf) |>
  mutate(status = lapply(seq_len(n()), function(i) {
    c(rep(1L, inf[i]), rep(0L, num[i] - inf[i]))
  })) |>
  unnest(status) |>
  select(house_id, member, status) |>
  as.data.frame()
