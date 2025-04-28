
hand_copied <-"inst/--developer/datasets/flu_1918_baltimora.csv"

flu_1918_baltimora <- read_csv(hand_copied) %>%
  pivot_longer(
    cols = starts_with("cases_"),
    names_to = "case_type",
    values_to = "n",
    values_drop_na = TRUE
  ) %>%
  mutate(inf = readr::parse_number(case_type)) %>%
  select(-case_type) %>%
  uncount(n) %>%
  mutate(house_id = seq_len(n())) %>%
  mutate(
    status = lapply(seq_len(n()), function(i) {
      c(rep(1L, inf[i]), rep(0L, size[i] - inf[i]))
    })
  ) %>%
  unnest(status) %>%
  select(-size, -inf)

usethis::use_data(flu_1918_baltimora, overwrite = TRUE)
