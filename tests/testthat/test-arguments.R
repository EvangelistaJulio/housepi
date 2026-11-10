library(testthat)

####  Definitions  ####

test_that("Settings are formatted as character", {
  # id_draws with multiple values
  result <- char_settings_input_(list(steps = 100, adapt = 10, thin = 1,
                                      chains = 1, id_draws = c("m1", "m2")))
  result |> expect_s3_class("data.frame")
  result$id_draws |> expect_identical("m1:m2")
  vapply(result, is.character, logical(1)) |> expect_all_true()
  # id_draws are dropped if NULL
  char_settings_input_(list(steps = 100L, id_draws = NULL)) |>
    getElement("id_draws") |>
    expect_identical("")
})

test_that("A count is a whole positive number", {
  assert_count_(3, "Number of replicates") |> expect_identical(3L)
  assert_count_(NA, "Number of replicates", allow_na = TRUE) |>
    expect_identical(NA_integer_)
  assert_count_(NA, "Number of replicates") |>
    expect_error("[Nn]umber of replicates should be a positive integer")
  assert_count_(2.5, "Number of replicates") |> expect_error("positive integer")
  assert_count_(c(1, 2), "Number of replicates", allow_na = TRUE) |>
    expect_error("positive integer")
})
