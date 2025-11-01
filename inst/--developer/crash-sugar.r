library(housepi)
library(testthat)

# Known infinity (values below minimum, i. e. proba 0)
hhmodel <- new_housepi_model("augmented_example") |> 
  add_data_wide(housepi::survey) |> 
  load_augmented_list(
    value = list(model = "uniform"),
    proba = list(func = "lognormal", 
                 func_args = c(x = ~ value, meanlog = 1, sdlog = 1), 
                 min = 1, max = 15, discr = TRUE)
  ) |> 
  evaluate_log_lik() |> 
  expect_equal(-Inf)

# Known NaN: qlnorm(15, meanlog = 4, sdlog = 0.01) ~= 0, impossible truncation
new_housepi_model("augmented_example") |> 
  add_data_wide(housepi::survey) |> 
  load_augmented_list(
    value = list(model = "uniform", args = c(min = 1, max = 10), discr = TRUE),
    proba = list(func = "lognormal", 
                 func_args = c(x = ~ value, meanlog = 4, sdlog = 0.01), 
                 min = 1, max = 15, discr = TRUE)
    ) |> 
  evaluate_log_lik() |> 
  expect_error(regexp = "Observation log-likelihood returned NaN")
