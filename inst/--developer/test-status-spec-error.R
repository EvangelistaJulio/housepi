# source loads
source("R/utils.R")
source("R/template.R")
source("R/augmented.R")
source("R/parameter.R")
source("R/infectivity.R")
source("R/group_rates.R")
source("R/expected_cases.R")
source("R/model_generator.R")
source("R/model_class.R")
source("R/sampler_class.R")
source("R/expose.R")

library(testthat)

test_that("Status-specific definitions are safe", {
  hhmodel <- new_housepi_model("augm-example") |>
    add_dataset(data.frame(
      "id" = c(1, 2, 3, 4),
      "status" = c(4, 2, 3, 1),
      "param" = c(1, 2, 3, 4)
    ))
  
  # first generation of definitions
  err_msg <- as.character(expect_error(
    hhmodel |>
      load_augmented_list(
        ..proba = list(model = "norm", status = 1, min = 0,
                       args = c(mean = "param@D")),
        ..proba = list(expr = "..layer", status = 2),
        ..layer = list(model = "lnorm", status = 2,
                       args = c(meanlog = "param@D")),
        #> NO ..proba provided for status 3 and 4
        log_proba = "log(..proba)"
      ) |>
      infer()
  ))
  expect_match(err_msg, "Missing definition")
  expect_match(err_msg, "status 3, 4")
  
  # second generation of definitions
  err_msg <- as.character(expect_error(
    hhmodel |>
      load_augmented_list(
        ..proba = list(model = "norm", status = 1, min = 0,
                       args = c(mean = "param@D")),
        ..proba = list(expr = "..layer", status = c(2, 3, 4)),
        ..layer = list(model = "lnorm", status = c(2, 3),
                       args = c(meanlog = "param@D")),
        #> NO ..layer provided for status 3, 4
        log_proba = "log(..proba)"
      ) |>
      infer()
  ))
  expect_match(err_msg, "Missing definition")
  expect_match(err_msg, "status 4")
})

