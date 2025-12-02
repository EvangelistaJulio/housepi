library(testthat)
library(housepi)

aux_add <- function(self, name, by, prefix = NA_character_,
                    reference = NULL, independent = TRUE,
                    prior = NA_character_, ...,
                    min = NULL, max = NULL,
                    proposal = NA_character_,
                    proposal_scale = NA_real_) {
  self$add_quo(name = name, by = by, prefix = prefix,
               reference = reference,
               independent = independent,
               prior = prior, rlang::enquos(...),
               min = rlang::enquo(min),
               max = rlang::enquo(max),
               proposal = proposal,
               proposal_scale = proposal_scale)
  invisible(self)
}

# a dummy dataset and relative information for factorization
dummy_df <- data.frame(age_group = c("young", "young", "old", "young"),
                       sex = c("M", "F", "M", "M"),
                       stringsAsFactors = FALSE)
dummy_factor_info <- list(age_group = list(levels = c("young", "old")),
                          sex = list(levels = c("M", "F")))

test_that("Group rates logic table computes combinations and counts", {

  tbl <- housepi:::group_rate_logic_table(dummy_df, by = c("age_group", "sex"))

  # counts should sum to number of rows and include expected frequencies
  expect_equal(sum(tbl$n), nrow(dummy_df))
  expect_true(all(sort(tbl$n) == sort(c(2L, 1L, 1L))))
  # ensure columns are present
  expect_true(all(c("age_group", "sex", "n") %in% colnames(tbl)))
})

test_that("Remove warns when removing unknown group rate", {
  grp <- housepi:::HousepiGroupRate$new()
  expect_warning(grp$remove(list("does_not_exist"), warn_removal = TRUE))
})

test_that("Associate param names decorates table and marks reference", {
  tbl <- housepi:::group_rate_logic_table(dummy_df, by = c("age_group", "sex"))
  res <- housepi:::associate_param_names(tbl, prefix = "g", reference = NULL)

  expect_true("param_name" %in% colnames(res))
  # exactly one .reference row (largest group)
  expect_equal(sum(res$param_name == ".reference"), 1L)
  # other param_names should start with provided prefix
  others <- res$param_name[res$param_name != ".reference"]
  expect_true(all(grepl("^g_", others)))
})

test_that("Group rates integrate with data and exposes to Cpp format", {
  grp <- housepi:::HousepiGroupRate$new()
  aux_add(grp, "contact", by = c("age_group", "sex"), prefix = "ct",
          independent = FALSE)
  # implement should build logic_tables without errors
  expect_silent(grp$implement(dummy_df, dummy_factor_info))

  expect_true("contact" %in% names(grp$logic_tables))
  tidy <- grp$tidy_table
  expect_true(is.data.frame(tidy))
  # expose_to_cpp convert factor columns into integers according to factors
  cpp_tbl <- grp$expose_to_cpp(dummy_factor_info)
  expect_true(is.data.frame(cpp_tbl))

  # columns corresponding to features are factorized
  expect_true(all(c("age_group", "sex") %in% colnames(cpp_tbl)))
  expect_true(all(cpp_tbl$age_group %in% c(1L, 2L)))
  expect_true(all(cpp_tbl$sex %in% c(1L, 2L)))

  # mcmc setup is default
  expect_equal(grp$parameters[[1]], housepi:::group_rates_default_mcmc_setup)
})

test_that("Parameter options work correctly", {
  grp <- HousepiGroupRate$new()
  aux_add(grp, "contact", by = c("age_group", "sex"), prefix = "ct",
          prior = "uniform", min = 0, max = 10,
          independent = FALSE)
  expect_silent(grp$implement(dummy_df, dummy_factor_info))
  expect_equal(grp$parameters[[1]],
               list(value = 1, prior = "uniform", min = 0, max = 10,
                    proposal = NA_character_, proposal_scale = NA_real_))

  grp <- HousepiGroupRate$new()
  aux_add(grp, "contact", by = c("age_group", "sex"), prefix = "ct",
          prior = "gamma", shape = 2, scale = 1, independent = FALSE)
  expect_silent(grp$implement(dummy_df, dummy_factor_info))
  expect_equal(grp$parameters[[1]],
               list(value = 1, prior = "gamma", shape = 2, scale = 1,
                    min = NA_real_, max = NA_real_,
                    proposal = NA_character_, proposal_scale = NA_real_))
})
