library(testthat)

aux_add <- function(self, name,
                    by = NA_character_,
                    prefix = NA_character_,
                    reference = NULL,
                    joint = FALSE,
                    logic_table = NULL,
                    append = FALSE,
                    params = NULL,
                    prior = NA_character_, ...,
                    min = NA_real_, max = NA_real_,
                    proposal = NA_character_,
                    proposal_scale = NA_real_,
                    .warn_over = TRUE) {
  self$add_quo(name = name, by = by, prefix = prefix,
               reference = reference, joint = joint,
               logic_table = logic_table,
               append = append, params = params,
               prior = prior, rlang::enquos(...),
               min = rlang::enquo(min),
               max = rlang::enquo(max),
               proposal = proposal,
               proposal_scale = proposal_scale,
               .warn_over = .warn_over)
  invisible(self)
}

# a dummy dataset and relative information for factorization
dummy_df <- data.frame(age_group = c("child", "child", "adult", "child"),
                       sex = c("M", "F", "M", "M"),
                       stringsAsFactors = FALSE)
dummy_factor_info <- list(age_group = list(levels = c("child", "adult")),
                          sex = list(levels = c("M", "F")))

# dummy group rate logic table
dummy_logic_table <- data.frame(isolation = c(FALSE, TRUE),
                                param_name = c(".reference", "inf_isol"))

test_that("Group rates logic table computes combinations and counts", {
  tbl <- group_rate_logic_table(dummy_df, c("age_group", "sex"))

  # counts should sum to number of rows and include expected frequencies
  expect_equal(sum(tbl$n), nrow(dummy_df))
  expect_true(all(sort(tbl$n) == sort(c(2L, 1L, 1L))))
  # ensure columns are present
  expect_true(all(c("age_group", "sex", "n") %in% colnames(tbl)))
})

test_that("Remove warns when removing unknown group rate", {
  grp <- HousepiGroupRate$new()
  expect_warning(grp$remove(list("does_not_exist"), warn_removal = TRUE))
})

test_that("Associate param names decorates table and marks reference", {
  tbl <- group_rate_logic_table(dummy_df, c("age_group", "sex"))
  res <- associate_param_names(tbl, prefix = "g", reference = NULL)

  expect_true("param_name" %in% colnames(res))
  # exactly one .reference row (largest group)
  expect_equal(sum(res$param_name == ".reference"), 1L)
  # other param_names should start with provided prefix
  others <- res$param_name[res$param_name != ".reference"]
  expect_true(all(grepl("^g_", others)))
})

test_that("Group rates integrate with data and exposes to Cpp format", {
  grp <- HousepiGroupRate$new()
  aux_add(grp, "contact", by = c("age_group", "sex"), prefix = "ct",
          joint = TRUE)
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
  default_mcmc_setup <- list(value = NA_real_, prior = "lognormal",
                             meanlog = 0, sdlog = 3, min = 0, max = 10,
                             proposal = "lognormal", proposal_scale = 1)
  expect_equal(grp$parameters[[1]], default_mcmc_setup)
})

test_that("Parameter options work correctly", {
  grp <- HousepiGroupRate$new()
  aux_add(grp, "contact", by = c("age_group", "sex"), prefix = "ct",
          prior = "uniform", min = 0, max = 10, joint = TRUE)
  expect_silent(grp$implement(dummy_df, dummy_factor_info))
  expect_equal(grp$parameters[[1]],
               list(value = 1, prior = "uniform", min = 0, max = 10,
                    proposal = NA_character_, proposal_scale = NA_real_))

  grp <- HousepiGroupRate$new()
  aux_add(grp, "contact", by = c("age_group", "sex"), prefix = "ct",
          prior = "gamma", shape = 2, scale = 1, joint = TRUE)
  expect_silent(grp$implement(dummy_df, dummy_factor_info))
  expect_equal(grp$parameters[[1]],
               list(value = 1, prior = "gamma", shape = "2|N", scale = "1|N",
                    min = NA_real_, max = NA_real_,
                    proposal = NA_character_, proposal_scale = NA_real_))
})

test_that("Assert-valid-logic-table validates input", {
  # must be a data.frame
  expect_error(assert_valid_logic_table(list(a = 1)),
               regexp = "data.frame")
  # must have param_name column
  expect_error(assert_valid_logic_table(data.frame(x = 1)),
               regexp = "param_name")
  # must have at least one row
  expect_error(assert_valid_logic_table(data.frame(param_name = character(0))),
               regexp = "at least one row")
  # reserved feature column names rejected
  expect_error(assert_valid_logic_table(data.frame(param_name = "a", id = 1)),
               regexp = "reserved")
})

test_that("Group rates deal with custom logic table correctly", {
  grp <- HousepiGroupRate$new()
  aux_add(grp, name = "infectivity", logic_table = dummy_logic_table,
          prior = "lognormal", meanlog = 0, sdlog = 2)

  # defs stored with mode sentinel
  expect_equal(grp$defs$infectivity$mode, "custom")
  expect_equal(grp$defs$infectivity$by, "isolation")

  # logic table stored
  expect_length(grp$logic_tables$infectivity, 1)
  stored <- grp$logic_tables$infectivity[[1]]
  expect_true("param_name" %in% colnames(stored))
  expect_equal(stored$param_name, c(".reference", "inf_isol"))

  # parameters generated for non-reference
  params <- grp$parameters
  expect_length(params, 1)
  expect_true("inf_isol" %in% names(params))

  # tidy_table
  tidy <- grp$tidy_table
  expect_true(is.data.frame(tidy))
  expect_true("group_rate" %in% colnames(tidy))
  expect_true(all(tidy$group_rate == "infectivity"))

  # remove logic
  grp$remove("infectivity")
  expect_false("infectivity" %in% names(grp$defs))
  expect_false("infectivity" %in% names(grp$logic_tables))
  expect_false("infectivity" %in% names(grp$param_overrides))
})

test_that("Group rate parameter overrides work", {
  grp <- HousepiGroupRate$new()
  lt <- data.frame(age = c("child", "adult", "elderly"),
                   param_name = c(".reference", "sus_old", "sus_elderly"))
  overrides <- list(
    sus_old = 1.73,
    sus_elderly = list(value = 0.5, prior = "lognormal", sdlog = 1.2,
                       min = NA_real_, max = 10,
                       proposal = "lognormal", proposal_scale = 1)
  )
  aux_add(grp, name = "susceptibility", logic_table = lt,
          params = overrides, prior = "gamma", shape = 2, scale = 1)

  expect_length(grp$parameters, 2)
  # overrides used instead of raw params
  expect_equal(grp$parameters$sus_old, 1.73)
  expect_equal(grp$parameters$sus_elderly$prior, "lognormal")
  expect_equal(grp$parameters$sus_elderly$sdlog, 1.2)

  # FIXME: params 'unknown' is not introduced in this definition
  expect_error(aux_add(
    grp, name = "susceptibility",
    logic_table = data.frame(sex = c("M", "F"),
                             param_name = c("something", ".reference")),
    append = TRUE, params = list(unknown = 2.4)
  ), regexp = "unknown")
})

test_that("Group Rate append logic table in custom mode works", {
  grp <- HousepiGroupRate$new()
  lt1 <- dummy_logic_table
  lt2 <- data.frame(vaccinated = c(FALSE, TRUE),
                    param_name = c(".reference", "inf_vax"))
  aux_add(grp, name = "infectivity", logic_table = lt1)
  aux_add(grp, name = "infectivity", logic_table = lt2, append = TRUE)

  # two tables stored under the same rate
  expect_length(grp$logic_tables$infectivity, 2)

  # parameters from both tables
  params <- grp$parameters
  expect_true("inf_isol" %in% names(params))
  expect_true("inf_vax" %in% names(params))

  expect_true(all(c("isolation", "vaccinated") %in% grp$data_required))
})

test_that("Group rate implementation preserves logic_table entries", {
  grp <- HousepiGroupRate$new()

  # add a standard auto-mode rate
  aux_add(grp, "community", by = "age_group", prefix = "com")

  # add a custom-mode rate
  aux_add(grp, name = "infectivity",
          logic_table = data.frame(age_group = c("adult", "child"),
                                   param_name = c("inf_old", ".reference")))

  # implement should build by-mode and preserve logic_table-mode
  grp$implement(dummy_df, dummy_factor_info)

  expect_true("infectivity" %in% names(grp$logic_tables))
  expect_true("community" %in% names(grp$logic_tables))
  # logic_table entry is unchanged
  expect_equal(grp$logic_tables$infectivity[[1]]$param_name,
               c("inf_old", ".reference"))

  expect_warning(aux_add(grp, "infectivity", by = "age_group"),
                 regexp = "infectivity")
  expect_warning(aux_add(grp, "community", logic_table = dummy_logic_table),
                 regexp = "community")
})

test_that("Group rate mixed auto and custom modes", {
  grp <- HousepiGroupRate$new()

  aux_add(grp, "community", by = "age_group", prefix = "com")
  grp$implement(dummy_df, dummy_factor_info)
  expect_true("community" %in% names(grp$logic_tables))

  aux_add(grp, name = "infectivity", logic_table = dummy_logic_table)
  expect_true("infectivity" %in% names(grp$logic_tables))

  # param_required should include both
  required <- grp$param_required
  expect_true("inf_isol" %in% required)
  expect_true(any(grepl("^com_", required)))

  # data required
  data_req <- grp$data_required
  expect_true(all(c("age_group", "isolation") %in% data_req))

  # trigger error (unknown feature in dummy_logic_table)
  expect_error(grp$implement(dummy_df, dummy_factor_info),
               regexp = "Unknown")
})

test_that("Load list works as expected", {
  grp <- HousepiGroupRate$new()
  helper_test_aux_load(
    grp,
    contact = list(by = c("age_group", "sex"), prefix = "ct",
                   prior = "uniform", min = 0, max = 10),
    susceptibility = list(by = c("age_group", "sex"),
                          prior = "lognormal", sdlog = 3)
  )
  expect_silent(grp$implement(dummy_df, dummy_factor_info))
  grp$parameters
})
