library(testthat)

aux_add <- function(self, name, by, prefix = NA_character_,
                    reference = NULL, joint = FALSE,
                    prior = NA_character_, ...,
                    min = NULL, max = NULL,
                    proposal = NA_character_,
                    proposal_scale = NA_real_) {
  self$add_quo(name = name, by = by, prefix = prefix,
               reference = reference,
               joint = joint,
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
  expect_equal(grp$parameters[[1]], group_rates_default_mcmc_setup)
})

test_that("Parameter options work correctly", {
  grp <- HousepiGroupRate$new()
  aux_add(grp, "contact", by = c("age_group", "sex"), prefix = "ct",
          prior = "uniform", min = 0, max = 10,
          joint = TRUE)
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
  # no duplicate parameters
  expect_error(assert_valid_logic_table(data.frame(param_name = rep("a", 2))),
               regexp = "Duplicate")
  # reserved feature column names rejected
  expect_error(assert_valid_logic_table(data.frame(param_name = "a", id = 1)),
               regexp = "reserved")
})

############################################################################
######################### Logic Table Mode Tests ###########################
############################################################################

test_that("HousepiGroupRate add_logic_table stores correctly", {
  grp <- HousepiGroupRate$new()
  lt <- data.frame(isolation = c(FALSE, TRUE),
                   param_name = c(".reference", "inf_isol"))
  grp$add_logic_table(name = "infectivity", logic_table = lt,
                      prior = "lognormal", meanlog = 0, sdlog = 2)

  # defs stored with mode sentinel
  expect_equal(grp$defs$infectivity$mode, "logic_table")
  expect_true(is.na(grp$defs$infectivity$by))

  # logic table stored
  expect_length(grp$logic_tables$infectivity, 1)
  stored <- grp$logic_tables$infectivity[[1]]
  expect_true("param_name" %in% colnames(stored))
  expect_true("n" %in% colnames(stored))
  expect_equal(stored$param_name, c(".reference", "inf_isol"))

  # parameters generated for non-reference
  params <- grp$parameters
  expect_length(params, 1)
  expect_true("inf_isol" %in% names(params))
})

test_that("HousepiGroupRate add_logic_table param_overrides work", {
  grp <- HousepiGroupRate$new()
  lt <- data.frame(age = c("young", "old", "elderly"),
                   param_name = c(".reference", "sus_old", "sus_elderly"))
  overrides <- list(
    sus_old = list(value = 1, prior = "uniform", min = 0, max = 5,
                   proposal = "normal", proposal_scale = 0.5),
    sus_elderly = list(value = 0.5, prior = "lognormal", meanlog = 0, sdlog = 1,
                       min = NA_real_, max = 10,
                       proposal = "lognormal", proposal_scale = 1)
  )
  grp$add_logic_table(
    name = "susceptibility", logic_table = lt,
    params = overrides, prior = "lognormal", meanlog = 0, sdlog = 3
  )

  params <- grp$parameters
  expect_length(params, 2)
  # overrides should be used, not shared raw_params
  expect_equal(params$sus_old$prior, "uniform")
  expect_equal(params$sus_elderly$prior, "lognormal")
  expect_equal(params$sus_elderly$sdlog, 1)
})

test_that("HousepiGroupRate add_logic_table append mode works", {
  grp <- HousepiGroupRate$new()
  lt1 <- data.frame(
    isolation = c(FALSE, TRUE),
    param_name = c(".reference", "inf_iso")
  )
  lt2 <- data.frame(
    vaccinated = c(FALSE, TRUE),
    param_name = c(".reference", "inf_vax")
  )
  grp$add_logic_table(name = "infectivity", logic_table = lt1)
  grp$add_logic_table(name = "infectivity", logic_table = lt2, append = TRUE)

  # two tables stored under the same rate
  expect_length(grp$logic_tables$infectivity, 2)

  # parameters from both tables
  params <- grp$parameters
  expect_true("inf_iso" %in% names(params))
  expect_true("inf_vax" %in% names(params))
})

test_that("HousepiGroupRate logic_table does not require data", {
  grp <- HousepiGroupRate$new()
  lt <- data.frame(
    isolation = c(FALSE, TRUE),
    param_name = c(".reference", "inf_iso")
  )
  grp$add_logic_table(name = "infectivity", logic_table = lt)

  # data_required should be empty
  expect_length(grp$data_required, 0)
})

test_that("HousepiGroupRate implement preserves logic_table entries", {
  grp <- HousepiGroupRate$new()

  # add a logic_table-mode rate
  lt <- data.frame(
    isolation = c(FALSE, TRUE),
    param_name = c(".reference", "inf_iso")
  )
  grp$add_logic_table(name = "infectivity", logic_table = lt)

  # add a standard by-mode rate
  aux_add(grp, "community", by = "age_group", prefix = "com")

  # implement should build by-mode and preserve logic_table-mode
  grp$implement(dummy_df, dummy_factor_info)

  expect_true("infectivity" %in% names(grp$logic_tables))
  expect_true("community" %in% names(grp$logic_tables))
  # logic_table entry is unchanged
  expect_equal(
    grp$logic_tables$infectivity[[1]]$param_name,
    c(".reference", "inf_iso")
  )
})

test_that("HousepiGroupRate remove cleans up logic_table entries", {
  grp <- HousepiGroupRate$new()
  lt <- data.frame(
    isolation = c(FALSE, TRUE),
    param_name = c(".reference", "inf_iso")
  )
  grp$add_logic_table(name = "infectivity", logic_table = lt)
  expect_true("infectivity" %in% names(grp$defs))

  grp$remove("infectivity")
  expect_false("infectivity" %in% names(grp$defs))
  expect_false("infectivity" %in% names(grp$logic_tables))
  expect_false("infectivity" %in% names(grp$param_overrides))
})

test_that("HousepiGroupRate tidy_table works with logic_table mode", {
  grp <- HousepiGroupRate$new()
  lt <- data.frame(
    isolation = c(FALSE, TRUE),
    param_name = c(".reference", "inf_iso")
  )
  grp$add_logic_table(name = "infectivity", logic_table = lt)
  tidy <- grp$tidy_table
  expect_true(is.data.frame(tidy))
  expect_true("group_rate" %in% colnames(tidy))
  expect_true(all(tidy$group_rate == "infectivity"))
})

test_that("HousepiGroupRate mixed mode: by + logic_table", {
  grp <- HousepiGroupRate$new()

  # standard by-mode
  aux_add(grp, "community", by = "age_group", prefix = "com")

  # logic_table-mode
  lt <- data.frame(
    isolation = c(FALSE, TRUE),
    param_name = c(".reference", "inf_iso")
  )
  grp$add_logic_table(name = "infectivity", logic_table = lt)

  # implement with data
  grp$implement(dummy_df, dummy_factor_info)

  # both rates should have logic tables
  expect_true("community" %in% names(grp$logic_tables))
  expect_true("infectivity" %in% names(grp$logic_tables))

  # param_required should include both
  required <- grp$param_required
  expect_true("inf_iso" %in% required)
  expect_true(any(grepl("^com_", required)))

  # data_required only from by-mode
  data_req <- grp$data_required
  expect_true("age_group" %in% data_req)
})

test_that("HousepiContactRate add_logic_table stores correctly", {
  cr <- HousepiContactRate$new()
  lt <- data.frame(
    from_age = c("adult", "adult", "child", "child"),
    to_age   = c("adult", "child", "adult", "child"),
    param_name = c(".reference", "cr_a2c", "cr_c2a", "cr_c2c")
  )
  cr$add_logic_table(logic_table = lt, prior = "lognormal", meanlog = 0,
                     sdlog = 3)

  # defs stored with mode sentinel
  expect_equal(cr$defs$mode, "logic_table")
  expect_true(is.na(cr$defs$by))

  # logic table stored
  expect_true(is.data.frame(cr$logic_table))
  expect_equal(nrow(cr$logic_table), 4)

  # parameters generated for non-reference

  params <- cr$parameters
  expect_length(params, 3)
  expect_true(all(c("cr_a2c", "cr_c2a", "cr_c2c") %in% names(params)))
})

test_that("HousepiContactRate add_logic_table param_overrides work", {
  cr <- HousepiContactRate$new()
  lt <- data.frame(
    from_age = c("adult", "child"),
    param_name = c(".reference", "cr_child")
  )
  overrides <- list(
    cr_child = list(value = 1, prior = "uniform", min = 0, max = 5,
                    proposal = "normal", proposal_scale = 0.5)
  )
  cr$add_logic_table(logic_table = lt, params = overrides)

  params <- cr$parameters
  expect_equal(params$cr_child$prior, "uniform")
  expect_equal(params$cr_child$min, 0)
})

test_that("HousepiContactRate logic_table skips implement", {
  cr <- HousepiContactRate$new()
  lt <- data.frame(
    from_age = c("adult", "child"),
    param_name = c(".reference", "cr_child")
  )
  cr$add_logic_table(logic_table = lt)

  # implement should not clear the logic table
  cr$implement(dummy_df, dummy_factor_info)
  expect_true(is.data.frame(cr$logic_table))
  expect_equal(cr$logic_table$param_name, c(".reference", "cr_child"))
})

test_that("HousepiContactRate tidy_table works with logic_table mode", {
  cr <- HousepiContactRate$new()
  lt <- data.frame(
    from_age = c("adult", "child"),
    param_name = c(".reference", "cr_child")
  )
  cr$add_logic_table(logic_table = lt)

  tidy <- cr$tidy_table
  expect_true(is.data.frame(tidy))
  expect_equal(nrow(tidy), 2)
})

############################################################################

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
