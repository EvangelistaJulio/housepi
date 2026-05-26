library(testthat)

aux_add <- function(self, name,
                    by = NA_character_,
                    reference = NULL,
                    joint = FALSE,
                    strata = NULL,
                    params = NULL,
                    prior = NA_character_, ...,
                    min = NA_real_, max = NA_real_,
                    proposal = NA_character_,
                    proposal_scale = NA_real_,
                    prefix = NA_character_,
                    .warn_over = TRUE) {
  self$add_quo(name = name, by = by,
               reference = reference, joint = joint,
               strata = strata, params = params,
               prior = prior, rlang::enquos(...),
               min = rlang::enquo(min),
               max = rlang::enquo(max),
               proposal = proposal,
               proposal_scale = proposal_scale,
               prefix = prefix,
               .warn_over = .warn_over)
  invisible(self)
}

# a dummy dataset and relative information for factorization
dummy_df <- data.frame(house_id = 1,
                       age = c(rep("child", 2), rep("adult", 2)),
                       sex = c("M", "F", "M", "M"))
# corresponding factors levels
dummy_factor_info <- list(age = list(levels = c("child", "adult")),
                          sex = list(levels = c("M", "F")))
# dummy group rate logic table
dummy_strata <- data.frame(isolation = c(FALSE, TRUE),
                           param_name = c(".reference", "inf_isol"))

####################################
## GROUP RATES : HELPER FUNCTIONS ##
####################################

test_that("Group rates logic table computes combinations and counts", {
  tbl <- group_rate_logic_table(dummy_df, c("age", "sex"))

  # counts should sum to number of rows and include expected frequencies
  expect_equal(sum(tbl$n), nrow(dummy_df))
  expect_true(all(sort(tbl$n) == sort(c(2L, 1L, 1L))))
  # ensure columns are present
  expect_true(all(c("age", "sex", "n") %in% colnames(tbl)))
})

test_that("Associate param names decorates table and marks reference", {
  tbl <- group_rate_logic_table(dummy_df, c("age", "sex"))
  res <- associate_group_param_names(tbl, prefix = "mgs", reference = NULL)

  expect_true("param_name" %in% colnames(res))
  # exactly one .reference row (largest group)
  expect_equal(sum(res$param_name == ".reference"), 1L)
  # other param_names should start with provided prefix
  others <- res$param_name[res$param_name != ".reference"]
  expect_true(all(grepl("^mgs_", others)))
})

test_that("Assert-valid-logic-table validates input", {
  # must be a data.frame
  assert_valid_logic_table(list(a = 1)) |>
    expect_error(regexp = "data.frame")
  # must have param_name column
  assert_valid_logic_table(data.frame(x = 1)) |>
    expect_error(regexp = "param_name")
  # must have at least one row
  assert_valid_logic_table(data.frame(param_name = character(0))) |>
    expect_error(regexp = "at least one row")
  # reserved feature column names rejected
  assert_valid_logic_table(data.frame(param_name = "a", id = 1)) |>
    expect_error(regexp = "reserved")
})

#################################
## GROUP RATES : ACTIVE FIELDS ##
#################################

test_that("Active fields are read-only", {
  HousepiGroupRate$new() |>
    helper_test_read_only(
      "parameters", "data_required",
      "param_required", "tidy_table"
    )
})

##################################
## GROUP RATES : STANDARD CALLS ##
##################################

test_that("Group rates integrate with data and exposes to Cpp format", {
  grp <- HousepiGroupRate$new()

  # remove warns when removing unknown group rate
  expect_warning(grp$remove(list("does_not_exist"), warn_removal = TRUE))

  # implement should build logic_tables without errors
  aux_add(grp, "foxhound", by = c("age", "sex"), joint = TRUE)
  expect_silent(grp$implement(dummy_df, dummy_factor_info))
  expect_true("foxhound" %in% names(grp$logic_tables))
  tidy <- grp$tidy_table
  expect_true(is.data.frame(tidy))

  # mcmc setup is default (Cauchemez et al., 2009)
  default_mcmc_setup <- list(prior = "lognormal", meanlog = 0, sdlog = 3,
                             min = 0, max = 10, proposal = "lognormal")
  expect_equal(grp$parameters[[1]], default_mcmc_setup)

  # expose_to_cpp convert factor columns into integers according to factors
  # UPDATE: this will change when cpp will allow for characters
  cpp_tbl <- grp$expose_to_cpp(dummy_factor_info)
  expect_true(is.data.frame(cpp_tbl))
  # columns corresponding to features are factorized
  expect_true(all(c("age", "sex") %in% colnames(cpp_tbl)))
  expect_true(all(cpp_tbl$age %in% c(1L, 2L)))
  expect_true(all(cpp_tbl$sex %in% c(1L, 2L)))
})

test_that("Group rates sentinels work correctly", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("chain_binomial")

  # empty rate definition
  hhmodel |>
    add_group_rate("susceptibility", by = NA) |>
    expect_error("[Ii]nvalid rate definition")

  # wrong type for covariate column choice
  hhmodel |>
    add_group_rate("susceptibility", by = 4) |>
    expect_error("character")

  # not implemented group
  hhmodel |>
    add_group_rate("cobras", by = "sex") |>
    expect_error("unused group rate.*cobras")

  # incomplete custom prior definition
  hhmodel |>
    add_group_rate("susceptibility", by = "sex", max = 3) |>
    expect_error("'prior' must be specified")

  # FIXME: this should fire!
  # hhmodel |>
  #   add_group_rate("susceptibility", by = "age")
})

test_that("Group rates reference support behaves as expected", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("chain_binomial")

  # missing reference
  strata <- data.frame(sex = c("F", "M"),
                       param_name = c("ladies", "gentlemen"))
  hhmodel |>
    add_group_rate("susceptibility", strata = strata) |>
    expect_error("param_name as '.reference'")
  # explicit no reference
  hhmodel |>
    add_group_rate("susceptibility", reference = ".none", strata = strata) |>
    expect_no_error()
})

test_that("Parameter options work correctly", {
  grp <- HousepiGroupRate$new()
  aux_add(grp, "peace", by = c("age", "sex"),
          prior = "uniform", min = 0, max = 10, joint = TRUE)
  expect_silent(grp$implement(dummy_df, dummy_factor_info))
  expect_equal(grp$parameters[[1]],
               list(prior = "uniform", min = 0, max = 10,
                    proposal = NA_character_, proposal_scale = NA_real_))

  grp <- HousepiGroupRate$new()
  aux_add(grp, "walker", by = c("age", "sex"),
          prior = "gamma", shape = 2, scale = 1, joint = TRUE)
  expect_silent(grp$implement(dummy_df, dummy_factor_info))
  expect_equal(grp$parameters[[1]],
               list(prior = "gamma", shape = "2|N", scale = "1|N",
                    min = NA_real_, max = NA_real_,
                    proposal = NA_character_, proposal_scale = NA_real_))
})

test_that("Load list works as expected", {
  grp <- HousepiGroupRate$new()
  helper_test_aux_load(
    grp,
    whatever = list(by = c("age", "sex"), prefix = "ct",
                    prior = "uniform", min = 0, max = 10),
    susceptibility = list(by = c("age", "sex"),
                          prior = "lognormal", sdlog = 3)
  )
  expect_silent(grp$implement(dummy_df, dummy_factor_info))
})

#############################################
## GROUP RATES : CUSTOM STRATA DEFINITIONS ##
#############################################

test_that("Group rates deal with custom logic table correctly", {
  grp <- HousepiGroupRate$new()
  aux_add(grp, name = "infectivity", strata = dummy_strata,
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
  strata <- data.frame(age = c("child", "adult", "elderly"),
                       param_name = c(".reference", "sus_old", "sus_elderly"))
  overrides <- list(
    sus_old = 1.73,
    sus_elderly = list(value = 0.5, prior = "lognormal", sdlog = 1.2,
                       min = NA_real_, max = 10,
                       proposal = "lognormal", proposal_scale = 1)
  )
  aux_add(grp, name = "susceptibility", strata = strata,
          params = overrides, prior = "gamma", shape = 2, scale = 1)

  expect_length(grp$parameters, 2)
  # overrides used instead of raw params
  expect_equal(grp$parameters$sus_old, 1.73)
  expect_equal(grp$parameters$sus_elderly$prior, "lognormal")
  expect_equal(grp$parameters$sus_elderly$sdlog, 1.2)

  aux_add(grp, name = "susceptibility",
          strata = data.frame(sex = c("M", "F"),
                              param_name = c("smtg", ".reference")),
          params = list(unknown = 2.4)) |>
    expect_error("unknown")
})

test_that("Group rate implementation preserves strata entries", {
  grp <- HousepiGroupRate$new()

  # add a standard auto-mode rate
  aux_add(grp, "community", by = "age", prefix = "com")

  # add a custom-mode rate
  aux_add(grp, name = "infectivity",
          strata = data.frame(age = c("adult", "child"),
                              param_name = c("inf_old", ".reference")))

  # implement should build by-mode and preserve strata-mode
  grp$implement(dummy_df, dummy_factor_info)

  expect_true("infectivity" %in% names(grp$logic_tables))
  expect_true("community" %in% names(grp$logic_tables))
  # strata entry is unchanged
  grp$logic_tables$infectivity[[1]]$param_name |>
    expect_equal(c("inf_old", ".reference"))

  aux_add(grp, "infectivity", by = "age") |>
    expect_warning("infectivity")
  aux_add(grp, "community", strata = dummy_strata) |>
    expect_warning("community")
})

test_that("Group rate mixed auto and custom modes", {
  grp <- HousepiGroupRate$new()

  aux_add(grp, "community", by = "age", prefix = "com")
  grp$implement(dummy_df, dummy_factor_info)
  expect_true("community" %in% names(grp$logic_tables))

  aux_add(grp, name = "infectivity", strata = dummy_strata)
  expect_true("infectivity" %in% names(grp$logic_tables))

  # param_required should include both
  required <- grp$param_required
  expect_true("inf_isol" %in% required)
  expect_true(any(grepl("^com_", required)))

  # data required
  data_req <- grp$data_required
  expect_true(all(c("age", "isolation") %in% data_req))

  # trigger error (unknown feature in dummy_strata)
  grp$implement(dummy_df, dummy_factor_info) |>
    expect_error("Unknown")
})
