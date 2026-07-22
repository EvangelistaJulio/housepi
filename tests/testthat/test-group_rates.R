library(testthat)

# a dummy dataset for the plain logic-table helpers
dummy_df <- data.frame(house_id = 1,
                       age = c(rep("child", 2), rep("adult", 2)),
                       sex = c("M", "F", "M", "M"))
# dummy group rate logic table
dummy_strata <- data.frame(isolation = c(FALSE, TRUE),
                           param_name = c(".reference", "inf_isol"))

# Model over the dummy features. Feature columns stay character on purpose:
# levels follow order of first appearance, which fixes the reference stratum.
mk_rate_model <- function(age = c("child", "child", "adult", "adult")) {
  data.frame(house_id = 1, status = 0,
             age = age,
             sex = c("M", "F", "M", "M"),
             isolation = c(FALSE, FALSE, TRUE, TRUE)) |>
    new_housepi_model("final_size")
}

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

test_that("Custom strata logic table is validated", {
  # must be a data.frame
  mk_rate_model() |>
    add_group_rate("susceptibility", strata = list(a = 1)) |>
    expect_error("data\\.frame")
  # must have param_name column
  mk_rate_model() |>
    add_group_rate("susceptibility", strata = data.frame(x = 1)) |>
    expect_error("param_name")
  # must have at least one row
  mk_rate_model() |>
    add_group_rate("susceptibility",
                   strata = data.frame(param_name = character(0))) |>
    expect_error("at least one row")
  # reserved feature column names rejected
  mk_rate_model() |>
    add_group_rate("susceptibility",
                   strata = data.frame(param_name = "a", id = 1)) |>
    expect_error("[Rr]eserved")
})

##################################
## GROUP RATES : STANDARD CALLS ##
##################################

test_that("Group rates integrate with data", {
  hhmodel <- mk_rate_model() |>
    add_group_rate("community", by = c("age", "sex"), joint = TRUE) |>
    expect_silent()

  tidy <- get_group_rates(hhmodel)
  expect_true(is.data.frame(tidy))
  expect_true(all(tidy$group_rate == "community"))
  # joint expansion: one rule group over both features
  expect_equal(length(unique(tidy$rule_group)), 1L)
  expect_true(all(c("age", "sex") %in% colnames(tidy)))
  expect_equal(sum(tidy$param_name == ".reference"), 1L)

  # mcmc setup is default (Cauchemez et al., 2009)
  params <- get_parameters(hhmodel)
  expect_true(all(params$prior == "lognormal"))
  expect_equal(params$args[[1]], c(meanlog = 0, sdlog = 3))
  expect_true(all(params$min == 0))
  expect_true(all(params$max == 10))
  expect_true(all(params$proposal == "rw_log"))
})

test_that("Group rates hold without data and re-derive on add_data_wide", {
  # rates ingest data: a data-less model holds them, no rate params yet
  metalgear <- new_housepi_model(model = "final_size") |>
    add_parameter_template() |>
    add_group_rate("susceptibility", by = "member")
  expect_silent(metalgear)
  expect_false(any(grepl("^sus_member", metalgear$parameters$name)))

  # attaching data materialises the rate parameters and runs
  diamond <- add_data_wide(metalgear, flu_2014_matsumoto)
  full_params <- diamond$parameters$name
  expect_true(any(grepl("^sus_member", full_params)))
  expect_true(is.finite(evaluate_log_lik(diamond)))

  # replacing data re-derives: a single-level subset drops the stale rate params
  student <- flu_2014_matsumoto[flu_2014_matsumoto$member == "student", ]
  mother_base <- add_data_wide(diamond, student)
  expect_false(any(grepl("^sus_member", mother_base$parameters$name)))
  expect_true(is.finite(evaluate_log_lik(mother_base)))

  # re-attaching the full data restores the exact original parameter set
  restored <- add_data_wide(mother_base, flu_2014_matsumoto)
  expect_setequal(restored$parameters$name, full_params)
})

test_that("Group rates sentinels work correctly", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("final_size")

  # empty rate definition
  hhmodel |>
    add_group_rate("susceptibility", by = NA) |>
    expect_error("[Ii]nvalid rate definition")

  # wrong type for covariate column choice
  hhmodel |>
    add_group_rate("susceptibility", by = 4) |>
    expect_error("[Aa]rgument .by. must be of type character")

  # not implemented group
  hhmodel |>
    add_group_rate("cobras", by = "sex") |>
    expect_error("[Uu]nused group rate.*cobras")

  # incomplete custom prior definition
  hhmodel |>
    add_group_rate("susceptibility", by = "sex", max = 3) |>
    expect_error("[Ss]pare prior arguments")

  hhmodel |>
    add_group_rate("susceptibility", by = "test_date") |>
    expect_error("[Mm]ust be a scalar")
})

test_that("Group rates reference support behaves as expected", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("final_size")

  # missing reference
  strata <- data.frame(sex = c("F", "M"),
                       param_name = c("ladies", "gentlemen"))
  hhmodel |>
    add_group_rate("susceptibility", strata = strata) |>
    expect_error("param_name as")
  # explicit no reference
  hhmodel |>
    add_group_rate("susceptibility", reference = ".none", strata = strata) |>
    expect_no_error()
})

test_that("Parameter options work correctly", {
  params <- mk_rate_model() |>
    add_group_rate("susceptibility", by = c("age", "sex"),
                   prior = "uniform", min = 0, max = 10, joint = TRUE) |>
    get_parameters()
  expect_true(all(params$prior == "uniform"))
  expect_true(all(params$min == 0))
  expect_true(all(params$max == 10))
  # a sampled rate parameter defaults to a random walk
  expect_true(all(params$proposal == "rw"))

  params <- mk_rate_model() |>
    add_group_rate("infectivity", by = c("age", "sex"),
                   prior = "gamma", shape = 2, scale = 1, joint = TRUE) |>
    get_parameters()
  expect_true(all(params$prior == "gamma"))
  expect_equal(params$args[[1]], c(shape = 2, scale = 1))
  expect_true(all(is.na(params$min)))
  expect_true(all(is.na(params$max)))
  expect_true(all(params$proposal == "rw"))
})

test_that("Load list works as expected", {
  hhmodel <- mk_rate_model() |>
    add_group_rate(
      community = list(by = c("age", "sex"), prefix = "ct",
                       prior = "uniform", min = 0, max = 10),
      susceptibility = list(by = c("age", "sex"),
                            prior = "lognormal", sdlog = 3)
    ) |>
    expect_silent()

  tidy <- get_group_rates(hhmodel)
  expect_setequal(unique(tidy$group_rate), c("community", "susceptibility"))
  expect_true(any(grepl("^ct_", tidy$param_name)))
  expect_true(any(grepl("^sus_", tidy$param_name)))
})

#############################################
## GROUP RATES : CUSTOM STRATA DEFINITIONS ##
#############################################

test_that("Group rates deal with custom logic table correctly", {
  hhmodel <- mk_rate_model() |>
    add_group_rate(name = "infectivity", strata = dummy_strata,
                   prior = "lognormal", meanlog = 0, sdlog = 2)

  tidy <- get_group_rates(hhmodel)
  expect_true(is.data.frame(tidy))
  expect_true(all(tidy$group_rate == "infectivity"))
  # the custom strata drives the feature column and the rows
  expect_true("isolation" %in% colnames(tidy))
  expect_setequal(tidy$param_name, c(".reference", "inf_isol"))

  # parameters generated for non-reference only
  params <- get_parameters(hhmodel)
  expect_equal(params$name, "inf_isol")
  expect_equal(params$prior, "lognormal")
  expect_equal(params$args[[1]], c(meanlog = 0, sdlog = 2))
})

test_that("Group rate parameter overrides work", {
  strata <- data.frame(age = c("child", "adult", "elderly"),
                       param_name = c(".reference", "sus_old", "sus_elderly"))
  overrides <- list(
    sus_old = 1.73,
    sus_elderly = list(value = 0.5, prior = "lognormal", sdlog = 1.2,
                       min = NA_real_, max = 10,
                       proposal = "rw_log", proposal_scale = 1)
  )
  hhmodel <- mk_rate_model(c("child", "adult", "elderly", "adult")) |>
    add_group_rate(name = "susceptibility", strata = strata,
                   params = overrides, prior = "gamma", shape = 2, scale = 1)

  params <- get_parameters(hhmodel)
  expect_setequal(params$name, c("sus_old", "sus_elderly"))
  # overrides used instead of raw params
  old <- subset(params, name == "sus_old")
  expect_equal(old$value, 1.73)
  expect_true(old$const)
  elderly <- subset(params, name == "sus_elderly")
  expect_equal(elderly$prior, "lognormal")
  expect_equal(elderly$args[[1]], c(sdlog = 1.2))

  hhmodel |>
    add_group_rate(name = "susceptibility",
                   strata = data.frame(
                     sex = c("M", "F"), param_name = c("smtg", ".reference")
                   ),
                   params = list(unknown = 2.4)) |>
    expect_error("[Uu]nknown")
})

test_that("Group rate implementation preserves strata entries", {
  hhmodel <- mk_rate_model() |>
    # add a standard auto-mode rate
    add_group_rate("community", by = "age", prefix = "com") |>
    # add a custom-mode rate
    add_group_rate(name = "infectivity",
                   strata = data.frame(
                     age = c("adult", "child"),
                     param_name = c("inf_old", ".reference")
                   ))

  tidy <- get_group_rates(hhmodel)
  expect_setequal(unique(tidy$group_rate), c("community", "infectivity"))
  # strata entry is unchanged
  subset(tidy, group_rate == "infectivity")$param_name |>
    expect_setequal(c("inf_old", ".reference"))
  expect_true(any(grepl("^com_", tidy$param_name)))

  add_group_rate(hhmodel, "infectivity", by = "age") |>
    expect_warning("[Oo]verwriting.*infectivity")
  add_group_rate(hhmodel, "community", strata = dummy_strata) |>
    expect_warning("[Oo]verwriting.*community")
})

test_that("Group rate mixed auto and custom modes", {
  hhmodel <- mk_rate_model() |>
    add_group_rate("community", by = "age", prefix = "com") |>
    add_group_rate(name = "infectivity", strata = dummy_strata)

  tidy <- get_group_rates(hhmodel)
  expect_setequal(unique(tidy$group_rate), c("community", "infectivity"))
  # both modes materialise their parameters
  param_names <- get_parameters(hhmodel)$name
  expect_true("inf_isol" %in% param_names)
  expect_true(any(grepl("^com_", param_names)))

  # a feature the data does not carry only fails at the wall
  data.frame(house_id = 1, status = 0, age = c("child", "adult")) |>
    new_housepi_model("final_size") |>
    add_group_rate(name = "infectivity", strata = dummy_strata) |>
    evaluate_log_lik() |>
    expect_error("missing data column")
})
