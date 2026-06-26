library(testthat)

test_that("Reject if calling on non-HousepiModel objects", {
  obj <- "not a model"
  err <- "class 'HousepiModel'"
  expect_error(add_data_long(obj, data.frame()), regexp = err)
  expect_error(get_data_wide(obj), regexp = err)
  expect_error(get_data_long(obj), regexp = err)
  expect_error(get_parameters(obj), regexp = err)
  expect_error(get_augmented(obj), regexp = err)
  expect_error(get_group_rates(obj), regexp = err)
  expect_error(get_infec_profile(obj), regexp = err)
  expect_error(drop_parameters(obj), regexp = err)
  expect_error(drop_augmented(obj), regexp = err)
  expect_error(drop_observations(obj), regexp = err)
  expect_error(drop_house_observ(obj), regexp = err)
  expect_error(drop_group_rates(obj), regexp = err)
  expect_error(drop_contact_rates(obj), regexp = err)
  expect_error(add_contact_rate(obj), regexp = err)
  expect_error(drop_infectivity_profile(obj), regexp = err)
  expect_error(drop_all(obj), regexp = err)
  expect_error(add_parameter(obj), regexp = err)
  expect_error(add_augmented(obj), regexp = err)
  expect_error(add_infec_profile(obj, "mk22"), regexp = err)
  expect_error(evaluate_log_lik(obj), regexp = err)
  expect_error(run_inference(obj), regexp = err)
  expect_error(run_simulation(obj), regexp = err)
  expect_error(add_parameter_template(obj), regexp = err)
  expect_error(add_augmented_template(obj), regexp = err)
  expect_error(add_group_rate(obj), regexp = err)
  expect_error(predicted_measures(obj, NULL), regexp = err)
  expect_error(predicted_cases_by_house_size(obj, NULL), regexp = err)
})

test_that("data_wide provided at construction is exposed via get_data_wide", {
  result <- new_housepi_model(flu_1918_baltimora, "blank")
  dw <- get_data_wide(result)
  expect_true(is.data.frame(dw))
  expect_true("house_id" %in% colnames(dw))
  expect_true("status" %in% colnames(dw))

  # id automatically added
  expect_identical(dw$id, as.character(seq_len(dim(dw)[1])))
  # house_id cast to character
  expect_identical(dw$house_id, as.character(flu_1918_baltimora$house_id))
  # statuses are integer
  expect_identical(dw$status, flu_1918_baltimora$status)
})

test_that("add_data_long and get_data_long work", {
  hhmodel <- new_housepi_model(model = "blank")
  # blank model has default data with id = "tommy", "gina"
  dl <- data.frame(id = hhmodel$data_wide$id,
                   measure = c(1.5, 2.5))
  result <- add_data_long(hhmodel, dl)
  out <- get_data_long(result)
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 2)
})

test_that("Wrappers work as expected on blank model", {
  hhmodel <- new_housepi_model(model = "blank")

  # parameters
  result <- add_parameter(hhmodel, name = "alpha", value = 0.5)
  pars <- get_parameters(result)
  expect_true(is.data.frame(pars))
  expect_true("alpha" %in% pars$name)

  result <- add_parameter(hhmodel, alpha = 1, beta = 2)
  pars <- get_parameters(result)
  expect_true(all(c("alpha", "beta") %in% pars$name))

  # infectivity profile
  expect_null(get_infec_profile(hhmodel))

  # augmented
  result <- add_augmented(hhmodel, name = "dummy", value = 3)
  aug <- get_augmented(result)
  expect_true(is.data.frame(aug))
  expect_true("dummy" %in% aug$name)

  result <- hhmodel |>
    add_augmented(dummy = list(prior = "uniform"))
  aug <- get_augmented(result)
  expect_true("dummy" %in% aug$name)
})

test_that("Bulk loading via `defs` mirrors inline name = value pairs", {
  hhmodel <- new_housepi_model(model = "blank")
  defs <- list(metal = 1, gear = list(prior = "uniform"))

  # an external list passed via `defs` equals the inline dotted form
  from_defs <- hhmodel |>
    add_parameter(import = defs) |>
    get_parameters()
  inline <- hhmodel |>
    add_parameter(metal = 1, gear = list(prior = "uniform")) |>
    get_parameters()
  expect_identical(from_defs, inline)

  # `defs` cannot be combined with a single name or with inline pairs
  hhmodel |>
    add_parameter("sokolov", import = defs) |>
    expect_error("Mixing")
  hhmodel |>
    add_parameter(nikolai = 1, import = defs) |>
    expect_error("Mixing")
})

############################################################################
## Drop wrappers
############################################################################

test_that("Drop wrappers return clean models", {
  hhmodel <- new_housepi_model(model = "blank") |>
    add_parameter(name = "alpha", value = 0.5) |>
    add_augmented(name = "dummy", value = 3)

  # drop_parameters returns new object without parameters
  dropped <- drop_parameters(hhmodel)
  expect_equal(nrow(get_parameters(dropped)), 0L)
  # original is unchanged

  expect_true("alpha" %in% get_parameters(hhmodel)$name)

  # drop_augmented
  dropped <- drop_augmented(hhmodel)
  expect_equal(nrow(get_augmented(dropped)), 0L)
  expect_true("dummy" %in% get_augmented(hhmodel)$name)

  # drop_group_rates
  expect_no_error(drop_group_rates(hhmodel))
  # drop_contact_rates
  expect_no_error(drop_contact_rates(hhmodel))

  # drop_infectivity_profile
  expect_no_error(drop_infectivity_profile(hhmodel))

  # drop_all
  dropped <- drop_all(hhmodel)
  expect_equal(nrow(get_parameters(dropped)), 0)
  expect_equal(nrow(get_augmented(dropped)), 0)
})

############################################################################
## Wrappers with data
############################################################################

test_that("Group rate wrappers", {
  small_gr <- data.frame(
    house_id = c(1, 1, 2, 2, 2),
    status = c(1, 0, 1, 0, 0),
    member = c("student", "other", "student", "mother", "father")
  )
  hhmodel <- new_housepi_model(small_gr, "chain_binomial")
  hhmodel <- add_group_rate(hhmodel, "community", by = "member")
  groups_1 <- get_group_rates(hhmodel)
  expect_true(is.data.frame(groups_1))
  expect_true(nrow(groups_1) > 0)

  expect_warning(
    result <- add_group_rate(hhmodel, community = "member"),
    "[Oo]verwriting.*community"
  )
  groups_2 <- get_group_rates(result)
  expect_identical(groups_1, groups_2)

  # load custom definitions
  hhmodel <- hhmodel |>
    add_group_rate(
      susceptibility = list(
        strata = data.frame(member = c("student", NA),
                            param_name = c("sus_student", ".reference")),
        params = list(sus_student = 1.23)
      )
    )
  expect_equal(subset(hhmodel$parameters, name == "sus_student")$value, 1.23)

  # check for non-existing custom features
  hhmodel <- drop_group_rates(hhmodel)
  hhmodel |>
    add_group_rate(
      "susceptibility",
      strata = data.frame(member = c("kaiju", NA),
                          param_name = c("sus_student", ".reference"))
    ) |>
    expect_warning("not found")

  hhmodel$group_rates
})

test_that("add_contact_rate wrapper input branches", {
  hhmodel <- data.frame(
    house_id = c(1, 1, 2, 2, 2),
    status = c(1, 0, 1, 0, 0),
    member = c("student", "other", "student", "mother", "father")
  ) |>
    new_housepi_model("chain_binomial")

  strata <- data.frame(
    from_member = c("student", NA),
    to_member = c("student", NA),
    param_name = c("ctr_kept_you_waiting", ".reference")
  )

  # 'by' and 'strata' together is ambiguous
  hhmodel |>
    add_contact_rate(by = "member", strata = strata) |>
    expect_error("[Ii]mpossible")

  # neither 'by' nor 'strata' is empty
  hhmodel |>
    add_contact_rate() |>
    expect_error("[Ii]nvalid")

  # first call succeeds and returns a NEW model (copy() semantic)
  result <- hhmodel |>
    add_contact_rate(strata = strata,
                     params = list(ctr_kept_you_waiting = 1.21))
  expect_false(identical(hhmodel, result))
  # original is untouched (tidy contact-rate table stays empty)
  expect_equal(nrow(hhmodel$contact_rates), 0L)
  expect_true(nrow(result$contact_rates) > 0L)
  # rebuild_ ran: param is exposed on the returned model
  expect_true("ctr_kept_you_waiting" %in% get_parameters(result)$name)

  # second call warns about overwrite
  expect_warning(
    add_contact_rate(result, strata = strata,
                     params = list(ctr_kept_you_waiting = 1.21)),
    regexp = "[Oo]verwriting"
  )

  # getter
  cr <- get_contact_rates(result)
  expect_equal(cr, result$contact_rates)
  expect_true(is.data.frame(cr))
  expect_true(nrow(cr) > 0)
})

test_that("add_infec_profile single and bulk forms work", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("augmented_infection_time")
  ip <- hhmodel |>
    add_infec_profile(func = "lognormal") |>
    get_infec_profile()
  expect_false(is.null(ip))
  expect_true("lognormal" %in% ip$func)

  # bulk form
  ip2 <- hhmodel |>
    add_infec_profile(ip1 = list(func = "lognormal")) |>
    get_infec_profile()
  expect_false(is.null(ip2))
  expect_true("lognormal" %in% ip2$func)
})
