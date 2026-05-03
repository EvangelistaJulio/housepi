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
  expect_error(drop_group_rates(obj), regexp = err)
  expect_error(drop_contact_rates(obj), regexp = err)
  expect_error(drop_infectivity_profile(obj), regexp = err)
  expect_error(drop_all(obj), regexp = err)
  expect_error(add_parameter(obj), regexp = err)
  expect_error(add_augmented(obj), regexp = err)
  expect_error(add_infec_profile(obj, "mk22"), regexp = err)
  expect_error(evaluate_log_lik(obj), regexp = err)
  expect_error(run_inference(obj), regexp = err)
  expect_error(run_simulation(obj), regexp = err)
  expect_error(add_parameter_list(obj), regexp = err)
  expect_error(add_augmented_list(obj), regexp = err)
  expect_error(add_parameter_template(obj), regexp = err)
  expect_error(add_augmented_template(obj), regexp = err)
  expect_error(add_group_rate_list(obj), regexp = err)
  expect_error(add_infec_profile_list(obj), regexp = err)
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

  # no templates
  expect_warning(add_parameter_template(hhmodel),
                 regexp = "No template parameters")
  expect_warning(add_augmented_template(hhmodel),
                 regexp = "No template augmented")

  # parameters
  result <- add_parameter(hhmodel, name = "alpha", value = 0.5)
  pars <- get_parameters(result)
  expect_true(is.data.frame(pars))
  expect_true("alpha" %in% pars$name)

  result <- add_parameter_list(hhmodel, alpha = 1, beta = 2)
  pars <- get_parameters(result)
  expect_true(all(c("alpha", "beta") %in% pars$name))

  # infectivity profile
  expect_null(get_infec_profile(hhmodel))

  # augmented
  result <- add_augmented(hhmodel, name = "dummy", value = 3)
  aug <- get_augmented(result)
  expect_true(is.data.frame(aug))
  expect_true("dummy" %in% aug$name)

  result <- add_augmented_list(hhmodel,
                                dummy = list(model = "uniform"))
  aug <- get_augmented(result)
  expect_true("dummy" %in% aug$name)
})

############################################################################
## Drop wrappers
############################################################################

test_that("drop wrappers return clean models", {
  hhmodel <- new_housepi_model(model = "blank")
  hhmodel <- add_parameter(hhmodel, name = "alpha", value = 0.5)
  hhmodel <- add_augmented(hhmodel, name = "dummy", value = 3)

  # drop_parameters returns new object without parameters
  dropped <- drop_parameters(hhmodel)
  expect_equal(nrow(get_parameters(dropped)), 0)
  # original is unchanged

  expect_true("alpha" %in% get_parameters(hhmodel)$name)

  # drop_augmented
  dropped <- drop_augmented(hhmodel)
  expect_equal(nrow(get_augmented(dropped)), 0)
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
## Inference/simulation wrappers (chain_binomial)
############################################################################

test_that("evaluate_log_lik wrapper returns numeric", {
  hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial")
  hhmodel$add_parameter_template()
  ll <- evaluate_log_lik(hhmodel)
  expect_true(is.numeric(ll))
  expect_length(ll, 1)
})

test_that("run_inference wrapper returns HousepiPosteriorDraws", {
  hhmodel <- new_housepi_model(flu_2014_matsumoto, "chain_binomial")
  hhmodel$add_parameter_template()
  result <- run_inference(hhmodel, steps = 5, progress = FALSE)
  expect_s3_class(result, "HousepiPosteriorDraws")
})

test_that("run_simulation wrapper returns data.frame", {
  hhmodel <- new_housepi_model(flu_2014_matsumoto, "chain_binomial")
  hhmodel$add_parameter_template()
  result <- run_simulation(hhmodel)
  expect_true(is.data.frame(result))
  expect_true(all(c("id", "house_id", "status") %in% colnames(result)))
})

test_that("Group rate wrappers", {
  hhmodel <- new_housepi_model(flu_2014_matsumoto, "chain_binomial")
  hhmodel <- add_group_rate(hhmodel, "community", by = "member")
  groups_1 <- get_group_rates(hhmodel)
  expect_true(is.data.frame(groups_1))
  expect_true(nrow(groups_1) > 0)

  result <- add_group_rate_list(hhmodel, community = "member")
  groups_2 <- get_group_rates(result)
  expect_identical(groups_1, groups_2)

  # load custom definitions
  hhmodel <- hhmodel |>
    add_group_rate_list(
      susceptibility = list(
        strata = data.frame(member = c("student", NA),
                            param_name = c("sus_student", ".reference")),
        params = list(sus_student = 1.23)
      )
    )
  expect_equal(subset(hhmodel$parameters, name == "sus_student")$value, 1.23)

  # check for non-existing custom features
  hhmodel <- drop_group_rates(hhmodel)
  expect_warning(add_group_rate(
    hhmodel, "susceptibility",
    strata = data.frame(member = c("kaiju", NA),
                        param_name = c("sus_student", ".reference"))
  ), regexp = "not found")

  hhmodel$group_rates
})

test_that("add_infec_profile and add_infec_profile_list work", {
  hhmodel <- new_housepi_model(sarscov2_2021_icofs,
                               "augmented_infection_time")
  result <- add_infec_profile(hhmodel, func = "lognormal")
  ip <- get_infec_profile(result)
  expect_false(is.null(ip))
  expect_true("lognormal" %in% ip$func)

  # add_infec_profile_list
  result2 <- add_infec_profile_list(hhmodel,
    ip1 = list(func = "lognormal")
  )
  ip2 <- get_infec_profile(result2)
  expect_false(is.null(ip2))
  expect_true("lognormal" %in% ip2$func)
})
