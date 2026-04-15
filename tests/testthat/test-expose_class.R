library(testthat)

test_that("Rejext if pipe calls get non-HousepiModel objects", {
  obj <- "not a model"
  err <- "class 'HousepiModel'"
  expect_error(add_data_wide(obj, data.frame()), regexp = err)
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
  expect_error(drop_infectivity_profile(obj), regexp = err)
  expect_error(drop_all(obj), regexp = err)
  expect_error(add_parameter(obj), regexp = err)
  expect_error(add_augmented(obj), regexp = err)
  expect_error(add_infec_profile(obj, "mk22"), regexp = err)
  expect_error(evaluate_log_lik(obj), regexp = err)
  expect_error(run_inference(obj), regexp = err)
  expect_error(run_simulation(obj), regexp = err)
  expect_error(load_parameter_list(obj), regexp = err)
  expect_error(load_augmented_list(obj), regexp = err)
  expect_error(load_parameter_template(obj), regexp = err)
  expect_error(load_augmented_template(obj), regexp = err)
  expect_error(load_group_rate_list(obj), regexp = err)
  expect_error(load_infec_profile_list(obj), regexp = err)
  expect_error(predicted_measures(obj, NULL), regexp = err)
  expect_error(predicted_cases_by_house_size(obj, NULL), regexp = err)
})

test_that("add_data_wide and get_data_wide work", {
  original <- helper_test_get_data("flu_1918_baltimora")
  hhmodel <- new_housepi_model("blank")
  result <- add_data_wide(hhmodel, original)
  dw <- get_data_wide(result)
  expect_true(is.data.frame(dw))
  expect_true("house_id" %in% colnames(dw))
  expect_true("status" %in% colnames(dw))

  # id automatically added
  expect_identical(dw$id, as.character(seq_len(dim(dw)[1])))
  # house_id cast to character
  expect_identical(dw$house_id, as.character(original$house_id))
  # statuses are integer
  expect_identical(dw$status, original$status)
})

test_that("add_data_long and get_data_long work", {
  hhmodel <- new_housepi_model("blank")
  # blank model has default data with id = "tommy", "gina"
  dl <- data.frame(id = c("tommy", "gina"),
                   measure = c(1.5, 2.5))
  result <- add_data_long(hhmodel, dl)
  out <- get_data_long(result)
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 2)
})

test_that("Wrappers work as expected on blank model", {
  hhmodel <- new_housepi_model("blank")

  #parameters
  result <- add_parameter(hhmodel, name = "alpha", value = 0.5)
  pars <- get_parameters(result)
  expect_true(is.data.frame(pars))
  expect_true("alpha" %in% pars$name)
  
  result <- load_parameter_list(hhmodel, alpha = 1, beta = 2)
  pars <- get_parameters(result)
  expect_true(all(c("alpha", "beta") %in% pars$name))
  
  #infectivity profile
  expect_null(get_infec_profile(hhmodel))

  # augmented
  result <- add_augmented(hhmodel, name = "dummy", value = 3)
  aug <- get_augmented(result)
  expect_true(is.data.frame(aug))
  expect_true("dummy" %in% aug$name)
  
  result <- load_augmented_list(hhmodel,
                                dummy = list(model = "uniform"))
  aug <- get_augmented(result)
  expect_true("dummy" %in% aug$name)
})


############################################################################
## Template loaders on blank model warn
############################################################################

test_that("load_parameter_template warns on blank model", {
  hhmodel <- new_housepi_model("blank")
  expect_warning(load_parameter_template(hhmodel),
                 regexp = "No template parameters")
})

test_that("load_augmented_template warns on blank model", {
  hhmodel <- new_housepi_model("blank")
  expect_warning(load_augmented_template(hhmodel),
                 regexp = "No template augmented")
})

############################################################################
## Drop wrappers
############################################################################

test_that("drop wrappers return clean models", {
  hhmodel <- new_housepi_model("blank")
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

  # drop_infectivity_profile
  expect_no_error(drop_infectivity_profile(hhmodel))

  # drop_all
  dropped <- drop_all(hhmodel)
  expect_equal(nrow(get_parameters(dropped)), 0)
  expect_equal(nrow(get_augmented(dropped)), 0)
})

############################################################################
## Summarise wrappers
############################################################################

test_that("summarise_cases and summarise_house_sizes work", {
  hhmodel <- new_housepi_model("chain_binomial")
  hhmodel$add_data_wide(helper_test_get_data("flu_1918_frederick"))
  expect_true(is.data.frame(summarise_cases(hhmodel)))
  expect_true(is.data.frame(summarise_house_sizes(hhmodel)))
})

############################################################################
## Inference/simulation wrappers (chain_binomial)
############################################################################

test_that("evaluate_log_lik wrapper returns numeric", {
  hhmodel <- new_housepi_model("chain_binomial")
  hhmodel$add_data_wide(helper_test_get_data("flu_1918_frederick"))
  hhmodel$load_parameter_template()
  ll <- evaluate_log_lik(hhmodel)
  expect_true(is.numeric(ll))
  expect_length(ll, 1)
})

test_that("run_inference wrapper returns HousepiPosteriorDraws", {
  hhmodel <- new_housepi_model("chain_binomial")
  hhmodel$add_data_wide(helper_test_get_data("flu_2014_matsumoto"))
  hhmodel$load_parameter_template()
  result <- run_inference(hhmodel, steps = 5, progress = FALSE)
  expect_s3_class(result, "HousepiPosteriorDraws")
})

test_that("run_simulation wrapper returns data.frame", {
  hhmodel <- new_housepi_model("chain_binomial")
  hhmodel$add_data_wide(helper_test_get_data("flu_2014_matsumoto"))
  hhmodel$load_parameter_template()
  result <- run_simulation(hhmodel)
  expect_true(is.data.frame(result))
  expect_true(all(c("id", "house_id", "status") %in% colnames(result)))
})

test_that("Group rate wrappers", {
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel$add_data_wide(helper_test_get_data("flu_2014_matsumoto"))
  hhmodel <- add_group_rate(hhmodel, "community", by = "member")
  groups_1 <- get_group_rates(hhmodel)
  expect_true(is.data.frame(groups_1))
  expect_true(nrow(groups_1) > 0)

  result <- load_group_rate_list(hhmodel, community = "member")
  groups_2 <- get_group_rates(result)
  expect_identical(groups_1, groups_2)
  
  # load custom definitions
  hhmodel <- hhmodel |>
    load_group_rate_list(
      susceptibility = list(
        strata = data.frame(member = c("student", NA),
                            param_name = c("sus_student", ".reference")),
        params = list(sus_student = 1.23)
      )
  )
  expect_equal(subset(hhmodel$parameters, name == "sus_student")$value, 1.23)
  
  # check for non-existing custom features
  hhmodel <- drop_group_rates(hhmodel)
  hhmodel <- add_group_rate(
    hhmodel, "susceptibility",
    strata = data.frame(member = c("kaiju", NA),
                        param_name = c("sus_student", ".reference"))
  )
  expect_warning(hhmodel$group_rates, regexp = "not found")
})

test_that("add_infec_profile and load_infec_profile_list work", {
  hhmodel <- new_housepi_model("augmented_infection_time")
  data_ait <- dplyr::mutate(helper_test_get_data("sarscov2_2021_icofs"))
  hhmodel$add_data_wide(data_ait)
  result <- add_infec_profile(hhmodel, func = "lognormal")
  ip <- get_infec_profile(result)
  expect_false(is.null(ip))
  expect_true("lognormal" %in% ip$func)

  # load_infec_profile_list
  result2 <- load_infec_profile_list(hhmodel,
    ip1 = list(func = "lognormal")
  )
  ip2 <- get_infec_profile(result2)
  expect_false(is.null(ip2))
  expect_true("lognormal" %in% ip2$func)
})
