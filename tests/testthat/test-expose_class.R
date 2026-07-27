library(testthat)

test_that("Reject if calling on non-HousepiModel objects", {
  obj <- "not a model"
  err <- "must be of class"
  add_data_long(obj, data.frame()) |> expect_error(err)
  add_data_wide(obj, data.frame()) |> expect_error(err)
  get_data_wide(obj) |> expect_error(err)
  get_data_long(obj) |> expect_error(err)
  get_parameters(obj) |> expect_error(err)
  get_augmented(obj) |> expect_error(err)
  get_group_risks(obj) |> expect_error(err)
  get_infectivity_profile(obj) |> expect_error(err)
  add_contact_risk(obj) |> expect_error(err)
  add_parameter(obj) |> expect_error(err)
  add_augmented(obj) |> expect_error(err)
  add_infectivity_profile(obj, "mk22") |> expect_error(err)
  evaluate_log_lik(obj) |> expect_error(err)
  run_inference(obj) |> expect_error(err)
  run_simulation(obj) |> expect_error(err)
  add_parameter_template(obj) |> expect_error(err)
  add_augmented_template(obj) |> expect_error(err)
  add_group_risk(obj) |> expect_error(err)
  run_custom_predictive_check(obj, NULL) |> expect_error(err)
  run_epi_predictive_check(obj, NULL) |> expect_error(err)
})

test_that("data_wide provided at construction is exposed via get_data_wide", {
  result <- new_housepi_model(flu_1918_baltimore, "test")
  dw <- get_data_wide(result)
  expect_true(is.data.frame(dw))
  expect_true("house_id" %in% colnames(dw))
  expect_true("status" %in% colnames(dw))

  # id automatically added
  expect_identical(dw$id, as.character(seq_len(dim(dw)[1])))
  # house_id cast to character
  expect_identical(dw$house_id, as.character(flu_1918_baltimore$house_id))
  # statuses are integer
  expect_identical(dw$status, flu_1918_baltimore$status)
})

test_that("add_data_long and get_data_long work", {
  hhmodel <- new_housepi_model(model = "test")
  # 'test' model has default data with id = "tommy", "gina"
  dl <- data.frame(id = hhmodel$data_wide$id,
                   measure = c(1.5, 2.5))
  result <- add_data_long(hhmodel, dl)
  out <- get_data_long(result)
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 2)
})

test_that("Wrappers work as expected on model 'test'", {
  hhmodel <- new_housepi_model(model = "test")

  # parameters
  result <- add_parameter(hhmodel, name = "alpha", value = 0.5)
  pars <- get_parameters(result)
  expect_true(is.data.frame(pars))
  expect_true("alpha" %in% pars$name)

  result <- add_parameter(hhmodel, alpha = 1, beta = 2)
  pars <- get_parameters(result)
  expect_true(all(c("alpha", "beta") %in% pars$name))

  # infectivity profile
  expect_null(get_infectivity_profile(hhmodel))

  # augmented
  result <- add_augmented(hhmodel, name = "dummy", value = 3)
  aug <- get_augmented(result)
  expect_true(is.data.frame(aug))
  expect_true("dummy" %in% aug$name)

  result <- hhmodel |>
    add_augmented(dummy = list(model = "uniform"))
  aug <- get_augmented(result)
  expect_true("dummy" %in% aug$name)
})

test_that("Inline weight labels do not collide with add_parameter formals", {
  hhmodel <- new_housepi_model(model = "test")

  # single-letter categorical weights c/d/p used to be partial-matched to the
  # const/discr/prior formals of add_parameter; they must reach the prior call
  # a collision would bind 'c' to 'const' and error ("Missing value for
  # constant parameter"); reaching a categorical 'snake' proves it did not
  result <- add_parameter(hhmodel, "snake", prior = "categorical",
                          a = 0.2, b = 0.1, c = 0.7)
  pars <- get_parameters(result)
  expect_true("snake" %in% pars$name)
  expect_equal(pars$prior[pars$name == "snake"], "categorical")
})

test_that("Bulk loading via `defs` mirrors inline name = value pairs", {
  hhmodel <- new_housepi_model(model = "test")
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

test_that("add_observation routes every input branch on `level`", {
  hhmodel <- data.frame(
    house_id = "outer_heaven", status = 1,
    onset_date = 3, test_date = 2
  ) |>
    new_housepi_model("test")

  # default level is member
  expect_identical(
    add_observation(hhmodel, "onset_date", model = "normal", mean = 0, sd = 1) |>
      get_observations() |>
      subset(name == "onset_date", "level", drop = TRUE),
    "member"
  )

  #> single definition
  single <- hhmodel |>
    add_observation("first_test", value = min(test_date@D), level = "house")
  #> multiple definitions
  bulk <- hhmodel |>
    add_observation(first_test = min(test_date@D), level = "house")
  #> import definitions
  input <- list(first_test = quote(min(test_date@D)))
  imported <- add_observation(hhmodel, import = input, level = "house")

  expect_identical(get_observations(bulk), get_observations(single))
  expect_identical(get_observations(imported), get_observations(single))
  expect_identical(
    get_observations(single)[["level"]], "house"
  )

  # an unknown level is rejected before anything is stored
  add_observation(hhmodel, big_boss = min(test_date@D), level = "snake") |>
    expect_error("should be one of")
})

############################################################################
## Wrappers with data
############################################################################

test_that("Group risk wrappers", {
  small_gr <- data.frame(
    house_id = c(1, 1, 2, 2, 2),
    status = c(1, 0, 1, 0, 0),
    member = c("student", "other", "student", "mother", "father")
  )
  hhmodel <- new_housepi_model(small_gr, "final_size")
  hhmodel <- add_group_risk(hhmodel, "community", by = "member")
  groups_1 <- get_group_risks(hhmodel)
  expect_true(is.data.frame(groups_1))
  expect_true(nrow(groups_1) > 0)

  expect_warning(
    result <- add_group_risk(hhmodel, community = "member"),
    "[Oo]verwriting.*community"
  )
  groups_2 <- get_group_risks(result)
  expect_identical(groups_1, groups_2)

  # load custom definitions
  hhmodel <- hhmodel |>
    add_group_risk(
      susceptibility = list(
        strata = data.frame(member = c("student", NA),
                            param_name = c("sus_student", ".reference")),
        params = list(sus_student = 1.23)
      )
    )
  expect_equal(subset(hhmodel$parameters, name == "sus_student")$value, 1.23)

  # check for non-existing custom features
  hhmodel |>
    add_group_risk(
      "susceptibility",
      strata = data.frame(member = c("kaiju", NA),
                          param_name = c("sus_student", ".reference")),
      append = FALSE
    ) |>
    expect_warning("[Uu]nknown.*member.*kaiju")

  # append = FALSE drops the existing group risks
  expect_true("susceptibility" %in% get_group_risks(hhmodel)[["group_risk"]])
  reset <- hhmodel |>
    add_group_risk(community = list(by = "member"), append = FALSE)
  expect_false("susceptibility" %in% get_group_risks(reset)[["group_risk"]])
  expect_setequal(get_group_risks(reset)[["group_risk"]], "community")
})

test_that("add_contact_risk wrapper input branches", {
  hhmodel <- data.frame(
    house_id = c(1, 1, 2, 2, 2),
    status = c(1, 0, 1, 0, 0),
    member = c("student", "other", "student", "mother", "father")
  ) |>
    new_housepi_model("final_size")

  strata <- data.frame(
    from_member = c("student", NA),
    to_member = c("student", NA),
    param_name = c("ctr_kept_you_waiting", ".reference")
  )

  # 'by' and 'strata' together is ambiguous
  hhmodel |>
    add_contact_risk(by = "member", strata = strata) |>
    expect_error("[Ii]mpossible")

  # neither 'by' nor 'strata' is empty
  hhmodel |>
    add_contact_risk() |>
    expect_error("[Ii]nvalid")

  # first call succeeds and returns a NEW model (copy() semantic)
  result <- hhmodel |>
    add_contact_risk(strata = strata,
                     params = list(ctr_kept_you_waiting = 1.21))
  expect_false(identical(hhmodel, result))
  # original is untouched (tidy contact-risk table stays empty)
  expect_equal(nrow(hhmodel$contact_risks), 0L)
  expect_true(nrow(result$contact_risks) > 0L)
  # rebuild_ ran: param is exposed on the returned model
  expect_true("ctr_kept_you_waiting" %in% get_parameters(result)$name)

  # second call warns about overwrite
  add_contact_risk(result, strata = strata,
                   params = list(ctr_kept_you_waiting = 1.21)) |>
    expect_warning("[Oo]verwriting")

  # getter
  cr <- get_contact_risks(result)
  expect_equal(cr, result$contact_risks)
  expect_true(is.data.frame(cr))
  expect_true(nrow(cr) > 0)
})

test_that("add_infectivity_profile single and bulk forms work", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained")
  ip <- hhmodel |>
    add_infectivity_profile(func = "lognormal") |>
    get_infectivity_profile()
  expect_false(is.null(ip))
  expect_true("lognormal" %in% ip$func)

  # bulk form
  ip2 <- hhmodel |>
    add_infectivity_profile(ip1 = list(func = "lognormal")) |>
    get_infectivity_profile()
  expect_false(is.null(ip2))
  expect_true("lognormal" %in% ip2$func)
})

test_that("add_infectivity_profile and list update the model directly", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |>
    add_infectivity_profile(func = "gamma", shape = 2, scale = 1)
  expect_true("gamma" %in% hhmodel$infectivity_profile$func)
  expect_warning(
    hhmodel <- add_infectivity_profile(hhmodel, func = "lognormal"),
    regexp = "Overwriting"
  )
  expect_true("lognormal" %in% hhmodel$infectivity_profile$func)
})
