library(testthat)

# model carrying no infectivity profile of its own
mk_infec <- local({
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |>
    add_parameter(import = "template") |>
    add_augmented(import = "template")
  list(model = hhmodel)
})

# same model, but in discrete-time
mk_infec_discr <- local({
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained", discr = TRUE, quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_augmented(import = "template")
  list(model = hhmodel)
})

############################################################################
################ INFECTIVITY PROFILE : DEFINITIONS #########################
############################################################################

test_that("error handling in infectivity profile definitions", {
  # missing function name
  mk_infec$model |>
    add_infectivity_profile(func = NA, discr = FALSE, a = 1) |>
    expect_error("[Ii]ncomplete definition")
  # unknown family
  mk_infec$model |>
    add_infectivity_profile(func = "fulton", x = 5) |>
    expect_error("[Uu]nknown infectivity profile")
  # mixing definition modes
  mk_infec$model |>
    add_infectivity_profile("sokolov", import = "template") |>
    expect_error("[Mm]ixing")
  # importing from unknown source
  mk_infec$model |>
    add_infectivity_profile(import = "templte") |>
    expect_error("[Uu]nknown import source")
})

test_that("infectivity profile definitions are correct", {
  defs <- mk_infec$model |>
    add_infectivity_profile(func = "gamma", shape = 1, scale = 2,
                            min = 2, max = 10) |>
    get_infectivity_profiles()

  expect_equal(defs$func, "gamma")
  expect_equal(defs$func_args[[1]], c(shape = 1, scale = 2))
  expect_equal(defs$status, NA_real_)
  expect_equal(defs$discr, FALSE)
  expect_equal(defs$min, 2)
  expect_equal(defs$max, 10)
})

test_that("infectivity profile with discrete passes on continuos time", {
  # add a discrete profile
  defs <- mk_infec$model |>
    add_infectivity_profile("ferretti20", discr = TRUE,
                            incub_period = 5, max = 20) |>
    get_infectivity_profiles()

  expect_equal(defs$func, "ferretti20")
  expect_equal(defs$func_args[[1]], c(incub_period = 5))
  expect_equal(defs$discr, TRUE)
  expect_equal(defs$max, 20)
})

test_that("status specific definition work fine", {
  defs <- mk_infec$model |>
    add_infectivity_profile("ferretti20", discr = TRUE,
                            incub_period = 5, max = 20) |>
    add_infectivity_profile(
      import = list(
        list(func = "lognormal", status = 1),
        list(func = "ferretti20", incub_period = 5, status = 2),
        v1 = "uniform"
      ),
      append = FALSE
    ) |>
    get_infectivity_profiles()

  expect_equal(defs$func, c("lognormal", "ferretti20", "uniform"))
  expect_equal(defs$func_args[[2]], c(incub_period = 5))
  expect_equal(defs$status, c(1, 2, NA_real_))
})

test_that("spurious names are dropped", {
  external_list <- list(
    ip1 = list(func = "gamma", status = NA,
               discr = TRUE, shape = ~ shape@P, scale = ~ scale@P),
    ip2 = list(func = "gamma", status = 2,
               discr = TRUE, shape = ~ shape@P, scale = ~ scale@P)
  )

  infec <- mk_infec$model |>
    add_parameter(shape = 4, scale = 2) |>
    add_infectivity_profile(import = external_list) |>
    get_infectivity_profiles()

  infec$status <- NULL
  expect_equal(nrow(unique(infec)), 1)
})

test_that("range minimum defaults to the first transmission day", {
  profile_min <- function(hhmodel, ...) {
    hhmodel |>
      add_infectivity_profile(...) |>
      get_infectivity_profiles() |>
      getElement("min")
  }

  # discrete profile change it even in continuous time
  profile_min(mk_infec$model, "gamma", mean = 4, sd = sqrt(3), discr = TRUE) |>
    expect_equal(1)
  # a user-supplied minimum wins in either build
  profile_min(mk_infec_discr$model, "gamma", mean = 4, sd = sqrt(3), min = 2) |>
    expect_equal(2)
  # template definitions return minimum according to model time choice
  profile_min(mk_infec$model, import = "template") |>
    expect_equal(0)
  profile_min(mk_infec_discr$model, import = "template") |>
    expect_equal(1)

  # 'loc' already places the origin: no further shift is stacked on top of it
  profile_min(mk_infec_discr$model, "gamma", mean = 4, sd = sqrt(3), loc = 2) |>
    expect_equal(0)
})

test_that("a discrete-time model refuses a profile declared continuous", {
  mk_infec_discr$model |>
    add_infectivity_profile("gamma", discr = FALSE, mean = 4, sd = sqrt(3)) |>
    expect_error("must be discrete")

  # the same declaration is honoured in continuous time
  mk_infec$model |>
    add_infectivity_profile("gamma", discr = FALSE, mean = 4, sd = sqrt(3)) |>
    get_infectivity_profiles() |>
    getElement("discr") |>
    expect_false()
})

test_that("a profile dependency reached through 'min' is reported", {
  # 'min' carries the only reference to the missing column
  mk_infec_discr$model |>
    add_infectivity_profile(
      "uniform", min = liquid_snake@D, max = infec_time
    ) |>
    evaluate_log_lik() |>
    expect_error("liquid_snake")
})

############################################################################
################ INFECTIVITY PROFILE : END-TO-END ##########################
############################################################################

test_that("declaration order does not matter", {
  # the profile folds into the augmented build, so it reads parameters that
  # only reach the engine on their own decorator
  before <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_infectivity_profile(import = "template") |>
    add_augmented(import = "template") |>
    add_parameter(import = "template")
  expect_silent(before)

  after <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_infectivity_profile(import = "template") |>
    add_augmented(import = "template")

  # both reach the wall: the log-likelihood is random, only finiteness is
  # comparable (see evaluate_log_lik)
  expect_true(is.finite(evaluate_log_lik(before)))
  expect_true(is.finite(evaluate_log_lik(after)))
})

test_that("a model needing a profile throws", {
  mk_infec$model |>
    get_infectivity_profiles() |>
    expect_null()
  mk_infec$model |>
    evaluate_log_lik() |>
    expect_error("[Mm]issing infectivity profile")
  mk_infec$model |>
    run_inference() |>
    expect_error("[Mm]issing infectivity profile")
  mk_infec$model |>
    run_simulation() |>
    expect_error("[Mm]issing infectivity profile")
})

test_that("a profile transmitting at generation time 0 is refused", {
  # distribution with non-zero at 0
  mk_infec$model |>
    add_infectivity_profile("gamma", mean = 4, sd = sqrt(3), loc = -1) |>
    evaluate_log_lik() |>
    expect_error("non-zero at generation")

  # FIXME: disabled assert in simulation because of clash with status-specific
  # same behavior must occur in simulation
  # mk_infec$model |>
  #   add_infectivity_profile("exponential") |>
  #   run_simulation() |>
  #   expect_error("non-zero at generation")

  # discretization resolves it automatically
  mk_infec$model |>
    add_infectivity_profile("exponential", discr = TRUE) |>
    evaluate_log_lik() |>
    expect_silent()

  # a registered profile considered as authored
  # NOTE: 'ferretti20' is a custom with non-zero pdf at origin
  mk_infec$model |>
    add_infectivity_profile("ferretti20", incub_period = 5) |>
    evaluate_log_lik() |>
    expect_silent()
})

test_that("non-zero profile for day-0 is always refused in discrete time", {
  mk_infec_discr$model |>
    add_infectivity_profile("gamma", mean = 4, sd = sqrt(3), min = 0) |>
    evaluate_log_lik() |>
    expect_error("non-zero at generation day 0")
  mk_infec_discr$model |>
    add_infectivity_profile("ferretti20", incub_period = 5, min = 0) |>
    evaluate_log_lik() |>
    expect_error("non-zero at generation day 0")
})

test_that("shifting a profile with 'loc' matches the explicit bounds", {
  # the support is loc + [min, max]: both spellings place it on
  # [.incub_period, .incub_period + 2]. The follow-up stops at day 6, so a
  # profile reaching past it loses the mass it puts outside
  outer_heaven <- data.frame(
    id = c("venom", "kaz", "ocelot"), house_id = "outer_heaven",
    status = c(1, 1, 0), start_date = 0, end_date = 6
  ) |>
    new_housepi_model("case_ascertained", discr = TRUE, quiet = TRUE) |>
    add_parameter(alpha = 0.001, beta = 0.07, gamma = 0.9, hhsize_ref = 4) |>
    # deterministic, so the likelihood is reproducible across builds
    add_augmented(infec_time = start_date@D, .incub_period = 3)

  shifted <- outer_heaven |>
    add_infectivity_profile("uniform", loc = .incub_period, max = 2) |>
    evaluate_log_lik()
  explicit <- outer_heaven |>
    add_infectivity_profile("uniform", min = .incub_period,
                            max = .incub_period + 2) |>
    evaluate_log_lik()
  expect_equal(shifted, explicit)

  # and 'max' still widens the support it is handed
  wider <- outer_heaven |>
    add_infectivity_profile("uniform", loc = .incub_period, max = 4) |>
    evaluate_log_lik()
  expect_false(isTRUE(all.equal(shifted, wider)))
})

test_that("an augmented parameterizing the profile drives transmission", {
  shadow_moses <- data.frame(
    id = c("ocelot", "raiden", "otacon"), house_id = "shadow_moses",
    status = c(1, 1, 0), start_date = 0, end_date = 20
  ) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0.001, beta = 0.07, gamma = 0.9, hhsize_ref = 4) |>
    add_augmented(
      infec_time = list(model = "uniform", min = start_date@D, max = end_date@D)
    )

  # the profile is resolved per member, so the augmented reaches the household
  # cache without ever passing through 'infec_time'
  short <- shadow_moses |>
    add_augmented(.incub_period = 3) |>
    add_infectivity_profile("ferretti20", incub_period = .incub_period) |>
    evaluate_log_lik()
  long <- shadow_moses |>
    add_augmented(.incub_period = 8) |>
    add_infectivity_profile("ferretti20", incub_period = .incub_period) |>
    evaluate_log_lik()
  expect_false(identical(short, long))
})

test_that("the household snapshot reports the profile of each member", {
  zanzibar <- data.frame(
    id = c("venom", "kaz", "ocelot"), house_id = "zanzibar",
    status = c(1, 1, 0), start_date = 0, end_date = 20, latency = c(2, 9, 4)
  ) |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_parameter(alpha = 0.001, beta = 0.07, gamma = 0.9, hhsize_ref = 4)

  # the model requires a profile, so the entry is kept but reads empty
  zanzibar$export_house_snapshot("zanzibar")$infectivity_profiles |>
    nrow() |>
    expect_equal(0L)

  # a generic profile repeats on every infected member, susceptibles read NA
  generic <- zanzibar |>
    add_augmented(infec_time = start_date@D) |>
    add_infectivity_profile("gamma", shape = 2, scale = 1) |>
    (\(x) x$export_house_snapshot("zanzibar")$infectivity_profiles)()

  expect_equal(names(generic),
               c("id", "discr", "func", "min", "max", "func_args"))
  expect_equal(generic$id, c("venom", "kaz", "ocelot"))
  expect_equal(generic$func, c("gamma", "gamma", NA_character_))
  expect_equal(generic$discr, c(FALSE, FALSE, NA))
  expect_equal(generic$func_args[[1]], c(shape = 2, scale = 1))
  expect_equal(generic$func_args[[2]], generic$func_args[[1]])
  # the susceptible carries no profile at all
  expect_length(generic$func_args[[3]], 0)

  # an augmented profile is resolved per member: same house, two incubations
  augmented <- zanzibar |>
    add_augmented(infec_time = start_date@D, .incub_period = latency@D) |>
    add_infectivity_profile("ferretti20", incub_period = .incub_period)

  # the augmented storage is not allocated until the model is evaluated
  augmented$export_house_snapshot("zanzibar")$infectivity_profiles$func |>
    expect_equal(rep(NA_character_, 3))

  invisible(evaluate_log_lik(augmented))
  snapshot <- augmented$export_house_snapshot("zanzibar")$infectivity_profiles

  expect_equal(snapshot$func, c("ferretti20", "ferretti20", NA_character_))
  expect_equal(unname(snapshot$func_args[[1]]["incub_period"]), 2)
  expect_equal(unname(snapshot$func_args[[2]]["incub_period"]), 9)
})

test_that("add_infectivity_profile single and bulk forms work", {
  ip <- mk_infec$model |>
    add_infectivity_profile(func = "lognormal") |>
    get_infectivity_profiles()
  expect_true("lognormal" %in% ip$func)

  # bulk form
  ip2 <- mk_infec$model |>
    add_infectivity_profile(ip1 = list(func = "lognormal")) |>
    get_infectivity_profiles()
  expect_true("lognormal" %in% ip2$func)
})

test_that("add_infectivity_profile and list update the model directly", {
  hhmodel <- mk_infec$model |>
    add_infectivity_profile(func = "gamma", shape = 2, scale = 1)
  expect_true("gamma" %in% hhmodel$infectivity_profiles$func)
  expect_warning(
    hhmodel <- add_infectivity_profile(hhmodel, func = "lognormal"),
    regexp = "Overwriting"
  )
  expect_true("lognormal" %in% hhmodel$infectivity_profiles$func)
})
