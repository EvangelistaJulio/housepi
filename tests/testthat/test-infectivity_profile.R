library(testthat)

####  Helpers  ####

# model carrying no infectivity profile of its own
mk_infec <- local({
  hh_model <- sarscov2_2021_icofs |>
    new_housepi_model("case-ascertained") |>
    add_parameter(import = "template") |>
    add_latent(import = "template")
  list(model = hh_model)
})

# same model, but in discrete-time
mk_infec_discr <- local({
  hh_model <- sarscov2_2021_icofs |>
    new_housepi_model("case-ascertained", discr = TRUE, quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_latent(import = "template")
  list(model = hh_model)
})

# three-member household with fixed transmission, followed until day 20
mk_house <- local({
  hh_model <- data.frame(
    id = c("m1", "m2", "m3"), hh_id = "h1",
    status = c(1, 1, 0), start_date = 0, end_date = 20, datum = c(2, 9, 4)
  ) |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_parameter(alpha = 0.001, beta = 0.07, gamma = 0.9, hh_size_ref = 4)
  list(model = hh_model)
})

# same household in discrete time, followed until day 6
mk_house_discr <- local({
  hh_model <- data.frame(
    id = c("m1", "m2", "m3"), hh_id = "h1",
    status = c(1, 1, 0), start_date = 0, end_date = 6
  ) |>
    new_housepi_model("case-ascertained", discr = TRUE, quiet = TRUE) |>
    add_parameter(alpha = 0.001, beta = 0.07, gamma = 0.9, hh_size_ref = 4)
  list(model = hh_model)
})

####  Definitions  ####

test_that("Infectivity profile definition errors are recognized", {
  # function name of the wrong type
  mk_infec$model |>
    add_infectivity_profile(func = 2, discr = FALSE, a = 1) |>
    expect_error(".func. of infectivity profile must be a non-empty string")
  # the name of a listed definition is a label, not its function
  mk_infec$model |>
    add_infectivity_profile(import = list(gamma = list(shape = 2))) |>
    expect_error("[Mm]issing .func. of infectivity profile .gamma.")
  # arguments without `func` are listed definitions naming no function
  mk_infec$model |>
    add_infectivity_profile(shape = 2, scale = 1) |>
    expect_error("[Mm]issing .func. of infectivity profile .shape.")
  # unknown family
  mk_infec$model |>
    add_infectivity_profile(func = "unknown", x = 5) |>
    expect_error("[Uu]nknown infectivity profile")
  # a categorical draws codes, not days
  mk_infec_discr$model |>
    add_infectivity_profile(func = "categorical", a = 1, b = 2) |>
    expect_error("[Uu]nknown infectivity profile .categorical.")
  # a name outside the naming convention
  mk_infec$model |>
    add_infectivity_profile(func = "gam|ma", shape = 2) |>
    expect_error("[Ff]unction .gam.ma. contained invalid")
  # mixing definition modes
  mk_infec$model |>
    add_infectivity_profile("gamma", import = "template") |>
    expect_error("[Mm]ixing")
  # importing from unknown source
  mk_infec$model |>
    add_infectivity_profile(import = "templte") |>
    expect_error("[Uu]nknown import source")
})

test_that("Infectivity profile definitions are correct", {
  defs <- mk_infec$model |>
    add_infectivity_profile(func = "gamma", shape = 1, scale = 2,
                            min = 2, max = 10) |>
    get_infectivity_profiles()

  defs$func |> expect_identical("gamma")
  defs$func_args[[1]] |> expect_equal(c(shape = 1, scale = 2))
  defs$status |> expect_equal(NA_real_)
  defs$discr |> expect_identical(FALSE)
  defs$min |> expect_equal(2)
  defs$max |> expect_equal(10)

  # a function taking no argument
  defs <- mk_infec$model |>
    add_infectivity_profile("uniform", max = 4) |>
    get_infectivity_profiles()
  defs$func |> expect_identical("uniform")
  defs$func_args[[1]] |> expect_all_na()
  defs$max |> expect_equal(4)

  # bulk form
  ## UPDATE: the bulk form names a definition only to discard the name: a
  ## profile has no name of its own, so `profile =` labels nothing
  mk_infec$model |>
    add_infectivity_profile(profile = list(func = "lognormal")) |>
    get_infectivity_profiles() |>
    getElement("func") |>
    expect_identical("lognormal")
})

test_that("Template import loads the profile shipped with the model", {
  defs <- mk_infec$model |>
    add_infectivity_profile(import = "template") |>
    get_infectivity_profiles()

  defs$func |> expect_identical("gamma")
  defs$func_args[[1]] |> expect_equal(c(mean = 4, sd = 1.7))
  defs$discr |> expect_identical(FALSE)
})

test_that("Status-specific definitions replace the generic one", {
  defs <- mk_infec$model |>
    add_infectivity_profile("ferretti20", discr = TRUE,
                            incub_period = 5, max = 20) |>
    add_infectivity_profile(
      import = list(
        list(func = "lognormal", status = 1),
        list(func = "ferretti20", incub_period = 5, status = 2),
        profile = "uniform"
      ),
      append = FALSE
    ) |>
    get_infectivity_profiles()

  defs$func |> expect_identical(c("uniform", "lognormal", "ferretti20"))
  defs$func_args[[3]] |>
    expect_equal(c(incub_period = 5, shape = 5.85, scale = 1.85,
                   c_loc = -4, mean_incub = 5.42))
  defs$status |> expect_equal(c(NA_real_, 1, 2))

  # one call over several statuses fills the defaults of every row it appends
  defs <- mk_infec$model |>
    add_infectivity_profile("gamma", shape = 2, status = c(1, 2)) |>
    get_infectivity_profiles()
  defs$func_args |> expect_equal(rep(list(c(shape = 2, scale = 1)), 2))
})

test_that("Spurious names of imported definitions are dropped", {
  external_list <- list(
    profile = list(func = "gamma", status = NA,
                   discr = TRUE, shape = ~ shape, scale = ~ scale),
    profile_2 = list(func = "gamma", status = 2,
                     discr = TRUE, shape = ~ shape, scale = ~ scale)
  )

  defs <- mk_infec$model |>
    add_parameter(shape = 4, scale = 2) |>
    add_infectivity_profile(import = external_list) |>
    get_infectivity_profiles()

  defs$status <- NULL
  unique(defs) |> expect_nrow(1L)
})

test_that("A new profile overwrites the previous one", {
  hh_model <- mk_infec$model |>
    add_infectivity_profile(func = "gamma", shape = 2, scale = 1)
  get_infectivity_profiles(hh_model)$func |> expect_identical("gamma")
  expect_warning(
    hh_model <- add_infectivity_profile(hh_model, func = "lognormal"),
    "Overwriting"
  )
  get_infectivity_profiles(hh_model)$func |> expect_identical("lognormal")
})

test_that("The domain minimum defaults to the first transmission day", {
  profile_min <- function(hh_model, ...) {
    hh_model |>
      add_infectivity_profile(...) |>
      get_infectivity_profiles() |>
      getElement("min")
  }

  # a discrete profile moves it even in continuous time
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

test_that("A discrete-time model refuses a profile declared continuous", {
  mk_infec_discr$model |>
    add_infectivity_profile("gamma", discr = FALSE, mean = 4, sd = sqrt(3)) |>
    expect_error("must be discrete")
})

test_that("A discrete-time model resolves a bare latent to the lattice", {
  # the model imposes discretization on its latents, no per-definition flag
  mk <- function(min) {
    mk_infec_discr$model |>
      add_latent(infec_time = list(min = !!min, max = 3), append = FALSE) |>
      add_infectivity_profile("uniform", max = 4) |>
      evaluate_log_lik()
  }
  mk(0) |> expect_all_finite()

  # the flat lattice law owns every integer value in the window, so a bound that
  # is not one has no reading, and the error names the latent that carries it
  mk(0.5) |> expect_error("min.*uniform.*discrete.*infec_time")
})

test_that("A profile dependency reached through 'min' is reported", {
  # 'min' carries the only reference to the missing column
  mk_infec_discr$model |>
    add_infectivity_profile("uniform", min = datum, max = infec_time) |>
    evaluate_log_lik() |>
    expect_error("datum")
})

####  End-to-end  ####

test_that("Declaration order does not matter", {
  # the profile folds into the latent build, so it reads parameters that
  # only reach the engine on their own decorator
  before <- sarscov2_2021_icofs |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_infectivity_profile(import = "template") |>
    add_latent(import = "template") |>
    add_parameter(import = "template")
  after <- add_infectivity_profile(mk_infec$model, import = "template")

  # both reach the wall: the log-likelihood is random, only finiteness is
  # comparable (see evaluate_log_lik)
  evaluate_log_lik(before) |> expect_all_finite()
  evaluate_log_lik(after) |> expect_all_finite()
})

test_that("A model needing a profile throws", {
  get_infectivity_profiles(mk_infec$model) |> expect_null()
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

test_that("A profile transmitting at generation time 0 is refused", {
  # distribution with non-zero at 0
  mk_infec$model |>
    add_infectivity_profile("gamma", mean = 4, sd = sqrt(3), loc = -1) |>
    evaluate_log_lik() |>
    expect_error("non-zero at generation")

  # discretization resolves it automatically
  mk_infec$model |>
    add_infectivity_profile("exponential", discr = TRUE) |>
    evaluate_log_lik() |>
    expect_all_finite()

  # a registered profile considered as authored
  # 'ferretti20' is a custom profile with non-zero pdf at origin
  mk_infec$model |>
    add_infectivity_profile("ferretti20", incub_period = 5) |>
    evaluate_log_lik() |>
    expect_all_finite()
})

test_that("A non-zero profile at day 0 is always refused in discrete time", {
  mk_infec_discr$model |>
    add_infectivity_profile("gamma", mean = 4, sd = sqrt(3), min = 0) |>
    evaluate_log_lik() |>
    expect_error("non-zero at generation day 0")
  mk_infec_discr$model |>
    add_infectivity_profile("ferretti20", incub_period = 5, min = 0) |>
    evaluate_log_lik() |>
    expect_error("non-zero at generation day 0")
  # simulation refuses it too, before any household is drawn
  mk_house_discr$model |>
    add_infectivity_profile("gamma", mean = 4, sd = sqrt(3), min = 0) |>
    run_simulation() |>
    expect_error("non-zero at generation day 0")
})

test_that("Shifting a profile with 'loc' matches the explicit bounds", {
  # the support is loc + [min, max]: both spellings place it on
  # [incub_period, incub_period + 2]. The follow-up stops at day 6, so a
  # profile reaching past it loses the mass it puts outside
  hh_model <- mk_house_discr$model |>
    # deterministic, so the likelihood is reproducible across builds
    add_latent(infec_time = start_date, incub_period = 3)

  shifted <- hh_model |>
    add_infectivity_profile("uniform", loc = incub_period, max = 2) |>
    evaluate_log_lik()
  hh_model |>
    add_infectivity_profile("uniform", min = incub_period,
                            max = incub_period + 2) |>
    evaluate_log_lik() |>
    expect_equal(shifted)

  # and 'max' still widens the support it is handed
  wider <- hh_model |>
    add_infectivity_profile("uniform", loc = incub_period, max = 4) |>
    evaluate_log_lik()
  isTRUE(all.equal(shifted, wider)) |> expect_false()
})

test_that("A profile shifted by a sampled latent is judged on its support", {
  # the incubation is drawn, so the shift holds no value while the guard runs:
  # the support it places starts at day 1 and clears the day of infection
  hh_model <- add_latent(mk_house_discr$model, infec_time = start_date)
  kept_away <- hh_model |>
    add_latent(
      "incub_period", model = "gamma", discr = TRUE, min = 1, mean = 3, sd = 2
    )
  kept_away |>
    add_infectivity_profile("uniform", loc = incub_period, max = 2) |>
    evaluate_log_lik() |>
    expect_all_finite()
  # the explicit spelling of the same support agrees
  kept_away |>
    add_infectivity_profile("uniform", min = incub_period,
                            max = incub_period + 2) |>
    evaluate_log_lik() |>
    expect_all_finite()

  # simulation judges the undrawn shift on the same support
  simulated <- kept_away |>
    add_infectivity_profile("uniform", loc = incub_period, max = 2) |>
    run_simulation()
  simulated$incub_period[simulated$status > 0] |> expect_all_finite()

  # dropping the incubation's own minimum makes day 0 attainable
  not_kept <- hh_model |>
    add_latent(
      "incub_period", model = "gamma", discr = TRUE, mean = 3, sd = 2
    ) |>
    add_infectivity_profile("uniform", loc = incub_period, max = 2)
  not_kept |>
    evaluate_log_lik() |>
    expect_error("not kept away from generation day 0")
  not_kept |>
    run_simulation() |>
    expect_error("not kept away from generation day 0")
})

test_that("A latent parameterizing the profile drives transmission", {
  hh_model <- mk_house$model |>
    add_latent("infec_time", min = start_date, max = end_date)

  # the profile is resolved per member, so the latent reaches the household
  # cache without ever passing through 'infec_time'
  short <- hh_model |>
    add_latent(incub_period = 3) |>
    add_infectivity_profile("ferretti20", incub_period = incub_period) |>
    evaluate_log_lik()
  long <- hh_model |>
    add_latent(incub_period = 8) |>
    add_infectivity_profile("ferretti20", incub_period = incub_period) |>
    evaluate_log_lik()
  identical(short, long) |> expect_false()
})

test_that("A profile reading a parameter follows each draw in simulation", {
  # the first case of each household infects the others one generation later
  hh_model <- data.frame(id = seq(150), hh_id = rep(seq(50), each = 3),
                         status = 1, start_date = 0, end_date = 40) |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_parameter(alpha = 1e-6, beta = 50, gamma = 0, hh_size_ref = 3,
                  rho = list(prior = "uniform", min = 1, max = 20)) |>
    add_latent("infec_time", min = start_date, max = end_date) |>
    add_infectivity_profile("gamma", mean = rho, sd = 0.5)
  # time from the first infection of each household to the other ones
  gaps <- function(simulated) {
    unlist(tapply(simulated$infec_time, simulated$hh_id, function(t) {
      t <- t[!is.na(t)]
      t[t > min(t)] - min(t)
    }))
  }

  # gaps below 6 for a mean of 2, above 6 for a mean of 12
  simulated <- run_simulation(hh_model, draws = data.frame(rho = c(2, 12)))
  gaps(simulated[[1]]) |> expect_all_lt(6)
  gaps(simulated[[2]]) |> expect_all_gt(6)
})

test_that("The household snapshot reports the profile of each member", {
  hh_model <- mk_house$model
  profiles_of <- \(h) inspect_state(h, "h1")$infectivity_profiles

  # the model requires a profile, so the entry is kept but reads empty
  profiles_of(hh_model) |> expect_nrow(0L)

  # a generic profile repeats on every infected member, susceptibles read NA
  generic <- hh_model |>
    add_latent(infec_time = start_date) |>
    add_infectivity_profile("gamma", shape = 2, scale = 1) |>
    profiles_of()

  names(generic) |>
    expect_identical(c("id", "discr", "func", "min", "max", "func_args"))
  generic$id |> expect_identical(c("m1", "m2", "m3"))
  generic$func |> expect_identical(c("gamma", "gamma", NA_character_))
  generic$discr |> expect_identical(c(FALSE, FALSE, NA))
  generic$func_args[[1]] |> expect_equal(c(shape = 2, scale = 1))
  generic$func_args[[2]] |> expect_identical(generic$func_args[[1]])
  # the susceptible carries no profile at all
  generic$func_args[[3]] |> expect_length(0)

  # a latent profile is resolved per member: same house, two incubations
  latent <- hh_model |>
    add_latent(infec_time = start_date, incub_period = datum) |>
    add_infectivity_profile("ferretti20", incub_period = incub_period)

  # the latent storage is not allocated until the model is evaluated
  profiles_of(latent)$func |> expect_identical(rep(NA_character_, 3))

  invisible(evaluate_log_lik(latent))
  snapshot <- profiles_of(latent)
  snapshot$func |>
    expect_identical(c("ferretti20", "ferretti20", NA_character_))
  unname(snapshot$func_args[[1]]["incub_period"]) |> expect_equal(2)
  unname(snapshot$func_args[[2]]["incub_period"]) |> expect_equal(9)
})
