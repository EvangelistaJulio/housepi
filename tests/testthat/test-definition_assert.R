library(testthat)

test_that("Helper format errors keeps full multi-line text of the pick", {
  # single error: returns the full (multi-line) text for the cli split
  helper_format_error(list("first line\nlast line")) |>
    expect_identical("first line\nlast line")

  # multiple errors: "Unused" or "Missing" are prioritized
  err <- list(a = "some error", b = "Unused argument x")
  helper_format_error(err) |> expect_identical("Unused argument x")
  err <- list(a = "some error", b = "Missing argument x")
  helper_format_error(err) |> expect_identical("Missing argument x")

  # "func" value in err triggers func branch (checks value, not key)
  err <- list(dist = "dist error", func = "func")
  helper_format_error(err) |> expect_identical(list(func = "func"))

  # fallback to first
  err <- list(a = "first", b = "second")
  helper_format_error(err) |> expect_identical("first")
})

test_that("Assert valid definition name rejects invalid names", {
  # empty string
  assert_valid_definition_name("", "parameter") |>
    expect_error("non-empty")
  # NA
  assert_valid_definition_name(NA, "parameter") |>
    expect_error("non-empty")
  # non-character
  assert_valid_definition_name(42, "parameter") |>
    expect_error("non-empty")
  # special characters
  assert_valid_definition_name("a b", "parameter") |>
    expect_error("invalid characters")
  assert_valid_definition_name("a@b", "parameter") |>
    expect_error("invalid characters")
  # reserved name pi
  assert_valid_definition_name("pi", "parameter") |>
    expect_error("reserved")
  # a leading dot is reserved to what housepi generates itself
  assert_valid_definition_name(".ration", "observation") |>
    expect_error("cannot start with a dot")
  assert_valid_definition_name(".ration", "latent") |>
    expect_error("cannot start with a dot")
  assert_valid_definition_name(".ration", "status process",
                               allow_reserved = TRUE) |>
    expect_null()
  # valid names pass
  assert_valid_definition_name("alpha", "parameter") |> expect_null()
  assert_valid_definition_name("beta_1", "parameter") |> expect_null()
})

test_that("Assert valid status rejects non-positive integers", {
  # valid
  assert_valid_status(c(NA, 1, 2, 3)) |> expect_null()
  # zero
  assert_valid_status(c(0, 1)) |>
    expect_error("positive integer")
  # negative
  assert_valid_status(c(-1, 2)) |>
    expect_error("positive integer")
  # decimal
  assert_valid_status(c(1.5)) |>
    expect_error("positive integer")
  # above the engine status storage
  assert_valid_status(254) |> expect_null()
  assert_valid_status(c(1, 255)) |>
    expect_error("maximum of 254")
})

test_that("Assert valid data columns holds the definition alphabet", {
  # a column is a key a definition reads: same alphabet as a name
  assert_valid_data_columns(c("id", "3", "if", "dose.2", "vacc_status")) |>
    expect_null()
  assert_valid_data_columns(c("id", "bad name")) |>
    expect_error("bad name")
  # the grammar delimiter never reaches the engine as a key
  assert_valid_data_columns("bad|name") |>
    expect_error("bad\\|name")

  # every offender is named at once, and the argument is the caller's
  assert_valid_data_columns(c("bad name", "ok.name", "bad!")) |>
    expect_error("bad name")
})

test_that("Reserved column lookalikes are named once per target", {
  # every spelling folds onto the reserved column it stands for
  assert_valid_data_columns(c("ID", "hhid", "Status")) |>
    expect_warning("ID.*individual identifier") |>
    expect_warning("hhid.*household identifier") |>
    expect_warning("Status.*infection status")
  assert_valid_data_columns(c("id", "hh.id", "Household_ID", "status")) |>
    expect_warning("hh.id.*Household_ID.*rename it to.*hh_id")

  # bare house and household stand for the household
  assert_valid_data_columns(c("id", "House", "status")) |>
    expect_warning("House.*household identifier")

  # the reserved column is there: a lookalike is a covariate
  assert_valid_data_columns(c("id", "hh_id", "status", "house_id", "ID")) |>
    expect_null()
  # an unrelated name, state included, is left alone
  assert_valid_data_columns(c("id", "hh_id", "status", "state", "member")) |>
    expect_null()

  # long data only reads the individual identifier
  assert_valid_data_columns(c("Status", "hhid", "measure"), "id") |>
    expect_null()
  assert_valid_data_columns(c("Id", "measure"), "id") |>
    expect_warning("Id.*individual identifier")
})

test_that("Definitions carry the defaults the user left out", {
  # a bare uniform prior reports its unit interval, a given bound is kept
  param <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    get_parameters()
  param$min |> expect_equal(c(0, 0, 0, NA))
  param$max |> expect_equal(c(1, 1, 3, NA))

  hhdata <- data.frame(hh_id = c("h1", "h1"), status = c(1, 0),
                       record = c(3L, 1L), record_2 = c(0.5, 2.1))
  hh_model <- hhdata |>
    new_housepi_model("final-size") |>
    add_parameter(
      rho = list(value = 0.3, prior = "gamma", shape = 2),
      rho_2 = list(value = 0.3, prior = "gamma", mean = 2, sd = 1),
      rho_3 = list(value = 0, prior = "normal", sd = 2, min = -1),
      rho_4 = list(func = "log", x = rho)
    ) |>
    add_observation(
      record = list(model = "binomial", size = 5, prob = rho),
      record_2 = list(model = "gamma", shape = 2, scale = 3)
    )
  param <- get_parameters(hh_model, name = c("rho", "rho_2", "rho_3", "rho_4"))
  # the parametrization in use lends its defaults, the other one has none
  param$args[[1]] |> expect_equal(c(shape = 2, scale = 1))
  param$args[[2]] |> expect_equal(c(mean = 2, sd = 1))
  param$args[[3]] |> expect_equal(c(sd = 2, mean = 0))
  # the open bound of `log` at 0 is passed back as it is declared
  param$min |> expect_equal(c(0, 0, -1, 0))

  obs <- get_observations(hh_model)
  obs$args[[2]] |> expect_equal(c(shape = 2, scale = 3))
  # the binomial upper bound reads `size`, so it is left to track it
  obs$min |> expect_equal(c(0, 0))
  obs$max |> expect_equal(c(NA_real_, NA_real_))
})

test_that("Constant arguments and bounds are checked at definition", {
  hh_model <- data.frame(hh_id = "h1", status = 1) |>
    new_housepi_model("test")

  # the latent build waits for the run, the definition does not
  hh_model |>
    add_latent(hidden = list(model = "categorical", 0, 0)) |>
    expect_error("weights summing to non-strictly")
  hh_model |>
    add_latent(hidden = list(model = "uniform", min = 2, max = 1)) |>
    expect_error("[Ii]nconsistent domain")

  # a sourced argument has no value at definition
  augm <- hh_model |>
    add_parameter(rho = list(prior = "gamma", shape = 2)) |>
    add_latent(hidden = list(model = "normal", sd = rho)) |>
    get_latent()
  augm$args[[1]] |> expect_identical(c(sd = "rho|P", mean = "0"))
})
