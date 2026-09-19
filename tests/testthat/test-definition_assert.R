library(testthat)

test_that("Helper format errors keeps full multi-line text of the pick", {
  # single error: returns the full (multi-line) text for the cli split
  expect_equal(helper_format_error(list("first line\nlast line")),
               "first line\nlast line")

  # multiple errors: "Unused" or "Missing" are prioritized
  err <- list(a = "some error", b = "Unused argument x")
  expect_equal(helper_format_error(err), "Unused argument x")
  err <- list(a = "some error", b = "Missing argument x")
  expect_equal(helper_format_error(err), "Missing argument x")

  # "func" value in err triggers func branch (checks value, not key)
  err <- list(dist = "dist error", func = "func")
  expect_equal(helper_format_error(err), list(func = "func"))

  # fallback to first
  err <- list(a = "first", b = "second")
  expect_equal(helper_format_error(err), "first")
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
  expect_no_error(assert_valid_definition_name(".ration", "status observation",
                                               allow_reserved = TRUE))
  # valid names pass
  expect_no_error(assert_valid_definition_name("alpha", "parameter"))
  expect_no_error(assert_valid_definition_name("beta_1", "parameter"))
})

test_that("Assert valid status rejects non-positive integers", {
  # valid
  expect_no_error(assert_valid_status(c(NA, 1, 2, 3)))
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
  expect_no_error(assert_valid_status(254))
  assert_valid_status(c(1, 255)) |>
    expect_error("maximum of 254")
})

test_that("Assert valid data columns holds the definition alphabet", {
  # a column is a key a definition sources with `@D`: same alphabet as a name
  expect_no_error(
    assert_valid_data_columns(c("id", "3", "if", "dose.2", "vacc_status"))
  )
  assert_valid_data_columns(c("id", "bad name")) |>
    expect_error("bad name")
  assert_valid_data_columns("bad-name") |>
    expect_error("Only letters, numbers")
  # the grammar delimiter never reaches the engine as a key
  assert_valid_data_columns("bad|name") |>
    expect_error("bad\\|name")

  # every offender is named at once, and the argument is the caller's
  assert_valid_data_columns(c("bad name", "ok.name", "bad!"), "data_long") |>
    expect_error("bad name")
  assert_valid_data_columns(c("bad name", "ok.name", "bad!"), "data_long") |>
    expect_error("data_long")
})

test_that("Definitions carry the defaults the user left out", {
  # a bare uniform prior reports its unit interval, a given bound is kept
  param <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter(import = "template") |>
    get_parameters()
  expect_equal(param$min, c(0, 0, 0, NA))
  expect_equal(param$max, c(1, 1, 3, NA))

  hhdata <- data.frame(house_id = c("h1", "h1"), status = c(1, 0),
                       record = c(3L, 1L), record_2 = c(0.5, 2.1))
  hhmodel <- hhdata |>
    new_housepi_model("final_size") |>
    add_parameter(
      rho = list(value = 0.3, prior = "gamma", shape = 2),
      rho_2 = list(value = 0.3, prior = "gamma", mean = 2, sd = 1),
      rho_3 = list(value = 0, prior = "normal", sd = 2, min = -1),
      rho_4 = list(func = "log", x = rho@P)
    ) |>
    add_observation(
      record = list(model = "binomial", size = 5, prob = rho@P),
      record_2 = list(model = "gamma", shape = 2, scale = 3)
    )
  param <- get_parameters(hhmodel, name = c("rho", "rho_2", "rho_3", "rho_4"))
  # the parametrization in use lends its defaults, the other one has none
  expect_equal(param$args[[1]], c(shape = 2, scale = 1))
  expect_equal(param$args[[2]], c(mean = 2, sd = 1))
  expect_equal(param$args[[3]], c(sd = 2, mean = 0))
  # the open bound of `log` at 0 is passed back as it is declared
  expect_equal(param$min, c(0, 0, -1, 0))

  obs <- get_observations(hhmodel)
  expect_equal(obs$args[[2]], c(shape = 2, scale = 3))
  # the binomial upper bound reads `size`, so it is left to track it
  expect_equal(obs$min, c(0, 0))
  expect_equal(obs$max, c(NA_real_, NA_real_))
})
