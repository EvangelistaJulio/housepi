library(testthat)

# data frame carrying one column of each storage type the engine encodes
mk_types <- local({
  hhdata <- data.frame(id = c("m1", "m2"),
                       hhid = c("h1", "h1"),
                       status = c(1, 0),
                       date_col = as.Date(c("2024-01-01", "2024-01-02")),
                       char_col = c("male", "female"),
                       logic_col = c(TRUE, FALSE),
                       int_col = c(3L, 5L),
                       whole_col = c(3, 5),
                       num_col = c(1.5, 2.5),
                       empty_int = c(NA_integer_, NA_integer_),
                       empty_col = c(NA_real_, NA_real_),
                       datum = c("a", NA),
                       datum_2 = c(TRUE, NA),
                       datum_3 = NA_character_,
                       stringsAsFactors = FALSE)
  list(data = hhdata)
})

# factorized data coerced to the numeric frame the engine hands back
as_cpp_frame <- function(factored) {
  cpp_data <- suppressWarnings(data.frame(lapply(factored$data, as.numeric)))
  colnames(cpp_data) <- colnames(factored$data)
  cpp_data
}

#########################
###  UTILS : HELPERS  ###
#########################

test_that("Overwrite is a no-op on the same path and raises on failure", {
  # no-op on identical paths
  path <- tempfile()
  file.create(path)
  overwrite(path, path) |> expect_invisible()

  # destination directory does not exist -> file.copy returns FALSE
  overwrite(path, file.path(tempfile(), "missing", "x")) |>
    expect_error("user permissions")
})

test_that("Rendered text keeps data braces literal through a later cli pass", {
  datum <- c("a{b}", "{c}")
  txt <- inline_cli_("Got {datum}.")
  cli::format_inline(txt) |> expect_identical("Got a{b} and {c}.")
  # the error names the brace-carrying option verbatim
  stop_h("Unknown entry.", find_closest_match("unrelated", datum)) |>
    expect_error("Available options.*a\\{b\\}")
})

test_that("Find closest match suggests close options", {
  # typo within distance 1
  result <- find_closest_match("alpa", c("alpha", "beta"))
  result |> expect_match("Did you mean", all = FALSE)
  result |> expect_match("alpha", all = FALSE)

  # min_dist = 0 disables suggestions
  find_closest_match("alpa", c("alpha"), min_dist = 0) |> expect_null()

  # no close match returns available options
  result <- find_closest_match("unknown", c("alpha", "beta"))
  result |> expect_match("Available option", all = FALSE)

  # empty options
  find_closest_match("x", character(0)) |> expect_identical(c())
})

test_that("Deep clone helper work as expected on default types", {
  # scalar passthrough
  helper_clone(42) |> expect_equal(42)
  helper_clone("hello") |> expect_identical("hello")

  # data.frame independence
  original <- data.frame(x = 1:3)
  cloned <- helper_clone(original)
  cloned$x[1] <- 99
  original$x[1] |> expect_equal(1)
  cloned$x[1] |> expect_equal(99)

  # list of data.frames recursively cloned
  frames <- list(a = data.frame(x = 1), b = data.frame(y = 2))
  cloned <- helper_clone(frames)
  cloned$a$x[1] <- 99
  frames$a$x[1] |> expect_equal(1)
  cloned$a$x[1] |> expect_equal(99)
})

test_that("Setup options are formatted as character", {
  # id_draws with multiple values
  result <- char_setup_input_(list(steps = 100, burnin = 10, thin = 1,
                                   chains = 1, id_draws = c("m1", "m2")))
  result |> expect_s3_class("data.frame")
  result$id_draws |> expect_identical("m1:m2")
  vapply(result, is.character, logical(1)) |> expect_all_true()
  # id_draws are dropped if NULL
  char_setup_input_(list(steps = 100L, id_draws = NULL)) |>
    getElement("id_draws") |>
    expect_identical("")
})

test_that("Char data input handles NULL and data.frame", {
  char_data_input_(NULL) |> expect_null()
  result <- char_data_input_(data.frame(x = 1:3, y = c(4.5, 5.5, 6.5)))
  result |> expect_s3_class("data.frame")
  vapply(result, is.character, logical(1)) |> expect_all_true()

  # column names are engine keys: they cross verbatim, never through make.names
  exotic <- data.frame(a = 1, check.names = FALSE)
  exotic[["3"]] <- 2
  exotic[["if"]] <- 3
  exotic[["dose.2"]] <- 4
  names(char_data_input_(exotic)) |> expect_identical(names(exotic))
})

test_that("Draws are thinned evenly, first and last kept", {
  helper_thin_even_(seq(10), 4) |> expect_identical(c(1L, 4L, 7L, 10L))
  hhdata <- data.frame(datum = seq(10))
  helper_thin_even_(hhdata, 4)$datum |> expect_identical(c(1L, 4L, 7L, 10L))
  # fewer draws than asked, or no target: every draw kept
  helper_thin_even_(seq(3), 4) |> expect_identical(seq(3))
  helper_thin_even_(hhdata, NA) |> expect_identical(hhdata)
})

test_that("A count is a whole positive number", {
  assert_count_(3, "Number of replicates") |> expect_identical(3L)
  assert_count_(NA, "Number of replicates", allow_na = TRUE) |>
    expect_identical(NA_integer_)
  assert_count_(NA, "Number of replicates") |>
    expect_error("[Nn]umber of replicates should be a positive integer")
  assert_count_(2.5, "Number of replicates") |> expect_error("positive integer")
  assert_count_(c(1, 2), "Number of replicates", allow_na = TRUE) |>
    expect_error("positive integer")
})

test_that("ESS thinning keeps a draw of every chain", {
  # four chains of 100 draws stuck apart: the bulk ESS is near the chain count
  rho <- rep(1:4, each = 100L) + stats::rnorm(400L, sd = 0.01)
  draws <- posterior::draws_df(rho = rho, .nchains = 4L)
  thinned <- helper_read_param_table_(draws)
  nrow(thinned) |> expect_lt(10L)
  round(thinned$rho) |> expect_setequal(1:4)
})

test_that("ESS thinning reads the chains apart", {
  # four chains each drifting up one ramp: bulk ESS is near 6 over the chains,
  # near 20 over the chains read end to end as one
  rho <- rep(seq_len(100L), 4L) + rep(1:4 / 10, each = 100L)
  draws <- posterior::draws_df(rho = rho, .nchains = 4L)
  helper_read_param_table_(draws) |>
    nrow() |>
    expect_lt(10L)
})

###############################
###  UTILS : FACTORIZATION  ###
###############################

test_that("Data factorization handles default types", {
  result <- cpp_data_factorize(mk_types$data)
  names(result) |> expect_contains(c("data", "info", "cols", "status"))

  # Date encoded as numeric
  result$info$date_col$class |> expect_identical("Date")
  # character encoded with levels
  result$info$char_col$class |> expect_identical("character")
  result$info$char_col$levels |> expect_identical(c("male", "female"))
  # logical encoded
  result$info$logic_col$class |> expect_identical("logical")

  # numeric column untouched, and id, hhid, status skipped (no info entry)
  (c("num_col", "id", "hhid", "status") %in% names(result$info)) |>
    expect_all_false()

  # status extraction
  result$status |> expect_identical(1L)
  # columns preserved
  result$cols |> expect_identical(colnames(mk_types$data))
})

test_that("Data refactor restores original types", {
  factored <- cpp_data_factorize(mk_types$data)
  restored <- cpp_data_refactor(as_cpp_frame(factored),
                                factored$info, factored$cols)

  restored$date_col |> expect_identical(mk_types$data$date_col)
  restored$char_col |> expect_identical(mk_types$data$char_col)
  restored$logic_col |> expect_identical(mk_types$data$logic_col)
  # column order matches original headers
  colnames(restored)[seq_along(factored$cols)] |>
    expect_identical(factored$cols)
})

test_that("Data factorization codes a labelled status column", {
  hhdata <- data.frame(id = paste0("m", 1:4),
                       hhid = rep(c("h1", "h2"), each = 2),
                       status = c("positive", rep("negative", 3)),
                       stringsAsFactors = FALSE)

  # the negative level is unnamed: every individual reads as a case
  cpp_data_factorize(hhdata) |> expect_warning("escaped label")
  cpp_data_factorize(hhdata) |> expect_warning('escaped = "negative"')

  # naming it silences the guard and recovers the intended coding
  result <- cpp_data_factorize(hhdata, escaped = "negative")
  as.integer(result$data$status) |> expect_identical(c(1L, 0L, 0L, 0L))
  result$status |> expect_identical(1L)
  result$status_labels |> expect_identical("positive")

  # a round-tripped factor declares the escaped level even with no escaped row
  round_trip <- hhdata
  round_trip$status <- factor(c("symptomatic", "symptomatic",
                                "asymptomatic", "symptomatic"),
                              levels = c("escaped", "symptomatic",
                                         "asymptomatic"))
  result <- cpp_data_factorize(round_trip)
  as.integer(result$data$status) |> expect_identical(c(1L, 1L, 2L, 1L))

  # no recognised negative label: the hint falls back to the generic line
  unknown <- hhdata
  unknown$status <- c("red", "blue", "blue", "blue")
  cpp_data_factorize(unknown) |>
    expect_warning("Name the value carried by individuals")

  # integer statuses follow the documented convention and are left untouched
  coded <- hhdata
  coded$status <- c(1L, 2L, 2L, 2L)
  result <- cpp_data_factorize(coded)
  result$status |> expect_identical(c(1L, 2L))
  result$status_labels |> expect_identical(character())

  # an all-missing column carries no outcome, so it says nothing
  missing <- hhdata
  missing$status <- NA_character_
  expect_no_warning(result <- cpp_data_factorize(missing))
  result$data$status |> expect_all_na()
})

test_that("Data factorization flags discrete-valued columns", {
  discr <- cpp_data_factorize(mk_types$data)$discr
  discr |> expect_type("logical")
  names(discr) |> expect_identical(colnames(mk_types$data))

  discr[c("date_col", "int_col", "empty_int")] |> expect_all_true()
  discr[c("whole_col", "num_col", "char_col", "empty_col")] |>
    expect_all_false()
})

test_that("Data factorization leaves a missing cell missing", {
  factored <- cpp_data_factorize(mk_types$data)

  # an empty cell is no category of its own: it carries no level, and it codes
  # to NA so the engine reads it as an empty cell
  factored$info$datum$levels |> expect_identical(c("a"))
  factored$info$datum_2$levels |> expect_identical(c("TRUE"))
  factored$data$datum |> expect_identical(c("1", NA))
  factored$data$datum_2 |> expect_identical(c("1", NA))

  # a column carrying nothing holds no level at all
  factored$info$datum_3$levels |> expect_identical(character(0))
  factored$data$datum_3 |> expect_all_na()
})

test_that("A data column R can only spell quoted, survives the round trip", {
  hhdata <- dummy_household
  hhdata[["3"]] <- c(42L, 7L)
  hhdata[["if"]] <- c(1L, 2L)

  hhmodel <- new_housepi_model(hhdata, "test")
  # the engine reports the column under the user's own spelling
  colnames(hhmodel$data_wide) |> expect_identical(colnames(hhdata))
  hhmodel$data_wide[["3"]] |> expect_equal(hhdata[["3"]])
  names(inspect_house_state(hhmodel, "house")$data_wide) |>
    expect_identical(colnames(hhdata))

  # a definition sources it by its name, quoted
  get_latent(add_latent(hhmodel, hidden = `if`))$value |>
    expect_identical("if|D")

  # data-long keys cross the same way
  long <- data.frame(id = hhdata$id, check.names = FALSE)
  long[["3"]] <- c(1.5, 2.5)
  add_data_long(hhmodel, long) |>
    get_data_long() |>
    getElement("3") |>
    expect_identical(long[["3"]])
})

test_that("The housepi ggplot theme is set", {
  skip_if_not_installed("ggplot2")
  set_housepi_theme()
  ggplot2::theme_get()$panel.background$fill |> expect_identical("#FEF9EC")
})
