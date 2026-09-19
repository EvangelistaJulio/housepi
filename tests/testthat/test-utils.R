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
  expect_invisible(overwrite(path, path))

  # destination directory does not exist -> file.copy returns FALSE
  overwrite(path, file.path(tempfile(), "missing", "x")) |>
    expect_error("user permissions")
})

test_that("Character to phrase produce phrases as expected", {
  expect_null(character_to_phrase(character(0)))
  expect_equal(character_to_phrase("alpha"), "alpha")
  expect_equal(character_to_phrase(c("alpha", "beta")), "alpha and beta")
  expect_equal(character_to_phrase(c("a", "b", "c")), "a, b and c")
})

test_that("Find closest match suggests close options", {
  # typo within distance 1
  result <- find_closest_match("alpa", c("alpha", "beta"))
  expect_true(any(grepl("Did you mean", result)))
  expect_true(any(grepl("alpha", result)))

  # min_dist = 0 disables suggestions
  expect_null(find_closest_match("alpa", c("alpha"), min_dist = 0))

  # no close match returns available options
  result <- find_closest_match("unknown", c("alpha", "beta"))
  expect_true(any(grepl("Available option", result)))

  # empty options
  expect_equal(find_closest_match("x", character(0)), c())
})

test_that("Deep clone helper work as expected on default types", {
  # scalar passthrough
  expect_equal(helper_clone(42), 42)
  expect_equal(helper_clone("hello"), "hello")

  # data.frame independence
  original <- data.frame(x = 1:3)
  cloned <- helper_clone(original)
  cloned$x[1] <- 99
  expect_equal(original$x[1], 1)
  expect_equal(cloned$x[1], 99)

  # list of data.frames recursively cloned
  frames <- list(a = data.frame(x = 1), b = data.frame(y = 2))
  cloned <- helper_clone(frames)
  cloned$a$x[1] <- 99
  expect_equal(frames$a$x[1], 1)
  expect_equal(cloned$a$x[1], 99)
})

test_that("Setup options are formatted as character", {
  # id_draws with multiple values
  result <- char_setup_input_(list(steps = 100, burnin = 10, thin = 1,
                                   chains = 1, id_draws = c("m1", "m2")))
  expect_s3_class(result, "data.frame")
  expect_equal(result$id_draws, "m1:m2")
  expect_true(all(vapply(result, is.character, logical(1))))
  # id_draws are dropped if NULL
  char_setup_input_(list(steps = 100L, id_draws = NULL)) |>
    getElement("id_draws") |>
    expect_identical("")
})

test_that("Char data input handles NULL and data.frame", {
  expect_null(char_data_input_(NULL))
  result <- char_data_input_(data.frame(x = 1:3, y = c(4.5, 5.5, 6.5)))
  expect_s3_class(result, "data.frame")
  expect_true(all(vapply(result, is.character, logical(1))))

  # column names are engine keys: they cross verbatim, never through make.names
  exotic <- data.frame(a = 1, check.names = FALSE)
  exotic[["3"]] <- 2
  exotic[["if"]] <- 3
  exotic[["dose.2"]] <- 4
  expect_equal(names(char_data_input_(exotic)), names(exotic))
})

###############################
###  UTILS : FACTORIZATION  ###
###############################

test_that("Data factorization handles default types", {
  result <- cpp_data_factorize(mk_types$data)
  expect_true(all(c("data", "info", "cols", "status") %in% names(result)))

  # Date encoded as numeric
  expect_equal(result$info$date_col$class, "Date")
  # character encoded with levels
  expect_equal(result$info$char_col$class, "character")
  expect_equal(result$info$char_col$levels, c("male", "female"))
  # logical encoded
  expect_equal(result$info$logic_col$class, "logical")

  # numeric column untouched, and id, hhid, status skipped (no info entry)
  any(c("num_col", "id", "hhid", "status") %in% names(result$info)) |>
    expect_false()

  # status extraction
  expect_equal(result$status, 1L)
  # columns preserved
  expect_equal(result$cols, colnames(mk_types$data))
})

test_that("Data refactor restores original types", {
  factored <- cpp_data_factorize(mk_types$data)
  restored <- cpp_data_refactor(as_cpp_frame(factored),
                                factored$info, factored$cols)

  expect_equal(restored$date_col, mk_types$data$date_col)
  expect_equal(restored$char_col, mk_types$data$char_col)
  expect_equal(restored$logic_col, mk_types$data$logic_col)
  # column order matches original headers
  expect_equal(colnames(restored)[seq_along(factored$cols)], factored$cols)
})

test_that("Data factorization codes a labelled status column", {
  hhdata <- data.frame(id = paste0("m", 1:4),
                       hhid = rep(c("h1", "h2"), each = 2),
                       status = c("positive", rep("negative", 3)),
                       stringsAsFactors = FALSE)

  # the negative level is unnamed: every individual reads as a case
  expect_warning(cpp_data_factorize(hhdata), "escaped label")
  expect_warning(cpp_data_factorize(hhdata), 'escaped = "negative"')

  # naming it silences the guard and recovers the intended coding
  result <- cpp_data_factorize(hhdata, escaped = "negative")
  expect_equal(as.integer(result$data$status), c(1L, 0L, 0L, 0L))
  expect_equal(result$status, 1L)
  expect_equal(result$status_labels, "positive")

  # a round-tripped factor declares the escaped level even with no escaped row
  round_trip <- hhdata
  round_trip$status <- factor(c("symptomatic", "symptomatic",
                                "asymptomatic", "symptomatic"),
                              levels = c("escaped", "symptomatic",
                                         "asymptomatic"))
  result <- cpp_data_factorize(round_trip)
  expect_equal(as.integer(result$data$status), c(1L, 1L, 2L, 1L))

  # no recognised negative label: the hint falls back to the generic line
  unknown <- hhdata
  unknown$status <- c("red", "blue", "blue", "blue")
  cpp_data_factorize(unknown) |>
    expect_warning("Name the value carried by individuals")

  # integer statuses follow the documented convention and are left untouched
  coded <- hhdata
  coded$status <- c(1L, 2L, 2L, 2L)
  result <- cpp_data_factorize(coded)
  expect_equal(result$status, c(1L, 2L))
  expect_equal(result$status_labels, character())

  # an all-missing column carries no outcome, so it says nothing
  missing <- hhdata
  missing$status <- NA_character_
  expect_no_warning(result <- cpp_data_factorize(missing))
  expect_true(all(is.na(result$data$status)))
})

test_that("Data factorization flags discrete-valued columns", {
  discr <- cpp_data_factorize(mk_types$data)$discr
  expect_type(discr, "logical")
  expect_equal(names(discr), colnames(mk_types$data))

  expect_true(all(discr[c("date_col", "int_col", "empty_int")]))
  expect_false(any(discr[c("whole_col", "num_col", "char_col", "empty_col")]))
})

test_that("Data factorization leaves a missing cell missing", {
  factored <- cpp_data_factorize(mk_types$data)

  # an empty cell is no category of its own: it carries no level, and it codes
  # to NA so the engine reads it as an empty cell
  expect_equal(factored$info$datum$levels, c("a"))
  expect_equal(factored$info$datum_2$levels, c("TRUE"))
  expect_equal(factored$data$datum, c("1", NA))
  expect_equal(factored$data$datum_2, c("1", NA))

  # a column carrying nothing holds no level at all
  expect_equal(factored$info$datum_3$levels, character(0))
  expect_true(all(is.na(factored$data$datum_3)))
})

test_that("A data column R can only spell quoted, survives the round trip", {
  hhdata <- dummy_household
  hhdata[["3"]] <- c(42L, 7L)
  hhdata[["if"]] <- c(1L, 2L)

  hhmodel <- new_housepi_model(hhdata, "test")
  # the engine reports the column under the user's own spelling
  expect_equal(colnames(hhmodel$data_wide), colnames(hhdata))
  expect_equal(hhmodel$data_wide[["3"]], hhdata[["3"]])
  expect_equal(names(inspect_house_state(hhmodel, "house")$data_wide),
               colnames(hhdata))

  # a definition sources it, quoted or bare, and finds the datum
  add_latent(hhmodel, hidden = 3@D) |> expect_no_error()
  add_latent(hhmodel, hidden = `3`@D) |> expect_no_error()

  # data-long keys cross the same way
  long <- data.frame(id = hhdata$id, check.names = FALSE)
  long[["3"]] <- c(1.5, 2.5)
  add_data_long(hhmodel, long) |>
    get_data_long() |>
    getElement("3") |>
    expect_equal(long[["3"]])
})

test_that("The housepi ggplot theme is set without errors", {
  expect_no_error(housepi_theme_set())
})
