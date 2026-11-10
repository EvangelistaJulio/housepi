library(testthat)

####  Helpers  ####

# data frame carrying one column of each storage type the engine encodes
mk_types <- local({
  hh_data <- data.frame(id = c("m1", "m2"),
                        hh_id = c("h1", "h1"),
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
  list(data = hh_data)
})

# factorized data coerced to the numeric frame the engine hands back
as_cpp_frame <- function(factored) {
  cpp_data <- suppressWarnings(data.frame(lapply(factored$data, as.numeric)))
  colnames(cpp_data) <- colnames(factored$data)
  cpp_data
}

####  Definitions  ####

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
  char_data_input_(exotic) |> expect_named(names(exotic))
})

test_that("Data factorization handles default types", {
  result <- cpp_data_factorize(mk_types$data)
  names(result) |> expect_contains(c("data", "info", "cols", "status"))

  # Date encoded as numeric
  result$info$date_col$class |> expect_identical("Date")
  # character encoded with levels
  result$info$char_col$class |> expect_identical("character")
  result$info$char_col$levels |> expect_identical(c("male", "female"))
  # logical encoded as FALSE 0, TRUE 1, with no levels
  result$info$logic_col$class |> expect_identical("logical")
  result$info$logic_col$levels |> expect_null()
  result$data$logic_col |> expect_identical(c("1", "0"))

  # numeric column untouched, and id, hh_id, status skipped (no info entry)
  (c("num_col", "id", "hh_id", "status") %in% names(result$info)) |>
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

test_that("A logical data column reads 0 and 1 in a definition", {
  hh_data <- data.frame(id = c("m1", "m2", "m3"), hh_id = "h1",
                        status = c(1, 0, 0), datum = c(TRUE, FALSE, TRUE),
                        datum_2 = 0)
  hh_model <- hh_data |>
    new_housepi_model("test", quiet = TRUE) |>
    add_observation(datum_2 = list(func = datum * 2))

  get_dataset(hh_model)$datum_2 |> expect_equal(c(2, 0, 2))
  get_dataset(hh_model)$datum |> expect_identical(hh_data$datum)
})

test_that("Data factorization codes a labelled status column", {
  hh_data <- data.frame(id = paste0("m", 1:4),
                        hh_id = rep(c("h1", "h2"), each = 2),
                        status = c("positive", rep("negative", 3)),
                        stringsAsFactors = FALSE)

  # the negative level is unnamed: every individual reads as a case
  cpp_data_factorize(hh_data) |> expect_warning("escaped label")
  cpp_data_factorize(hh_data) |> expect_warning('escaped = "negative"')

  # naming it silences the guard and recovers the intended coding
  result <- cpp_data_factorize(hh_data, escaped = "negative")
  as.integer(result$data$status) |> expect_identical(c(1L, 0L, 0L, 0L))
  result$status |> expect_identical(1L)
  result$status_labels |> expect_identical("positive")

  # a round-tripped factor declares the escaped level even with no escaped row
  round_trip <- hh_data
  round_trip$status <- factor(c("symptomatic", "symptomatic",
                                "asymptomatic", "symptomatic"),
                              levels = c("escaped", "symptomatic",
                                         "asymptomatic"))
  result <- cpp_data_factorize(round_trip)
  as.integer(result$data$status) |> expect_identical(c(1L, 1L, 2L, 1L))

  # no recognised negative label: the hint falls back to the generic line
  unknown <- hh_data
  unknown$status <- c("red", "blue", "blue", "blue")
  cpp_data_factorize(unknown) |>
    expect_warning("Name the value carried by individuals")

  # integer statuses follow the documented convention and are left untouched
  coded <- hh_data
  coded$status <- c(1L, 2L, 2L, 2L)
  result <- cpp_data_factorize(coded)
  result$status |> expect_identical(c(1L, 2L))
  result$status_labels |> expect_identical(character())

  # an all-missing column carries no outcome, so it says nothing
  missing <- hh_data
  missing$status <- NA_character_
  expect_no_warning(result <- cpp_data_factorize(missing))
  result$data$status |> expect_all_na()
})

test_that("Data factorization flags discrete-valued columns", {
  discr <- cpp_data_factorize(mk_types$data)$discr
  discr |> expect_type("logical")
  discr |> expect_named(colnames(mk_types$data))

  discr[c("date_col", "int_col", "logic_col", "empty_int")] |>
    expect_all_true()
  discr[c("whole_col", "num_col", "char_col", "empty_col")] |>
    expect_all_false()
})

test_that("Data factorization leaves a missing cell missing", {
  factored <- cpp_data_factorize(mk_types$data)

  # an empty cell is no category of its own: it carries no level, and it codes
  # to NA so the engine reads it as an empty cell
  factored$info$datum$levels |> expect_identical(c("a"))
  factored$data$datum |> expect_identical(c("1", NA))
  factored$data$datum_2 |> expect_identical(c("1", NA))

  # a column carrying nothing holds no level at all
  factored$info$datum_3$levels |> expect_identical(character(0))
  factored$data$datum_3 |> expect_all_na()
})

test_that("A data column R can only spell quoted, survives the round trip", {
  hh_data <- dummy_household
  hh_data[["3"]] <- c(42L, 7L)
  hh_data[["if"]] <- c(1L, 2L)

  hh_model <- new_housepi_model(hh_data, "test")
  # the engine reports the column under the user's own spelling
  colnames(get_dataset(hh_model)) |> expect_identical(colnames(hh_data))
  get_dataset(hh_model)[["3"]] |> expect_equal(hh_data[["3"]])
  names(inspect_state(hh_model, "house")$dataset) |>
    expect_identical(colnames(hh_data))

  # a definition sources it by its name, quoted
  get_latent(add_latent(hh_model, hidden = `if`))$func |>
    expect_identical("if|D")

  # long dataset keys cross the same way
  long <- data.frame(id = hh_data$id, check.names = FALSE)
  long[["3"]] <- c(1.5, 2.5)
  add_long_dataset(hh_model, long) |>
    get_long_dataset() |>
    getElement("3") |>
    expect_identical(long[["3"]])
})
