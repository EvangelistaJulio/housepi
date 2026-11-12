library(testthat)

####  Helpers  ####

# data frame carrying one column of each storage type the engine encodes
mk_types <- local({
  hh_data <- data.frame(id = c("m1", "m2"),
                        hh_id = c("h1", "h1"),
                        status = c(1, 0),
                        date_col = as.Date(c("2024-01-01", "2024-01-02")),
                        char_col = c("male", "female"),
                        fct_col = factor(c("b", "a"),
                                         levels = c("c", "b", "a")),
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

# factorized data bound into the frame the engine hands back
as_cpp_frame <- function(factored) cbind(factored$data$ids, factored$data$cells)

####  Definitions  ####

test_that("Data factorization handles default types", {
  result <- cpp_data_factorize(mk_types$data)
  names(result) |> expect_contains(c("data", "info", "cols", "status"))

  # Date encoded as numeric
  result$info$date_col$class |> expect_identical("Date")
  # character encoded with levels, sorted whatever the row order
  result$info$char_col$class |> expect_identical("character")
  result$info$char_col$levels |> expect_identical(c("female", "male"))
  # logical encoded as FALSE 0, TRUE 1, with no levels
  result$info$logic_col$class |> expect_identical("logical")
  result$info$logic_col$levels |> expect_null()
  result$data$cells$logic_col |> expect_identical(c(1, 0))

  # double column untouched, character id and hh_id, status (no info entry)
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
  # a column a definition writes keeps the type the engine gives it
  cpp_data_refactor(as_cpp_frame(factored), factored$info, factored$cols,
                    produced = "int_col")$int_col |>
    expect_type("double")
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

test_that("A categorical column reads its code whatever the row order", {
  hh_data <- data.frame(id = c("m1", "m2", "m3"), hh_id = "h1",
                        status = c(1, 0, 0), datum = c("b", "c", "a"),
                        datum_2 = 0)
  code_of <- function(hh_data) {
    hh_model <- hh_data |>
      new_housepi_model("test", quiet = TRUE) |>
      add_observation(datum_2 = list(func = datum * 1))
    stats::setNames(get_dataset(hh_model)$datum_2, hh_data$datum)
  }
  code_of(hh_data[3:1, ])[hh_data$datum] |>
    expect_identical(code_of(hh_data))
})

test_that("Status labels take one order whatever the row order", {
  hh_data <- data.frame(hh_id = "h1", datum = 1,
                        status = c("symptomatic", "asymptomatic", "escaped"))
  # character labels -> sorted, after the escaped label
  for (rows in list(1:3, 3:1)) {
    hh_model <- new_housepi_model(hh_data[rows, ], "final-size")
    levels(get_dataset(hh_model)$status) |>
      expect_identical(c("escaped", "asymptomatic", "symptomatic"))
  }
})

test_that("Data round-trips through the engine unchanged", {
  hh_data <- data.frame(id = 1:4, hh_id = c(1L, 1L, 2L, 2L),
                        status = c(1L, 0L, 0L, 0L), datum = c(3L, 1L, NA, 2L),
                        datum_2 = c(0.1 + 0.2, 1 / 3, NA, 1e-20 / 3),
                        datum_3 = factor(c("b", "a", NA, "b"),
                                         levels = c("c", "b", "a")))
  hh_model <- new_housepi_model(hh_data, "test", quiet = TRUE)
  get_dataset(hh_model) |> expect_identical(hh_data)

  # a factor identifier keeps its levels, a character one stays character
  hh_data_2 <- transform(hh_data, id = paste0("m", id),
                         hh_id = factor(hh_id, levels = c(2, 1)))
  new_housepi_model(hh_data_2, "test", quiet = TRUE) |>
    get_dataset() |>
    expect_identical(hh_data_2)
})

test_that("A Date column stored as integers reads back as a Date", {
  # integer storage, as data.table::fread() reads dates
  hh_data <- data.frame(id = c("m1", "m2"), hh_id = "h1", status = c(1L, 0L),
                        datum = structure(c(19723L, 19724L), class = "Date"))
  new_housepi_model(hh_data, "test", quiet = TRUE) |>
    get_dataset() |>
    getElement("datum") |>
    expect_identical(as.Date(c("2024-01-01", "2024-01-02")))
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
  result$data$cells$status |> expect_identical(c(1, 0, 0, 0))
  result$status |> expect_identical(1L)
  result$status_labels |> expect_identical("positive")

  # a round-tripped factor declares the escaped level even with no escaped row
  round_trip <- hh_data
  round_trip$status <- factor(c("symptomatic", "symptomatic",
                                "asymptomatic", "symptomatic"),
                              levels = c("escaped", "symptomatic",
                                         "asymptomatic"))
  result <- cpp_data_factorize(round_trip)
  result$data$cells$status |> expect_identical(c(1, 1, 2, 1))

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
  # a fractional status is no code
  coded$status <- c(1.5, 2, 2, 2)
  cpp_data_factorize(coded) |> expect_error("1.5.*not a whole number")

  # an all-missing column carries no outcome, so it says nothing
  missing <- hh_data
  missing$status <- NA_character_
  expect_no_warning(result <- cpp_data_factorize(missing))
  result$data$cells$status |> expect_all_na()
})

test_that("Data factorization flags discrete-valued columns", {
  discr <- cpp_data_factorize(mk_types$data)$discr
  discr |> expect_type("logical")
  discr |> expect_named(colnames(mk_types$data))

  # a categorical column reads integer codes
  discr[c("date_col", "int_col", "logic_col", "empty_int", "char_col",
          "fct_col")] |>
    expect_all_true()
  discr[c("whole_col", "num_col", "empty_col")] |>
    expect_all_false()
})

test_that("Data factorization leaves a missing cell missing", {
  factored <- cpp_data_factorize(mk_types$data)

  # an empty cell is no category of its own: it carries no level, and it codes
  # to NA so the engine reads it as an empty cell
  factored$info$datum$levels |> expect_identical(c("a"))
  factored$data$cells$datum |> expect_identical(c(1, NA))
  factored$data$cells$datum_2 |> expect_identical(c(1, NA))

  # a column carrying nothing holds no level at all
  factored$info$datum_3$levels |> expect_identical(character(0))
  factored$data$cells$datum_3 |> expect_all_na()
})

test_that("A data column R can only spell quoted, survives the round trip", {
  hh_data <- dummy_household
  hh_data[["3"]] <- c(42L, 7L)
  hh_data[["if"]] <- c(1L, 2L)

  hh_model <- new_housepi_model(hh_data, "test")
  # the engine reports the column under the user's own spelling
  colnames(get_dataset(hh_model)) |> expect_identical(colnames(hh_data))
  get_dataset(hh_model)[["3"]] |> expect_identical(hh_data[["3"]])
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
