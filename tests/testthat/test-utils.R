library(testthat)

test_that("Character to string list formats vectors correctly", {
  expect_equal(housepi:::character_to_string_list("x"), "'x'.")
  str_list <- letters[1:5]
  # truncation at max_print
  expect_equal(housepi:::character_to_string_list(str_list, max_print = 2),
               "'a', 'b', ...")
  # default max_print=4, five elements triggers ellipsis
  expect_equal(housepi:::character_to_string_list(str_list),
               "'a', 'b', 'c', 'd', ...")
  # max_print larger than vector
  expect_equal(housepi:::character_to_string_list(str_list, max_print = 10),
               "'a', 'b', 'c', 'd', 'e'.")
})

test_that("Character to phrase joins words with and", {
  expect_null(housepi:::character_to_phrase(character(0)))
  expect_equal(housepi:::character_to_phrase("alpha"), "alpha")
  expect_equal(housepi:::character_to_phrase(c("alpha", "beta")),
               "alpha and beta")
  expect_equal(housepi:::character_to_phrase(c("a", "b", "c")),
               "a, b and c")
})

test_that("Find closest match suggests close options", {
  # typo within distance 1
  result <- housepi:::find_closest_match("alpa", c("alpha", "beta"))
  expect_true(any(grepl("Did you mean", result)))
  expect_true(any(grepl("alpha", result)))

  # min_dist = 0 disables suggestions
  expect_null(housepi:::find_closest_match("alpa", c("alpha"), min_dist = 0))

  # no close match returns available options
  result <- housepi:::find_closest_match("zzzzz", c("alpha", "beta"))
  expect_true(any(grepl("Available option", result)))

  # empty options
  expect_equal(housepi:::find_closest_match("x", character(0)), c())
})

test_that("helper_clone deep clones data.frame, list, and scalars", {
  # data.frame independence
  original <- data.frame(x = 1:3)
  cloned <- housepi:::helper_clone(original)
  cloned$x[1] <- 99
  expect_equal(original$x[1], 1)

  # list of data.frames recursively cloned
  lst <- list(a = data.frame(x = 1), b = data.frame(y = 2))
  cloned_lst <- housepi:::helper_clone(lst)
  cloned_lst$a$x[1] <- 99
  expect_equal(lst$a$x[1], 1)

  # scalar passthrough
  expect_equal(housepi:::helper_clone(42), 42)
  expect_equal(housepi:::helper_clone("hello"), "hello")
})

test_that("get_val_ retrieves from list with default", {
  lst <- list(a = 1, b = "two")
  expect_equal(housepi:::get_val_(lst, "a"), 1)
  expect_equal(housepi:::get_val_(lst, "b"), "two")
  expect_true(is.na(housepi:::get_val_(lst, "c")))
  expect_equal(housepi:::get_val_(lst, "c", default = 0), 0)
})

test_that("cpp_data_factorize handles dates, factors, characters, logicals", {
  df <- data.frame(
    id = c("p1", "p2"),
    house_id = c("h1", "h1"),
    status = c(1, 0),
    date_col = as.Date(c("2024-01-01", "2024-01-02")),
    char_col = c("male", "female"),
    logic_col = c(TRUE, FALSE),
    num_col = c(1.5, 2.5),
    stringsAsFactors = FALSE
  )

  result <- housepi:::cpp_data_factorize(df)

  # structure
  expect_true(is.list(result))
  expect_true(all(c("data", "info", "cols", "status") %in% names(result)))

  # Date encoded as numeric
  expect_equal(result$info$date_col$class, "Date")
  expect_true(is.numeric(as.numeric(result$data$date_col)))

  # character encoded with levels
  expect_equal(result$info$char_col$class, "character")
  expect_equal(result$info$char_col$levels, c("male", "female"))

  # logical encoded
  expect_equal(result$info$logic_col$class, "logical")

  # numeric column untouched (no info entry)
  expect_false("num_col" %in% names(result$info))

  # id, house_id, status are skipped
  expect_false("id" %in% names(result$info))
  expect_false("house_id" %in% names(result$info))
  expect_false("status" %in% names(result$info))

  # status extraction
  expect_equal(result$status, 1L)

  # columns preserved
  expect_equal(result$cols, colnames(df))
})

test_that("cpp_data_refactor restores original types", {
  df <- data.frame(
    id = c("p1", "p2"),
    house_id = c("h1", "h1"),
    status = c(1, 0),
    date_col = as.Date(c("2024-01-01", "2024-01-02")),
    char_col = c("male", "female"),
    logic_col = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  factored <- housepi:::cpp_data_factorize(df)
  # simulate what C++ returns: numeric data.frame
  cpp_data <- suppressWarnings(data.frame(lapply(factored$data, as.numeric)))
  colnames(cpp_data) <- colnames(factored$data)

  restored <- housepi:::cpp_data_refactor(cpp_data, factored$info, factored$cols)

  expect_equal(restored$date_col, df$date_col)
  expect_equal(restored$char_col, df$char_col)
  expect_equal(restored$logic_col, df$logic_col)
  # column order matches original headers
  expect_equal(colnames(restored)[seq_along(factored$cols)], factored$cols)
})

test_that("stop_h and warning_h format multi-line messages", {
  expect_error(housepi:::stop_h("line1", "line2"), regexp = "line1\nline2")
  expect_warning(housepi:::warning_h("w1", "w2"), regexp = "w1\nw2")
})

test_that("rnd_string produces string of correct length", {
  expect_equal(nchar(housepi:::rnd_string()), 8)
  expect_equal(nchar(housepi:::rnd_string(12)), 12)
  # different calls produce different strings (with high probability)
  expect_false(housepi:::rnd_string() == housepi:::rnd_string())
})
