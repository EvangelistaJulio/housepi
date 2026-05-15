library(testthat)

test_that("Character to string list formats vectors correctly", {
  expect_equal(character_to_string_list("x"), "'x'.")
  str_list <- letters[1:5]
  # truncation at max_print
  expect_equal(character_to_string_list(str_list, max_print = 2),
               "'a', 'b', ...")
  # default max_print=4, five elements triggers ellipsis
  expect_equal(character_to_string_list(str_list),
               "'a', 'b', 'c', 'd', ...")
  # max_print larger than vector
  expect_equal(character_to_string_list(str_list, max_print = 10),
               "'a', 'b', 'c', 'd' and 'e'.")
})

test_that("Character to phrase produce phrases as expected", {
  expect_null(character_to_phrase(character(0)))
  expect_equal(character_to_phrase("alpha"), "alpha")
  expect_equal(character_to_phrase(c("alpha", "beta")),
               "alpha and beta")
  expect_equal(character_to_phrase(c("a", "b", "c")),
               "a, b and c")
})

test_that("Find closest match suggests close options", {
  # typo within distance 1
  result <- find_closest_match("alpa", c("alpha", "beta"))
  expect_true(any(grepl("Did you mean", result)))
  expect_true(any(grepl("alpha", result)))

  # min_dist = 0 disables suggestions
  expect_null(find_closest_match("alpa", c("alpha"), min_dist = 0))

  # no close match returns available options
  result <- find_closest_match("ocelot", c("alpha", "beta"))
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

  # list of data.frames recursively cloned
  lst <- list(a = data.frame(x = 1), b = data.frame(y = 2))
  cloned_lst <- helper_clone(lst)
  cloned_lst$a$x[1] <- 99
  expect_equal(lst$a$x[1], 1)
})

# UPDATE: I will move factorization in cpp
df_factory <- function() {
  data.frame(id = c("p1", "p2"),
             house_id = c("h1", "h1"),
             status = c(1, 0),
             date_col = as.Date(c("2024-01-01", "2024-01-02")),
             char_col = c("male", "female"),
             logic_col = c(TRUE, FALSE),
             num_col = c(1.5, 2.5),
             stringsAsFactors = FALSE)
}

test_that("cpp data factorization handles default types", {
  df <- df_factory()
  result <- cpp_data_factorize(df)

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

test_that("cpp data refactor restores original types", {
  df <- df_factory()
  factored <- cpp_data_factorize(df)
  cpp_data <- suppressWarnings(data.frame(lapply(factored$data, as.numeric)))
  colnames(cpp_data) <- colnames(factored$data)

  restored <- cpp_data_refactor(cpp_data, factored$info, factored$cols)

  expect_equal(restored$date_col, df$date_col)
  expect_equal(restored$char_col, df$char_col)
  expect_equal(restored$logic_col, df$logic_col)
  # column order matches original headers
  expect_equal(colnames(restored)[seq_along(factored$cols)], factored$cols)
})

test_that("Char setup input formats setup options correctly", {
  # id_draws with multiple values
  setup <- list(steps = 100, burnin = 10, thin = 1,
                chains = 1, id_draws = c("id1", "id2"))
  result <- char_setup_input_(setup)
  expect_true(is.data.frame(result))
  expect_equal(result$id_draws, "id1:id2")

  # id_draws NULL
  setup2 <- list(steps = 100L, id_draws = NULL)
  result2 <- char_setup_input_(setup2)
  expect_true(is.na(result2$id_draws))

  # all columns are character
  expect_true(all(sapply(result, is.character)))
})

# UPDATE: once data becomes variant, no need to cast to character
test_that("Char data input handles NULL and data.frame", {
  expect_null(char_data_input_(NULL))
  df <- data.frame(x = 1:3, y = c(4.5, 5.5, 6.5))
  result <- char_data_input_(df)
  expect_true(is.data.frame(result))
  expect_true(all(sapply(result, is.character)))
})

test_that("housepi_theme_set is called without errors", {
  expect_no_error(housepi_theme_set())
})
