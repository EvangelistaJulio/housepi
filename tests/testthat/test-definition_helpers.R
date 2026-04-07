library(testthat)

test_that("Helper extract definition work as expected", {
  # two adjacent blocks
  result <- helper_extract_defin_blocks("{a}{b}")
  expect_equal(result, list("{a}", "{b}"))

  # nested braces count as one block
  result <- helper_extract_defin_blocks("{a{b}c}")
  expect_equal(result, list("{a{b}c}"))

  # no braces returns empty list
  result <- helper_extract_defin_blocks("abc")
  expect_equal(result, list())
})

test_that("Helper evaluate definition block evaluates conditions by type", {
  # block with mixed types: extract matching
  result <- helper_evaluate_defin_block("{{alpha|A}:{beta|P}}", "A")
  expect_equal(result, "alpha")

  result <- helper_evaluate_defin_block("{{alpha|A}:{beta|P}}", "P")
  expect_equal(result, "beta")

  # non-matching type returns empty
  result <- helper_evaluate_defin_block("{alpha|A}", "P")
  expect_equal(result, c())
})

test_that("housepi_walk_definition handles recursive and terminal nodes", {
  # terminal match
  result <- housepi_walk_definition("alpha|A", "A")
  expect_equal(result, "alpha")

  # terminal non-match
  result <- housepi_walk_definition("alpha|A", "P")
  expect_null(result)

  # block-based definition
  result <- housepi_walk_definition("{alpha|A}:{beta|A}", "A")
  expect_true("alpha" %in% result || "beta" %in% result)

  # NA / NULL / empty string return NULL
  expect_null(housepi_walk_definition(NA, "A"))
  expect_null(housepi_walk_definition(NULL, "A"))
  expect_null(housepi_walk_definition("", "A"))

  # deduplication
  result <- housepi_walk_definition("{alpha|A}:{alpha|A}", "A")
  expect_equal(length(result), 1)
  expect_equal(result, "alpha")
})
