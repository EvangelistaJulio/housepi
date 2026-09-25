library(testthat)

test_that("Helper extract definition work as expected", {
  # two adjacent blocks
  result <- helper_extract_defin_blocks("[a][b]")
  result |> expect_identical(list("[a]", "[b]"))

  # nested brackets count as one block
  result <- helper_extract_defin_blocks("[a[b]c]")
  result |> expect_identical(list("[a[b]c]"))

  # no brackets returns empty list
  result <- helper_extract_defin_blocks("abc")
  result |> expect_identical(list())
})

test_that("Helper evaluate definition block evaluates conditions by type", {
  # block with mixed types: extract matching
  result <- helper_evaluate_defin_block("[[alpha|A]:[beta|P]]", "A")
  result |> expect_identical("alpha")

  result <- helper_evaluate_defin_block("[[alpha|A]:[beta|P]]", "P")
  result |> expect_identical("beta")

  # non-matching type returns empty
  result <- helper_evaluate_defin_block("[alpha|A]", "P")
  result |> expect_identical(c())
})

test_that("housepi_walk_definition handles recursive and terminal nodes", {
  # terminal match
  result <- housepi_walk_definition("alpha|A", "A")
  result |> expect_identical("alpha")

  # terminal non-match
  result <- housepi_walk_definition("alpha|A", "P")
  result |> expect_null()

  # block-based definition
  result <- housepi_walk_definition("[alpha|A]:[beta|A]", "A")
  ("alpha" %in% result || "beta" %in% result) |> expect_true()

  # NA / NULL / empty string return NULL
  housepi_walk_definition(NA, "A") |> expect_null()
  housepi_walk_definition(NULL, "A") |> expect_null()
  housepi_walk_definition("", "A") |> expect_null()

  # deduplication
  result <- housepi_walk_definition("[alpha|A]:[alpha|A]", "A")
  length(result) |> expect_identical(1L)
  result |> expect_identical("alpha")
})
