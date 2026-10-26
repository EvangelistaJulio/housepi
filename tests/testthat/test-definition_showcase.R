library(testthat)

test_that("Showcaser extract sourced values", {
  # matching source type returns just the name
  showcaser_extract("alpha|A", "A") |> expect_identical("alpha")
  # non-matching source type returns original
  showcaser_extract("alpha|A", "P") |> expect_identical("alpha|A")
  # NA input
  showcaser_extract(NA, "A") |> expect_identical(NA_character_)
  # numeric passthrough
  showcaser_extract(42, "A") |> expect_equal(42)
  # vectorized
  showcaser_extract(c("alpha|A", "beta|P", NA), "A") |>
    expect_identical(c("alpha", "beta|P", NA_character_))
})

test_that("Showcaser extract numeric handles double and pi constant", {
  showcaser_extract_numeric("3|N") |> expect_equal(3)
  showcaser_extract_numeric("-1.5|N") |> expect_equal(-1.5)
  # reserved constant pi
  showcaser_extract_numeric("pi|N") |> expect_equal(pi)
  # non-N type returns original string
  showcaser_extract_numeric("alpha|A") |> expect_identical("alpha|A")
  # a missing datum returns NA_real_
  showcaser_extract_numeric(NA_character_) |> expect_equal(NA_real_)
  # an unresolved flag stays logical, so the column keeps its type
  showcaser_extract_numeric(c(NA, FALSE)) |> expect_identical(c(NA, FALSE))
  # numeric passthrough
  showcaser_extract_numeric(42) |> expect_equal(42)
})

test_that("Showcaser extract func extracts function names", {
  # 3+ pipe segments: returns function name
  showcaser_extract_func("log|x|[alpha|A]") |> expect_identical("log")
  showcaser_extract_func("sum|x:x2|[a|A]:[b|A]") |> expect_identical("sum")
  # 2 segments: returns original
  showcaser_extract_func("alpha|A") |> expect_identical("alpha|A")
  # NA returns NA_real_
  showcaser_extract_func(NA) |> expect_equal(NA_real_)
  # numeric
  showcaser_extract_func(42) |> expect_equal(42)
})

test_that("helper_split_external_semicol splits at depth 0 colons", {
  # simple split
  helper_split_external_semicol("a:b") |> expect_identical(c("a", "b"))
  # nested brackets preserve inner colon
  helper_split_external_semicol("[a:b]:c") |> expect_identical(c("a:b", "c"))
  # single element (no colon)
  helper_split_external_semicol("abc") |> expect_identical("abc")
  # outer brackets stripped
  helper_split_external_semicol("[abc]") |> expect_identical("abc")
  # multiple nesting levels
  helper_split_external_semicol("[a:[b:c]]:d") |>
    expect_identical(c("a:[b:c]", "d"))
})

test_that("Showcase extract args extracts function arguments", {
  # encoded function with arguments
  result <- showcase_extract_args(list("log|x|[3|N]"))
  result[[1]] |> expect_equal(c(x = 3))

  # multi-argument function
  result <- showcase_extract_args(list("sum|x:x2|[1|N]:[2|N]"))
  result[[1]] |> expect_equal(c(x = 1, x2 = 2))

  # NA input
  result <- showcase_extract_args(list(NA))
  result[[1]] |> expect_equal(NA_real_)

  # non-function (2 pipe segments)
  result <- showcase_extract_args(list("alpha|A"))
  result[[1]] |> expect_equal(NA_real_)

  # non-character, non-NA element: returns c(element, NA_real_)
  result <- showcase_extract_args(list(7))
  result[[1]] |> expect_equal(c(7, NA_real_))
})

test_that("Housepi defininition showcaser converts encoded defs to readable", {
  # pure numeric column
  defs <- data.frame(const = "3|N", stringsAsFactors = FALSE)
  result <- housepi_defin_showcaser(defs)
  result$const |> expect_equal(3)

  # function column gets split into name + args
  defs <- data.frame(func = "log|x|[alpha|A]", min = "0|N",
                     stringsAsFactors = FALSE)
  result <- housepi_defin_showcaser(defs)
  result$func |> expect_identical("log")
  colnames(result) |> expect_contains("func_args")
  result$min |> expect_equal(0)
})
