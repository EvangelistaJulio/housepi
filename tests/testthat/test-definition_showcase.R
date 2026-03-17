library(testthat)

test_that("Showcaser extract sourced values", {
  # matching source type returns just the name
  expect_equal(housepi:::showcaser_extract("alpha|A", "A"), "alpha")
  # non-matching source type returns original
  expect_equal(housepi:::showcaser_extract("alpha|A", "P"), "alpha|A")
  # NA input
  expect_equal(housepi:::showcaser_extract(NA, "A"), NA_character_)
  # numeric passthrough
  expect_equal(housepi:::showcaser_extract(42, "A"), 42)
  # vectorized
  expect_equal(
    housepi:::showcaser_extract(c("alpha|A", "beta|P", NA), "A"),
    c("alpha", "beta|P", NA_character_)
  )
})

test_that("Showcaser extract numeric handles double and pi constant", {
  expect_equal(housepi:::showcaser_extract_numeric("3|N"), 3)
  expect_equal(housepi:::showcaser_extract_numeric("-1.5|N"), -1.5)
  # reserved constant pi
  expect_equal(housepi:::showcaser_extract_numeric("pi|N"), pi)
  # non-N type returns original string
  expect_equal(housepi:::showcaser_extract_numeric("alpha|A"), "alpha|A")
  # NA returns NA_real_
  expect_equal(housepi:::showcaser_extract_numeric(NA), NA_real_)
  # numeric passthrough
  expect_equal(housepi:::showcaser_extract_numeric(42), 42)
})

test_that("Showcaser extract func extracts function names", {
  # 3+ pipe segments: returns function name
  expect_equal(housepi:::showcaser_extract_func("log|x|{alpha|A}"), "log")
  expect_equal(housepi:::showcaser_extract_func("sum|x:x2|{a|A}:{b|A}"), "sum")
  # 2 segments: returns original
  expect_equal(housepi:::showcaser_extract_func("alpha|A"), "alpha|A")
  # NA returns NA_real_
  expect_equal(housepi:::showcaser_extract_func(NA), NA_real_)
})

test_that("helper_split_external_semicol splits at depth 0 colons", {
  # simple split
  expect_equal(housepi:::helper_split_external_semicol("a:b"), c("a", "b"))
  # nested braces preserve inner colon
  expect_equal(housepi:::helper_split_external_semicol("{a:b}:c"),
               c("a:b", "c"))
  # single element (no colon)
  expect_equal(housepi:::helper_split_external_semicol("abc"), "abc")
  # outer braces stripped
  expect_equal(housepi:::helper_split_external_semicol("{abc}"), "abc")
  # multiple nesting levels
  expect_equal(
    housepi:::helper_split_external_semicol("{a:{b:c}}:d"),
    c("a:{b:c}", "d")
  )
})

test_that("Showcase extract args extracts function arguments", {
  # encoded function with arguments
  result <- housepi:::showcase_extract_args(list("log|x|{3|N}"))
  expect_equal(result[[1]], c(x = 3))

  # multi-argument function
  result <- housepi:::showcase_extract_args(list("sum|x:x2|{1|N}:{2|N}"))
  expect_equal(result[[1]], c(x = 1, x2 = 2))

  # NA input
  result <- housepi:::showcase_extract_args(list(NA))
  expect_equal(result[[1]], NA_real_)

  # non-function (2 pipe segments)
  result <- housepi:::showcase_extract_args(list("alpha|A"))
  expect_equal(result[[1]], NA_real_)
})

test_that("Housepi defininition showcaser converts encoded defs to readable", {
  # pure numeric column
  defs <- data.frame(value = "3|N", stringsAsFactors = FALSE)
  result <- housepi:::housepi_defin_showcaser(defs)
  expect_equal(result$value, 3)

  # function column gets split into name + args
  defs <- data.frame(func = "log|x|{alpha|A}", min = "0|N",
                     stringsAsFactors = FALSE)
  result <- housepi:::housepi_defin_showcaser(defs)
  expect_equal(result$func, "log")
  expect_true("func_args" %in% colnames(result))
  expect_equal(result$min, 0)
})
