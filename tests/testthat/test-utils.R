library(testthat)

test_that("Rlang standard functions work as expected", {
  expect_error(housepi:::DefEncoder$new(~log ( alpha @ B )),
               regexp = "Invalid type 'B'")
  expect_equal(housepi:::DefEncoder$new(~log(alpha))$x,
               "log|x|{alpha|A}")
  expect_equal(housepi:::DefEncoder$new(~log(alpha @ P))$x,
               "log|x|{alpha|P}")
})

test_that("Rlang basic identity assignments work as expected", {
  expect_equal(housepi:::DefEncoder$new(+3)$x, "3|N")
  expect_equal(housepi:::DefEncoder$new(-1)$x, "-1|N")

  expect_equal(housepi:::DefEncoder$new(~ alpha)$x, "alpha|A")
  expect_equal(housepi:::DefEncoder$new(~ alpha @ P)$x, "alpha|P")
  expect_equal(housepi:::DefEncoder$new(~ +alpha)$x, "alpha|A")
  expect_equal(housepi:::DefEncoder$new(~ +alpha@P)$x, "alpha|P")
  expect_equal(housepi:::DefEncoder$new(~ -alpha)$x, "minus|x|{alpha|A}")
  expect_equal(housepi:::DefEncoder$new(~ -alpha@D)$x, "minus|x|{alpha|D}")
})

test_that("Rlang arithmetic operations work as expected", {
  expect_equal(housepi:::DefEncoder$new(~ alpha + beta)$x,
               "sum|x:x2|{alpha|A}:{beta|A}")
  expect_equal(housepi:::DefEncoder$new(~ alpha @ P+ beta)$x,
               "sum|x:x2|{alpha|P}:{beta|A}")
  expect_equal(housepi:::DefEncoder$new(~ alpha + beta @ P)$x,
               "sum|x:x2|{alpha|A}:{beta|P}")
  expect_equal(housepi:::DefEncoder$new(~ alpha @ P + beta @ P)$x,
               "sum|x:x2|{alpha|P}:{beta|P}")
  expect_equal(housepi:::DefEncoder$new(~ (alpha) -(beta))$x,
               "sub|x:x2|{alpha|A}:{beta|A}")
  expect_equal(housepi:::DefEncoder$new(~ alpha * (beta))$x,
               "prod|x:x2|{alpha|A}:{beta|A}")
  expect_equal(housepi:::DefEncoder$new(~ alpha * (beta @ P))$x,
               "prod|x:x2|{alpha|A}:{beta|P}")
  expect_equal(housepi:::DefEncoder$new(~ alpha @P * beta @ P)$x,
               "prod|x:x2|{alpha|P}:{beta|P}")
})

test_that("Rlang nested expression", {
  expect_equal(housepi:::DefEncoder$new(~a + 3 * b)$x,
               "sum|x:x2|{a|A}:{prod|x:x2|{3|N}:{b|A}}")
  expect_equal(housepi:::DefEncoder$new(~a + b - c)$x,
               "sub|x:x2|{sum|x:x2|{a|A}:{b|A}}:{c|A}")
  expect_equal(housepi:::DefEncoder$new(~log ( + alpha@P  ))$x,
               "log|x|{alpha|P}")
})
