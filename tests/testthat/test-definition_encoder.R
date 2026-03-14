library(testthat)

test_that("Rlang standard functions work as expected", {
  expect_true(is.na(housepi:::housepi_extractor(NA)))
  expect_true(is.null(housepi:::housepi_extractor(NULL)))

  expect_error(housepi:::housepi_extractor(log(alpha@B)),
               regexp = "Invalid type 'B'")
  expect_equal(housepi:::housepi_extractor(log(alpha)),
               "log|x|{alpha|A}")
  expect_equal(housepi:::housepi_extractor(log(alpha@P)),
               "log|x|{alpha|P}")
})

test_that("Rlang basic identity assignments work as expected", {
  expect_equal(housepi:::housepi_extractor(+3), "3|N")
  expect_equal(housepi:::housepi_extractor(-1), "-1|N")

  expect_equal(housepi:::housepi_extractor(alpha), "alpha|A")
  expect_equal(housepi:::housepi_extractor(alpha@P), "alpha|P")
  expect_equal(housepi:::housepi_extractor(+alpha), "alpha|A")
  expect_equal(housepi:::housepi_extractor(+alpha@P), "alpha|P")
  expect_equal(housepi:::housepi_extractor(-alpha), "minus|x|{alpha|A}")
  expect_equal(housepi:::housepi_extractor(-alpha@D), "minus|x|{alpha|D}")
})

test_that("Rlang arithmetic operations work as expected", {
  expect_equal(housepi:::housepi_extractor(alpha + beta),
               "sum|x:x2|{alpha|A}:{beta|A}")
  expect_equal(housepi:::housepi_extractor(alpha@P + beta),
               "sum|x:x2|{alpha|P}:{beta|A}")
  expect_equal(housepi:::housepi_extractor(alpha + beta@P),
               "sum|x:x2|{alpha|A}:{beta|P}")
  expect_equal(housepi:::housepi_extractor(alpha@P + beta@P),
               "sum|x:x2|{alpha|P}:{beta|P}")
  expect_equal(housepi:::housepi_extractor((alpha) - (beta)),
               "sub|x:x2|{alpha|A}:{beta|A}")
  expect_equal(housepi:::housepi_extractor(alpha * (beta)),
               "prod|x:x2|{alpha|A}:{beta|A}")
  expect_equal(housepi:::housepi_extractor(alpha * (beta@P)),
               "prod|x:x2|{alpha|A}:{beta|P}")
  expect_equal(housepi:::housepi_extractor(alpha@P * beta@P),
               "prod|x:x2|{alpha|P}:{beta|P}")
})

test_that("Variadic sum and prod are flattened", {
  expect_equal(housepi:::housepi_extractor(2 * john * tina * cesar),
               "prod|x:x2:x3:x4|{2|N}:{john|A}:{tina|A}:{cesar|A}")
  expect_equal(housepi:::housepi_extractor(a + b + c),
               "sum|x:x2:x3|{a|A}:{b|A}:{c|A}")
  expect_equal(housepi:::housepi_extractor(a + b + c + d),
               "sum|x:x2:x3:x4|{a|A}:{b|A}:{c|A}:{d|A}")
  expect_equal(housepi:::housepi_extractor(a + 3 * b * c),
               "sum|x:x2|{a|A}:{prod|x:x2:x3|{3|N}:{b|A}:{c|A}}")
})

test_that("Rlang nested expression", {
  expect_equal(housepi:::housepi_extractor(a + 3 * b),
               "sum|x:x2|{a|A}:{prod|x:x2|{3|N}:{b|A}}")
  expect_equal(housepi:::housepi_extractor(a + b - c),
               "sub|x:x2|{sum|x:x2|{a|A}:{b|A}}:{c|A}")
  # UPDATE: I still can't flatten this
  expect_equal(housepi:::housepi_extractor(a - b - c),
               "sub|x:x2|{sub|x:x2|{a|A}:{b|A}}:{c|A}")
  expect_equal(housepi:::housepi_extractor(log(+ alpha@P)),
               "log|x|{alpha|P}")
  expect_equal(housepi:::housepi_extractor(log(a * 3 * b)),
               "log|x|{prod|x:x2:x3|{a|A}:{3|N}:{b|A}}")
})
