library(testthat)

test_that("Rlang standard functions work as expected", {
  expect_true(is.na(housepi_extractor(NA)))
  expect_true(is.null(housepi_extractor(NULL)))

  housepi_extractor(log(alpha@B)) |>
    expect_error("Invalid type")
  expect_equal(housepi_extractor(log(alpha)),
               "log|x|[alpha|?]")
  expect_equal(housepi_extractor(log(alpha@P)),
               "log|x|[alpha|P]")
})

test_that("A bare name is encoded unresolved", {
  # the encoder knows no model: only the sigil settles a source here
  expect_equal(housepi_extractor(alpha), "alpha|?")
  expect_equal(housepi_extractor(alpha@P), "alpha|P")
  # numbers, numeric strings and math constants stay univocal
  expect_equal(housepi_extractor(3), "3|N")
  expect_equal(housepi_extractor("3"), "3|N")
  expect_equal(housepi_extractor(pi), "pi|N")
  expect_true(is.na(housepi_extractor(NA)))
})

test_that("Resolution rewrites the unresolved terminals only", {
  dict <- c(alpha = "P", beta = "D", both = "", solo = "A")

  expect_equal(housepi_resolve_definition("alpha|?", dict), "alpha|P")
  # an explicit sigil and a number are never touched
  expect_equal(housepi_resolve_definition("alpha|A", dict), "alpha|A")
  expect_equal(housepi_resolve_definition("3|N", dict), "3|N")
  # a name two holders claim, and a name nothing declares, keep the sentinel
  expect_equal(housepi_resolve_definition("both|?", dict), "both|?")
  expect_equal(housepi_resolve_definition("nowhere|?", dict), "nowhere|?")
  # nested blocks are rewritten in place, terminal by terminal
  expect_equal(
    housepi_resolve_definition("sum|x:x2|[alpha|?]:[prod|x:x2|[beta|?]:[3|N]]",
                               dict),
    "sum|x:x2|[alpha|P]:[prod|x:x2|[beta|D]:[3|N]]"
  )
  # an empty dictionary settles nothing
  expect_equal(housepi_resolve_definition("alpha|?", character()), "alpha|?")
  expect_true(is.na(housepi_resolve_definition(NA_character_, dict)))
})

test_that("Rlang paranthesis are interpreted correctly", {
  expect_equal(housepi_extractor(((alpha / beta))),
               housepi_extractor(alpha / beta))
  expect_equal(housepi_extractor((delta * (alpha / beta) ^ gamma)),
               housepi_extractor(delta * pow(alpha / beta, gamma)))
})

test_that("Rlang basic identity assignments work as expected", {
  expect_equal(housepi_extractor(+3), "3|N")
  expect_equal(housepi_extractor(-1), "-1|N")

  expect_equal(housepi_extractor(alpha), "alpha|?")
  expect_equal(housepi_extractor(alpha@P), "alpha|P")
  expect_equal(housepi_extractor(+alpha), "alpha|?")
  expect_equal(housepi_extractor(+alpha@P), "alpha|P")
  expect_equal(housepi_extractor(-alpha), "minus|x|[alpha|?]")
  expect_equal(housepi_extractor(-alpha@D), "minus|x|[alpha|D]")
})

test_that("Rlang arithmetic operations work as expected", {
  expect_equal(housepi_extractor(alpha + beta),
               "sum|x:x2|[alpha|?]:[beta|?]")
  expect_equal(housepi_extractor(alpha@P + beta),
               "sum|x:x2|[alpha|P]:[beta|?]")
  expect_equal(housepi_extractor(alpha + beta@P),
               "sum|x:x2|[alpha|?]:[beta|P]")
  expect_equal(housepi_extractor(alpha@P + beta@P),
               "sum|x:x2|[alpha|P]:[beta|P]")
  expect_equal(housepi_extractor((alpha) - (beta)),
               "sub|x:x2|[alpha|?]:[beta|?]")
  expect_equal(housepi_extractor(alpha * (beta)),
               "prod|x:x2|[alpha|?]:[beta|?]")
  expect_equal(housepi_extractor(alpha * (beta@P)),
               "prod|x:x2|[alpha|?]:[beta|P]")
  expect_equal(housepi_extractor(alpha@P * beta@P),
               "prod|x:x2|[alpha|P]:[beta|P]")
})

test_that("Variadic sum and prod are flattened", {
  expect_equal(housepi_extractor(2 * john * tina * cesar),
               "prod|x:x2:x3:x4|[2|N]:[john|?]:[tina|?]:[cesar|?]")
  expect_equal(housepi_extractor(a + b + c),
               "sum|x:x2:x3|[a|?]:[b|?]:[c|?]")
  expect_equal(housepi_extractor(a + b + c + d),
               "sum|x:x2:x3:x4|[a|?]:[b|?]:[c|?]:[d|?]")
  expect_equal(housepi_extractor(a + 3 * b * c),
               "sum|x:x2|[a|?]:[prod|x:x2:x3|[3|N]:[b|?]:[c|?]]")
})

test_that("Injected Date literals encode as numeric", {
  d_min <- as.Date(32)
  d_max <- as.Date(45)

  housepi_extractor(!!d_min) |>
    expect_equal(housepi_extractor(!!as.numeric(d_min)))

  # a self-contained constant expression folds to its numeric value; a call
  # over housepi nodes keeps its folded Date leaves but stays symbolic
  expect_equal(housepi_extractor(!!d_min - 5), "27|N")
  expect_equal(housepi_extractor(alpha + !!d_min),
               "sum|x:x2|[alpha|?]:[32|N]")

  m_date <- new_housepi_model(model = "test") |>
    add_latent("gray_fox", model = "uniform",
               min = !!d_min, max = !!d_max)
  m_num <- new_housepi_model(model = "test") |>
    add_latent("gray_fox", model = "uniform",
               min = !!as.numeric(d_min), max = !!as.numeric(d_max))

  expect_identical(m_date$latent$min, m_num$latent$min)
  expect_identical(m_date$latent$max, m_num$latent$max)
})

test_that("Unquoted constant calls fold to their numeric value", {
  # base-R coercions are not housepi functions: they evaluate to a constant
  expect_equal(housepi_extractor(as.Date("2021-07-01")),
               housepi_extractor(!!as.numeric(as.Date("2021-07-01"))))
  expect_equal(housepi_extractor(as.integer(5.3)), "5|N")
  expect_equal(housepi_extractor(as.numeric("3")), "3|N")

  # constant leaves fold even nested inside a symbolic housepi call
  expect_equal(housepi_extractor(diamond_dog + as.integer(5.3)),
               "sum|x:x2|[diamond_dog|?]:[5|N]")

  # a call over housepi nodes is not evaluable and stays symbolic
  expect_equal(housepi_extractor(floor(diamond_dog)),
               "floor|x|[diamond_dog|?]")
})

test_that("Rlang nested expression", {
  expect_equal(housepi_extractor(a + 3 * b),
               "sum|x:x2|[a|?]:[prod|x:x2|[3|N]:[b|?]]")
  expect_equal(housepi_extractor(a + b - c),
               "sub|x:x2|[sum|x:x2|[a|?]:[b|?]]:[c|?]")
  # UPDATE: I still don't know how to flatten differences
  expect_equal(housepi_extractor(a - b - c),
               "sub|x:x2|[sub|x:x2|[a|?]:[b|?]]:[c|?]")
  expect_equal(housepi_extractor(log(+ alpha@P)),
               "log|x|[alpha|P]")
  expect_equal(housepi_extractor(log(a * 3 * b)),
               "log|x|[prod|x:x2:x3|[a|?]:[3|N]:[b|?]]")
})

test_that("Rlang sourced unary drops the quoting backticks", {
  # backticks are R syntax for a name it cannot spell bare, not part of the key
  expect_equal(housepi_extractor(`3`@D), housepi_extractor(3@D))
  expect_equal(housepi_extractor(`3`@D), "3|D")
  expect_equal(housepi_extractor(`if`@D), "if|D")
  # an unquoted source is untouched
  expect_equal(housepi_extractor(revolver_ocelot@D), "revolver_ocelot|D")
})
