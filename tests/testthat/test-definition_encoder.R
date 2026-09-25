library(testthat)

test_that("Rlang standard functions work as expected", {
  housepi_extractor(NULL) |> expect_null()
  housepi_extractor(NA) |> expect_all_na()

  housepi_extractor(log(alpha@B)) |>
    expect_error("Invalid type")
  housepi_extractor(log(alpha)) |>
    expect_identical("log|x|[alpha|?]")
  housepi_extractor(log(alpha@P)) |>
    expect_identical("log|x|[alpha|P]")
})

test_that("Invalid definition values are refused at encoding", {
  # the reserved columns are no data source
  housepi_extractor(id@D) |>
    expect_error("[Cc]annot.*id.*as data.*reserved")
  housepi_extractor(hhid@D) |>
    expect_error("[Cc]annot.*hhid.*as data.*reserved")
  housepi_extractor(status@D) |>
    expect_error("[Cc]annot.*status.*as data.*reserved")
  # a vector is no single source
  housepi_extractor(!!c(1, 2)) |>
    expect_error("[Nn]ot a valid source")
  housepi_extractor(!!sum) |>
    expect_error("[Nn]ot a valid source")
  # a call already carries its arguments
  housepi_combine_to_call("log|x|[a|A]", list(x = "1|N")) |>
    expect_error("[Ii]mpossible to associate arguments to a function")
})

test_that("A bare name is encoded unresolved", {
  # the encoder knows no model: only the sigil settles a source here, while
  # numbers, numeric strings and math constants stay univocal
  housepi_extractor(3) |>
    expect_identical("3|N")
  housepi_extractor("3") |>
    expect_identical("3|N")
  housepi_extractor(pi) |>
    expect_identical("pi|N")
})

test_that("Resolution rewrites the unresolved terminals only", {
  dict <- c(alpha = "P", beta = "D", both = "", solo = "A")

  housepi_resolve_definition("alpha|?", dict) |>
    expect_identical("alpha|P")
  # an explicit sigil and a number are never touched
  housepi_resolve_definition("alpha|A", dict) |>
    expect_identical("alpha|A")
  housepi_resolve_definition("3|N", dict) |>
    expect_identical("3|N")
  # a name two holders claim, and a name nothing declares, keep the sentinel
  housepi_resolve_definition("both|?", dict) |>
    expect_identical("both|?")
  housepi_resolve_definition("nowhere|?", dict) |>
    expect_identical("nowhere|?")
  # nested blocks are rewritten in place, terminal by terminal
  housepi_resolve_definition("sum|x:x2|[alpha|?]:[prod|x:x2|[beta|?]:[3|N]]",
                             dict) |>
    expect_identical("sum|x:x2|[alpha|P]:[prod|x:x2|[beta|D]:[3|N]]")
  # an empty dictionary settles nothing
  housepi_resolve_definition("alpha|?", character()) |>
    expect_identical("alpha|?")
  is.na(housepi_resolve_definition(NA_character_, dict)) |>
    expect_true()
})

test_that("Rlang parentheses are interpreted correctly", {
  housepi_extractor(((alpha / beta))) |>
    expect_equal(housepi_extractor(alpha / beta))
  housepi_extractor((delta * (alpha / beta) ^ gamma)) |>
    expect_equal(housepi_extractor(delta * pow(alpha / beta, gamma)))
})

test_that("Rlang basic identity assignments work as expected", {
  housepi_extractor(+3) |> expect_identical("3|N")
  housepi_extractor(-1) |> expect_identical("-1|N")

  housepi_extractor(alpha) |> expect_identical("alpha|?")
  housepi_extractor(alpha@P) |> expect_identical("alpha|P")
  housepi_extractor(+alpha) |> expect_identical("alpha|?")
  housepi_extractor(+alpha@P) |> expect_identical("alpha|P")
  housepi_extractor(-alpha) |> expect_identical("minus|x|[alpha|?]")
  housepi_extractor(-alpha@D) |> expect_identical("minus|x|[alpha|D]")
})

test_that("Rlang arithmetic operations work as expected", {
  housepi_extractor(alpha + beta) |>
    expect_identical("+|x:x2|[alpha|?]:[beta|?]")
  housepi_extractor(alpha@P + beta) |>
    expect_identical("+|x:x2|[alpha|P]:[beta|?]")
  housepi_extractor(alpha + beta@P) |>
    expect_identical("+|x:x2|[alpha|?]:[beta|P]")
  housepi_extractor(alpha@P + beta@P) |>
    expect_identical("+|x:x2|[alpha|P]:[beta|P]")
  housepi_extractor((alpha) - (beta)) |>
    expect_identical("+|x:x2|[alpha|?]:[minus|x|[beta|?]]")
  housepi_extractor(alpha * (beta)) |>
    expect_identical("*|x:x2|[alpha|?]:[beta|?]")
  housepi_extractor(alpha * (beta@P)) |>
    expect_identical("*|x:x2|[alpha|?]:[beta|P]")
  housepi_extractor(alpha@P * beta@P) |>
    expect_identical("*|x:x2|[alpha|P]:[beta|P]")
})

test_that("Operator chains are flattened, explicit reductions are not", {
  housepi_extractor(2 * rho * rho_2 * rho_3) |>
    expect_identical("*|x:x2:x3:x4|[2|N]:[rho|?]:[rho_2|?]:[rho_3|?]")
  housepi_extractor(a + b + c) |>
    expect_identical("+|x:x2:x3|[a|?]:[b|?]:[c|?]")
  housepi_extractor(a + b + c + d) |>
    expect_identical("+|x:x2:x3:x4|[a|?]:[b|?]:[c|?]:[d|?]")
  housepi_extractor(a + 3 * b * c) |>
    expect_identical("+|x:x2|[a|?]:[*|x:x2:x3|[3|N]:[b|?]:[c|?]]")
  # an explicit reduction keeps its arguments as written
  housepi_extractor(sum(a + b)) |>
    expect_identical("sum|x|[+|x:x2|[a|?]:[b|?]]")
  housepi_extractor(prod(a, b * c)) |>
    expect_identical("prod|x:x2|[a|?]:[*|x:x2|[b|?]:[c|?]]")
})

test_that("Injected Date literals encode as numeric", {
  d_min <- as.Date(32)
  d_max <- as.Date(45)

  housepi_extractor(!!d_min) |>
    expect_identical(housepi_extractor(!!as.numeric(d_min)))

  # a fully numeric constant expression folds to its numeric value
  housepi_extractor(!!d_min - 5) |>
    expect_identical("27|N")
  # a Date folds to numeric
  housepi_extractor(alpha + !!d_min) |>
    expect_identical("+|x:x2|[alpha|?]:[32|N]")

  m_date <- new_housepi_model(model = "test") |>
    add_latent("hidden", model = "uniform", min = !!d_min, max = !!d_max)
  m_num <- new_housepi_model(model = "test") |>
    add_latent("hidden", model = "uniform",
               min = !!as.numeric(d_min), max = !!as.numeric(d_max))
  m_date$latent$min |> expect_identical(m_num$latent$min)
  m_date$latent$max |> expect_identical(m_num$latent$max)
})

test_that("Unquoted constant calls fold to their numeric value", {
  # base-R coercions are not housepi functions: they evaluate to a constant
  housepi_extractor(as.Date("2021-07-01")) |>
    expect_identical(housepi_extractor(!!as.numeric(as.Date("2021-07-01"))))
  housepi_extractor(as.integer(5.3)) |>
    expect_identical("5|N")
  housepi_extractor(as.numeric("3")) |>
    expect_identical("3|N")

  # constant leaves fold even nested inside a symbolic housepi call
  housepi_extractor(rho + as.integer(5.3)) |>
    expect_identical("+|x:x2|[rho|?]:[5|N]")

  # a call over housepi nodes is not evaluable and stays symbolic
  housepi_extractor(floor(rho)) |>
    expect_identical("floor|x|[rho|?]")
})

test_that("Rlang nested expression", {
  housepi_extractor(a + 3 * b) |>
    expect_identical("+|x:x2|[a|?]:[*|x:x2|[3|N]:[b|?]]")
  housepi_extractor(a + b - c) |>
    expect_identical("+|x:x2:x3|[a|?]:[b|?]:[minus|x|[c|?]]")
  housepi_extractor(a - b - c) |>
    expect_identical("+|x:x2:x3|[a|?]:[minus|x|[b|?]]:[minus|x|[c|?]]")
  housepi_extractor(log(+ alpha@P)) |>
    expect_identical("log|x|[alpha|P]")
  housepi_extractor(log(a * 3 * b)) |>
    expect_identical("log|x|[*|x:x2:x3|[a|?]:[3|N]:[b|?]]")
})

test_that("Rlang sourced unary drops the quoting backticks", {
  housepi_extractor(`3`@D) |>
    expect_identical(housepi_extractor(3@D))
  housepi_extractor(`3`@D) |>
    expect_identical("3|D")
  housepi_extractor(`if`@D) |>
    expect_identical("if|D")
  housepi_extractor(datum@D) |>
    expect_identical("datum|D")
})
