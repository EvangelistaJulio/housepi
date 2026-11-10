library(testthat)

####  Definitions  ####

test_that("Durations round up to the smallest time scale shown", {
  secs <- c(0, 0.2, 59.5, 61, 3599.5, 3661, 86399, 90001)
  vapply(secs, format_secs_, "") |>
    expect_identical(c("0s", "1s", "1m00s", "1m01s",
                       "1h00m", "1h02m", "1d00h", "1d02h"))
  # missing or infinite duration -> placeholder
  vapply(c(NA, Inf), format_secs_, "") |> expect_identical(c("?s", "?s"))
})

test_that("Rendered text keeps data braces literal through a later cli pass", {
  datum <- c("a{b}", "{c}")
  txt <- inline_cli_("Got {datum}.")
  cli::format_inline(txt) |> expect_identical("Got a{b} and {c}.")
  # the error names the brace-carrying option verbatim
  stop_h("Unknown entry.", find_closest_match("a{b", datum)) |>
    expect_error("[Dd]id you mean.*a\\{b\\}")
})

test_that("Find closest match suggests close options", {
  # typo within distance 1
  result <- find_closest_match("alpa", c("alpha", "beta"))
  result |> expect_match("Did you mean", all = FALSE)
  result |> expect_match("alpha", all = FALSE)

  # a short input tolerates fewer edits
  find_closest_match("sex", c("age", "beta")) |> expect_null()

  # min_dist = 0 disables suggestions
  find_closest_match("alpa", c("alpha"), min_dist = 0) |> expect_null()

  # no close match never lists the options, it returns the fallback hint
  find_closest_match("unknown", c("alpha", "beta")) |> expect_null()
  find_closest_match("unknown", c("alpha", "beta"), otherwise = "See docs.") |>
    expect_identical("See docs.")

  # empty options
  find_closest_match("x", character(0)) |> expect_identical(c())
})

####  End-to-end  ####

test_that("A cli condition from an engine callback keeps its bullets", {
  hh_model <- data.frame(hh_id = "h1", status = 0) |>
    new_housepi_model("test")
  warner <- function(x) {
    cli::cli_warn(c("Reported.", i = "A hint."))
    x
  }
  aborter <- \(x) cli::cli_abort(c("Failed.", i = "A hint."))

  # the condition reaches the caller as the callback raised it
  relayed <- hh_model |>
    add_reporting_process(dataset = warner) |>
    run_simulation() |>
    expect_warning("Reported")
  conditionMessage(relayed) |>
    expect_identical(conditionMessage(tryCatch(warner(), warning = identity)))
  relayed <- hh_model |>
    add_reporting_process(dataset = aborter) |>
    run_simulation() |>
    expect_error("Failed")
  conditionMessage(relayed) |>
    expect_identical(conditionMessage(tryCatch(aborter(), error = identity)))
})
