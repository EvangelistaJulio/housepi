library(testthat)

test_that("set_source_name_ assigns names to vector", {
  result <- set_source_name_(c("a", "b"), "param")
  expect_equal(names(result), c("param", "param"))
  expect_equal(unname(result), c("a", "b"))

  # NA name: unchanged
  result <- set_source_name_(c("a", "b"), NA)
  expect_null(names(result))

  # NULL input: returns NULL
  expect_null(set_source_name_(NULL, "param"))

  # empty string name: unchanged
  result <- set_source_name_(c("a"), "")
  expect_null(names(result))
})

test_that("Missing template mechanism work as expected", {
  required <- c(A = "a", B = "b", C = "c")
  defined <- c("a", "b")

  # partial overlap
  expect_equal(template_missing_find(required, defined), "c")
  # all defined
  expect_identical(template_missing_find(defined, defined), character(0))
  # no required
  expect_null(template_missing_find(character(0), c("a")))
  # no defined: returns all required sorted
  expect_equal(template_missing_find(rev(defined), character(0)), defined)

  result <- template_missing_stop(required, defined, "parameter")
  expect_true(length(result) > 0)
  expect_true(any(grepl("missing", result)))

  # no missing: empty output
  result <- template_missing_stop(c("a"), c("a"), "parameter")
  expect_equal(result, c())
})

test_that("Not Used Template mechanisms work as expected", {
  required <- c(A = "a")
  defined <- c("a", "b", "c")
  # partial overlap
  expect_equal(template_not_used_find(required, defined), c("b", "c"))
  # all used: returns empty (character(0))
  expect_length(template_not_used_find(defined, defined), 0)
  # no required: returns all defined sorted
  expect_equal(template_not_used_find(character(0), rev(defined)), defined)
  # no defined
  expect_null(template_not_used_find(c("a"), character(0)))

  # warning for unused
  expect_warning(template_not_used_warn(c("a"), defined, "venom"),
                 regexp = "unused")
  # all used: no warning
  expect_no_warning(template_not_used_warn(defined, defined, "gekko"))
})

test_that("template_overwrite_warn emits warning for overwrites", {
  expect_warning(template_overwrite_warn(c("x"), "parameter"),
                 regexp = "Overwriting")
  # empty: no warning
  expect_no_warning(template_overwrite_warn(character(0), "parameter"))
})

test_that("fmt_assert_cli formats messages correctly", {
  # with messages (non-minor): rewrites "were" to "are"
  result <- fmt_assert_cli(
    "Parameters", c("There were 2 missing parameter definitions: 'a', 'b'.")
  )
  expect_true(length(result) > 0)
  expect_true(any(grepl("Parameters", result)))
  expect_true(any(grepl("are", result)))

  # empty messages: success
  result <- fmt_assert_cli("Parameters", c())
  expect_true(any(grepl("complete", result)))

  # minor mode: rewrites "is required by" to "can be specified for"
  result <- fmt_assert_cli(
    "Parameters",
    c("There was 1 definition. It was required by augmented definitions."),
    minor = TRUE
  )
  expect_true(any(grepl("can be specified", result)))
})

test_that("Template household model branches inform correctly", {
  # print existing templates
  expect_snapshot(template_household_model())
  # non-character name
  expect_error(template_household_model(42),
               regexp = "not a character")
  # blank with no path: warning about no description
  expect_warning(template_household_model("blank"),
                 regexp = "no description")
  # unknown name
  expect_error(template_household_model("foxhound"),
               regexp = "not a known template")
})

test_that("Template name normalization", {
  model_name <- "chain binomial"
  expect_equal(spaced_template_choice("Chain-Binomial"),
               model_name)
  expect_equal(spaced_template_choice("chain_binomial"),
               model_name)
  expect_equal(spaced_template_choice(" chain binoMial "),
               model_name)
  result <- spaced_template_choice("CHAIN BINOMIAL")
  expect_equal(result, model_name)

  result <- find_closest_template_name("chain binomia")
  expect_true(any(grepl("Did you mean", result)))

  result2 <- find_closest_template_name("")
  expect_true(any(grepl("template_household_model", result2)))
})

test_that("Template extract keys from model.cpp", {
  model_dir <- get_template_dir("chain_binomial")
  model_cpp <- file.path(model_dir, "model.cpp")
  # DataKey should return empty character vector
  result <- template_extract_keys(model_cpp, "DataKey")
  expect_true(is.character(result))

  # non-existent key returns empty
  result2 <- template_extract_keys(model_cpp, "Raiden")
  expect_equal(result2, character(0))
})

test_that("Template overwrite warn handles many elements", {
  expect_warning(
    template_overwrite_warn(letters[1:10], "Miller", max_print = 3),
    regexp = "more than 3"
  )
})
