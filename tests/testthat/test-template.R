library(testthat)

test_that("set_source_name_ assigns names to vector", {
  result <- housepi:::set_source_name_(c("a", "b"), "param")
  expect_equal(names(result), c("param", "param"))
  expect_equal(unname(result), c("a", "b"))

  # NA name: unchanged
  result <- housepi:::set_source_name_(c("a", "b"), NA)
  expect_null(names(result))

  # NULL input: returns NULL
  expect_null(housepi:::set_source_name_(NULL, "param"))

  # empty string name: unchanged
  result <- housepi:::set_source_name_(c("a"), "")
  expect_null(names(result))
})

test_that("template_missing_find detects missing definitions", {
  # partial overlap
  expect_equal(
    housepi:::template_missing_find(c("a", "b", "c"), c("a", "b")),
    "c"
  )
  # all defined: returns empty (character(0))
  expect_length(housepi:::template_missing_find(c("a", "b"), c("a", "b")), 0)
  # no required
  expect_null(housepi:::template_missing_find(character(0), c("a")))
  # no defined: returns all required sorted
  expect_equal(
    housepi:::template_missing_find(c("b", "a"), character(0)),
    c("a", "b")
  )
})

test_that("template_not_used_find detects unused definitions", {
  # partial overlap
  expect_equal(
    housepi:::template_not_used_find(c("a"), c("a", "b", "c")),
    c("b", "c")
  )
  # all used: returns empty (character(0))
  expect_length(housepi:::template_not_used_find(c("a", "b"), c("a", "b")), 0)
  # no required: returns all defined sorted
  expect_equal(
    housepi:::template_not_used_find(character(0), c("b", "a")),
    c("a", "b")
  )
  # no defined
  expect_null(housepi:::template_not_used_find(c("a"), character(0)))
})

test_that("template_missing_stop formats error messages", {
  required <- c(a = "a", b = "b", c = "c")
  defined <- c("a", "b")
  result <- housepi:::template_missing_stop(required, defined, "parameter")
  expect_true(length(result) > 0)
  expect_true(any(grepl("missing", result)))

  # no missing: empty output
  result <- housepi:::template_missing_stop(c("a"), c("a"), "parameter")
  expect_equal(result, c())
})

test_that("template_not_used_warn emits warning for unused", {
  expect_warning(
    housepi:::template_not_used_warn(c("a"), c("a", "b"), "parameter"),
    regexp = "unused"
  )
  # all used: no warning
  expect_no_warning(
    housepi:::template_not_used_warn(c("a", "b"), c("a", "b"), "parameter")
  )
})

test_that("template_overwrite_warn emits warning for overwrites", {
  expect_warning(
    housepi:::template_overwrite_warn(c("x"), "parameter"),
    regexp = "Overwriting"
  )
  # empty: no warning
  expect_no_warning(
    housepi:::template_overwrite_warn(character(0), "parameter")
  )
})

test_that("template_blank_data_wide_ returns expected structure", {
  result <- housepi:::template_blank_data_wide_()
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true(all(c("id", "house_id", "status") %in% colnames(result)))
  expect_equal(result$status, c(1, 0))
})

test_that("fmt_assert_cli formats messages correctly", {
  # with messages (non-minor): rewrites "were" to "are"
  result <- housepi:::fmt_assert_cli(
    "Parameters", c("There were 2 missing parameter definitions: 'a', 'b'.")
  )
  expect_true(length(result) > 0)
  expect_true(any(grepl("Parameters", result)))
  expect_true(any(grepl("are", result)))

  # empty messages: success
  result <- housepi:::fmt_assert_cli("Parameters", c())
  expect_true(any(grepl("complete", result)))

  # minor mode: rewrites "is required by" to "can be specified for"
  result <- housepi:::fmt_assert_cli(
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
  expect_equal(housepi:::spaced_template_choice("Chain-Binomial"),
               model_name)
  expect_equal(housepi:::spaced_template_choice("chain_binomial"),
               model_name)
  expect_equal(housepi:::spaced_template_choice(" chain binoMial "),
               model_name)
  result <- housepi:::spaced_template_choice("CHAIN BINOMIAL")
  expect_equal(result, model_name)

  result <- housepi:::find_closest_template_name("chain binomia")
  expect_true(any(grepl("Did you mean", result)))

  result2 <- housepi:::find_closest_template_name("")
  expect_true(any(grepl("template_household_model", result2)))
})

test_that("Template extract keys from model.cpp", {
  model_dir <- housepi:::get_template_dir("chain_binomial")
  model_cpp <- file.path(model_dir, "model.cpp")
  # DataKey should return empty character vector
  result <- housepi:::template_extract_keys(model_cpp, "DataKey")
  expect_true(is.character(result))

  # non-existent key returns empty
  result2 <- housepi:::template_extract_keys(model_cpp, "Raiden")
  expect_equal(result2, character(0))
})

test_that("Template overwrite warn handles many elements", {
  expect_warning(
    housepi:::template_overwrite_warn(
      letters[1:10], "parameter", max_print = 3
    ),
    regexp = "more than 3"
  )
})
