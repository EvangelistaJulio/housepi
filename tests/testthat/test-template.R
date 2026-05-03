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

test_that("browse_housepi_templates branches inform correctly", {
  # print existing templates
  expect_snapshot(browse_housepi_templates())
  # non-character name
  expect_error(browse_housepi_templates(42),
               regexp = "not a character")
  # unknown name: error mentions "not a known template" and a hint
  expect_error(browse_housepi_templates("foxhound"),
               regexp = "not a known template")
  # unknown name close to a real one: hint should suggest "Did you mean"
  expect_error(browse_housepi_templates("chain binomia"),
               regexp = "Did you mean")
})

test_that("browse_housepi_templates triggers help preview when path is NA", {
  # mock the internal helper so the test never spawns a real previewer
  called_with <- NULL
  with_mocked_bindings(
    open_model_helper_ = function(path) {
      called_with <<- path
      invisible()
    },
    code = {
      browse_housepi_templates("chain_binomial")
    }
  )
  expect_true(!is.null(called_with))
  expect_true(grepl("chain_binomial$", called_with))
})

test_that("open_model_helper_ errors when info.Rd is missing", {
  empty_dir <- tempfile("ocelot_")
  dir.create(empty_dir)
  on.exit(unlink(empty_dir, recursive = TRUE), add = TRUE)
  expect_error(open_model_helper_(empty_dir),
               regexp = "Impossible to show model description")
})

test_that("browse_housepi_templates copies template to a new path", {
  dest <- tempfile("solid_snake_")
  on.exit(unlink(dest, recursive = TRUE), add = TRUE)

  # first copy: creates the model.cpp
  browse_housepi_templates("chain_binomial", path = dest)
  initial_mtime <- file.info(file.path(dest, "model.cpp"))$mtime

  expect_true(dir.exists(dest))
  expect_true(file.exists(file.path(dest, "model.cpp")))
  expect_true(file.exists(file.path(dest, "param.R")))
  expect_true(file.exists(file.path(dest, "info.Rd")))
  # chain_binomial has no augmented file: nothing should be copied
  expect_false(file.exists(file.path(dest, "augm.R")))

  # second copy without replace: warns and does NOT overwrite
  expect_warning(
    browse_housepi_templates("chain_binomial", path = dest),
    regexp = "already exists"
  )
  expect_equal(file.info(file.path(dest, "model.cpp"))$mtime, initial_mtime)

  # with replace = TRUE: overwrites silently
  Sys.sleep(1.1)  # ensure mtime resolution captures the rewrite
  expect_no_warning(
    browse_housepi_templates("chain_binomial", path = dest, replace = TRUE)
  )
  expect_gt(
    as.numeric(file.info(file.path(dest, "model.cpp"))$mtime),
    as.numeric(initial_mtime)
  )
})

test_that("browse_housepi_templates rejects '.' as destination path", {
  sandbox <- tempfile("mantis_")
  dir.create(sandbox)
  on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)
  withr::with_dir(sandbox, {
    expect_warning(
      browse_housepi_templates("chain_binomial", path = "."),
      regexp = "current directory"
    )
    # helper redirects "." to a folder named after the template
    expect_true(dir.exists("chain_binomial"))
  })
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
  expect_true(any(grepl("browse_housepi_templates", result2)))
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
