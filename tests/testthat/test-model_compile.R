library(testthat)

####  Construction  ####

test_that("assert_household_data validates inputs", {
  assert_household_data("not a frame", NA_character_) |>
    expect_error("data.frame")
  assert_household_data("test", NA_character_) |>
    expect_error("new_housepi_model")
  assert_household_data(data.frame(), NA_character_) |>
    expect_error("empty")
  # NULL data is deferred for non-test models (attach later via add_data_wide)
  assert_household_data(NULL, NA_character_) |>
    expect_null()
  assert_household_data(NULL, "test") |>
    expect_s3_class("data.frame")
})

test_that("Inconsistent model choices are caught", {
  # helper suggests correct template name
  new_housepi_model(model = "case ascertaned") |>
    expect_error("you mean.*case ascertained")

  # model name imposed redundancy
  c("case ascertained", "case-ascertained", "case-ascertained") |>
    vapply(\(x) new_housepi_model(model = x)$name, character(1)) |>
    unique() |>
    expect_identical("case ascertained")
})

test_that("Model choice helper handles malformed input", {
  # NA falls back to 'test' with a warning
  expect_warning(
    hh_model <- new_housepi_model(model = NA_character_),
    "test"
  )
  hh_model$name |> expect_identical("test")

  # non-character input is rejected
  new_housepi_model(model = 42) |> expect_error("non-empty")

  # length != 1 character input is rejected
  new_housepi_model(model = c("test", "final-size")) |>
    expect_error("non-empty")

  # '.' as a model choice is rejected with a clear message
  new_housepi_model(model = ".") |> expect_error("different folder")

  # unknown template suggests closest match (none for unrelated input)
  new_housepi_model(model = "unknown.model") |>
    expect_error("not a folder nor a template")
})

test_that("Escaped label is validated at construction", {
  # non-character input is rejected
  new_housepi_model(model = "test", escaped = 42) |>
    expect_error("non-empty")

  # length != 1 character input is rejected
  new_housepi_model(model = "test", escaped = c("healthy", "negative")) |>
    expect_error("non-empty")
})

test_that("Compile flags are validated against model source", {
  # unknown user flag
  new_housepi_model(model = "test", flags = "_UNKNOWN_FLAG") |>
    expect_error("No visible option")
  # discrete-time flag is only allowed when the model declares it
  new_housepi_model(model = "test", discr = TRUE) |>
    expect_error("discrete time")
})

test_that("as_template is refused for a built-in template", {
  expect_warning(
    hh_model <- new_housepi_model(model = "test", as_template = TRUE),
    "browse_templates"
  )
  # construction still succeeds
  hh_model$name |> expect_identical("test")
})

test_that("as_template is dropped on the 'test' fallback", {
  # a missing model falls back to 'test': only fallback warning should fire
  warns <- character()
  withCallingHandlers(
    hh_model <- new_housepi_model(as_template = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  warns |> expect_match("No model specified", all = FALSE)
  warns |> expect_no_match("Cannot store template")
  hh_model$name |> expect_identical("test")
})

test_that("A custom model directory without model.cpp is rejected", {
  path <- file.path(tempdir(), "empty_model")
  dir.create(path, showWarnings = FALSE)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)

  new_housepi_model(model = path) |> expect_error("empty_model.*was not found")
})

test_that("A build cache that cannot be created is reported", {
  path <- file.path(tempfile(), "squatted_model")
  dir.create(path, recursive = TRUE)
  # the build cache is named after the model folder, under tempdir()
  cache <- file.path(tempdir(), "squatted_model")
  file.create(cache)
  on.exit(unlink(c(dirname(path), cache), recursive = TRUE), add = TRUE)

  new_housepi_model(model = path) |> expect_error("[Ff]ailed to generate")
})

####  Runner  ####

test_that("Model has no leaks", {
  mem_used_mb <- \() gc()["Vcells", "(Mb)"]

  hhmodel <- new_housepi_model(model = "case ascertained")
  rm(hhmodel)
  invisible(gc())
  mem_before <- mem_used_mb()

  # allocation loop
  for (i in seq_len(50)) {
    hhmodel <- new_housepi_model(model = "case ascertained")
    rm(hhmodel)
    if (i %% 10 == 0) gc()
  }
  invisible(gc())
  (mem_used_mb() - mem_before) |> expect_lt(5)
})

####  Runtime compiles  ####

test_that("A model compiling at run time requires a compiler", {
  path <- tempfile("no_compiler")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  system.file("examples", "custom", "model.cpp", package = "housepi") |>
    file.copy(path)
  withr::local_options(pkgbuild.has_compiler = FALSE)

  new_housepi_model(model = path, force = TRUE) |>
    expect_error("[Nn]o.*compiler found")
})

test_that("A session template compiles user flags and warns when they change", {
  # a fresh name each run: a template loaded by an earlier run stays loaded
  path <- tempfile("session_template")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  browse_templates("case-ascertained", path, replace = TRUE)
  # dll of an earlier session, not loaded in this one: replaced, no warning
  writeLines("stale", template_dll_file(path))

  # build + load a discrete template: user paths, flags and libs reach the
  # compile command as one token each, spaces included, and the flags are
  # recorded in a sidecar
  inc <- file.path(tempdir(), "include dir")
  lib <- file.path(tempdir(), "lib dir")
  log <- suppressMessages(capture.output(
    new_housepi_model(model = path, discr = TRUE, as_template = TRUE,
                      force = TRUE, paths = inc, libs = lib),
    type = "output"
  )) |>
    paste(collapse = "\n")
  log |> expect_match(shQuote(paste0("-I", inc)), fixed = TRUE)
  log |>
    expect_match(paste0("-I", system.file("include", package = "housepi")),
                 fixed = TRUE)
  log |> expect_match("-D_HOUSEPI_DISCRETE_TIME", fixed = TRUE)
  # the link command is echoed with its quotes already removed
  log |> expect_match(paste0("-L", lib), fixed = TRUE)
  readLines(dll_flags_file(template_dll_file(path))) |>
    expect_identical("_HOUSEPI_DISCRETE_TIME")

  # a status template in the model folder is read by `import = "template"`
  dput(list(symptomatic = 0.6, asymptomatic = 0.4), file.path(path, "status.R"))
  status <- new_housepi_model(model = path, discr = TRUE, as_template = TRUE,
                              quiet = TRUE) |>
    add_status_process(import = "template") |>
    get_status_process()
  status$status |> expect_identical(c("symptomatic", "asymptomatic"))
  status$weight |> expect_equal(c(0.6, 0.4))

  # the loaded discrete dll cannot be replaced by the continuous variant
  new_housepi_model(model = path, as_template = TRUE, quiet = TRUE) |>
    expect_warning("replace existing template")

  # a dll whose flags sidecar is lost counts as unknown: rebuild is attempted
  unlink(dll_flags_file(template_dll_file(path)))
  new_housepi_model(model = path, discr = TRUE, as_template = TRUE,
                    quiet = TRUE) |>
    expect_warning("replace existing template")
})

test_that("Custom models share function registries", {
  browse_functions("growth") |> expect_error("Unknown function")

  # compile example, skipping any cached build
  hh_model <- flu_1918_frederick |>
    new_housepi_model(system.file("examples", "custom", package = "housepi"),
                      force = TRUE, quiet = TRUE)

  # no templates available for custom example
  hh_model |>
    add_latent(import = "template") |>
    expect_error("file.*not found")
  # no help page for custom example
  hh_model$help() |> expect_error("not a template")

  # Growth_ overrides only raw_function_ and raw_primitive_
  # other raw_ virtuals fall back to raise_undefined_call_err_
  x <- seq(0, 3, 0.5)
  func <- browse_functions("growth", initial = 2, rate = 0.5)
  func$func(x) |> expect_equal(2 * exp(0.5 * x))
  func$prim(x) |> expect_equal(2 * exp(0.5 * x) / 0.5)
  func$log_func(x) |> expect_error("undefined")
})
