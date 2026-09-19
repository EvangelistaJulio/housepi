library(testthat)

# Over budget: two DLL compiles, one of a session template carrying user
# paths, flags and libs, one of the custom example model.

# generator of the 'test' model, whose engine surface no model wraps
mk_generator <- local({
  list(generator = HousepiGenerator$new("test"))
})

######################################
## MODEL GENERATOR : CONSTRUCTION   ##
######################################

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

test_that("Model generator inconsistent calls are caught", {
  # helper suggests correct template name
  new_housepi_model(model = "case ascertaned") |>
    expect_error("you mean.*case ascertained")

  # model name imposed redundancy
  c("case ascertained", "case-ascertained", "case-ascertained") |>
    vapply(function(x) new_housepi_model(model = x)$name, character(1)) |>
    unique() |>
    expect_equal("case ascertained")
})

test_that("Model choice helper handles malformed input", {
  # NA falls back to 'test' with a warning
  expect_warning(
    hhmodel <- new_housepi_model(model = NA_character_),
    regexp = "test"
  )
  expect_equal(hhmodel$name, "test")

  # non-character input is rejected
  new_housepi_model(model = 42) |>
    expect_error("non-empty")

  # length != 1 character input is rejected
  new_housepi_model(model = c("test", "final-size")) |>
    expect_error("non-empty")

  # '.' as a model choice is rejected with a clear message
  new_housepi_model(model = ".") |>
    expect_error("different folder")

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
    hhmodel <- new_housepi_model(model = "test", as_template = TRUE),
    regexp = "browse_templates"
  )
  # construction still succeeds
  expect_equal(hhmodel$name, "test")
})

test_that("as_template is dropped on the 'test' fallback", {
  # a missing model falls back to 'test': only fallback warning should fire
  warns <- character()
  withCallingHandlers(
    hhmodel <- new_housepi_model(as_template = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("No model specified", warns)))
  expect_false(any(grepl("Cannot store template", warns)))
  expect_equal(hhmodel$name, "test")
})

test_that("A custom model directory without model.cpp is rejected", {
  path <- file.path(tempdir(), "empty_model")
  dir.create(path, showWarnings = FALSE)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)

  new_housepi_model(model = path) |>
    expect_error("empty_model.*was not found")
})

######################################
## MODEL GENERATOR : ENGINE SURFACE ##
######################################

test_that("Active fields expose state", {
  generator <- mk_generator$generator

  expect_equal(generator$raw_name, "test")
  expect_equal(generator$name, "test")
  expect_true(generator$is_test_model())
  expect_true(dir.exists(generator$source_dir))
  expect_identical(generator$cols_data_wide, character())
  expect_null(generator$info_data_wide)
  expect_identical(generator$infec_status_set, integer())
  expect_identical(generator$status_data_labels, character())
  expect_identical(generator$escaped, "escaped")
  expect_false(is.null(generator$mcmc_ptr))
})

test_that("Required keys reflect the underlying template", {
  generator <- mk_generator$generator
  expect_length(generator$data_required, 0L)
  expect_length(generator$group_risks_required, 0L)
  expect_length(generator$latent_required, 0L)

  generator <- HousepiGenerator$new("case ascertained")
  generator$data_required |>
    expect_setequal(c("start_date", "end_date"))
  generator$param_required |>
    expect_setequal(c("alpha", "beta", "gamma", "hhsize_ref"))
  generator$group_risks_required |>
    expect_setequal(c("community", "susceptibility", "infectivity"))
  generator$latent_required |>
    expect_setequal("infec_time")

  # every entry is tagged as coming from the transmission model
  generator$param_required |>
    names() |>
    unique() |>
    expect_equal("model")
})

test_that("Risk anchors reflect the baselines the model declares", {
  anchors <- HousepiGenerator$new("final-size")$risk_anchors
  expect_setequal(names(anchors), c("alpha", "beta"))

  # the community risk scales alpha alone, no contact involved
  expect_equal(anchors$alpha$risks, "community")
  expect_false(anchors$alpha$contacts)

  # susceptibility and infectivity scale beta, alongside contact risks
  expect_setequal(anchors$beta$risks, c("susceptibility", "infectivity"))
  expect_true(anchors$beta$contacts)

  # a model without risks declares no anchor
  expect_length(mk_generator$generator$risk_anchors, 0L)
})

test_that("Wide data is factorized and refactored by the generator", {
  generator <- HousepiGenerator$new("test")

  generator$factorize_data_wide(dummy_household)
  expect_setequal(generator$cols_data_wide, colnames(dummy_household))
  expect_equal(generator$infec_status_set, 1L)
  expect_type(generator$info_data_wide, "list")

  # round-trip refactor keeps the original character columns intact
  refactored <- generator$refactor_data_wide(dummy_household)
  expect_equal(refactored$id, dummy_household$id)
  expect_equal(refactored$hhid, dummy_household$hhid)

  # data long format
  generator$import_data_long(list(a = 1)) |>
    expect_error("data.frame")
})

test_that("init_pointers is a no-op when pointers already exist", {
  generator <- HousepiGenerator$new("test")
  ptr_before <- generator$mcmc_ptr
  expect_false(is.null(ptr_before))
  generator$init_pointers()
  expect_identical(generator$mcmc_ptr, ptr_before)
})

test_that("Model has no leaks", {
  mem_used_mb <- function() gc()["Vcells", "(Mb)"]

  generator <- HousepiGenerator$new("case ascertained")
  rm(generator)
  invisible(gc())
  mem_before <- mem_used_mb()

  # allocation loop
  for (i in seq_len(50)) {
    generator <- HousepiGenerator$new("case ascertained")
    rm(generator)
    if (i %% 10 == 0) gc()
  }
  invisible(gc())
  expect_lt(mem_used_mb() - mem_before, 5)
})

########################################
## MODEL GENERATOR : RUNTIME COMPILES ##
########################################

test_that("A session template compiles user flags and warns when they change", {
  path <- file.path(tempdir(), "session_template")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  browse_templates("case-ascertained", path, replace = TRUE)

  # build + load a discrete template: user paths, flags and libs reach the
  # compile command, and the flags are recorded in a sidecar
  log <- suppressMessages(capture.output(
    new_housepi_model(model = path, discr = TRUE, as_template = TRUE,
                      force = TRUE, paths = "/tmp/include_dir",
                      libs = "/tmp/lib_dir"),
    type = "output"
  )) |>
    paste(collapse = "\n")
  expect_match(log, "-I/tmp/include_dir", fixed = TRUE)
  expect_match(log,
               paste0("-I", system.file("include", package = "housepi")),
               fixed = TRUE)
  expect_match(log, "-D_HOUSEPI_DISCRETE_TIME", fixed = TRUE)
  expect_match(log, "-L/tmp/lib_dir", fixed = TRUE)
  expect_identical(readLines(dll_flags_file(template_dll_file(path))),
                   "_HOUSEPI_DISCRETE_TIME")

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
  browse_functions("growth") |>
    expect_error("Unknown function")

  # compile example, skipping any cached build
  hhmodel <- flu_1918_frederick |>
    new_housepi_model(system.file("examples", "custom", package = "housepi"),
                      force = TRUE, quiet = TRUE)

  # no templates available for custom example
  hhmodel |>
    add_latent(import = "template") |>
    expect_error("file.*not found")

  # Growth_ overrides only raw_function_ and raw_primitive_
  # other raw_ virtuals fall back to raise_undefined_call_err_
  x <- seq(0, 3, 0.5)
  func <- browse_functions("growth", initial = 2, rate = 0.5)
  expect_equal(func$func(x), 2 * exp(0.5 * x))
  expect_equal(func$prim(x), 2 * exp(0.5 * x) / 0.5)
  func$log_func(x) |>
    expect_error("undefined")
  # default raw_log_ratio_ returns 0
  func$log_ratio(x, x + 0.1) |>
    expect_equal(rep(0L, length(x)))
})
