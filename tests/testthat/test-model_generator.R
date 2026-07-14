library(testthat)

test_that("Compilation flags are ok", {
  expect_false(system.file("include", package = "BH") == "")
})

test_that("assert_household_data validates inputs", {
  assert_household_data("not a frame", NA_character_) |>
    expect_error("data.frame")
  assert_household_data("test", NA_character_) |>
    expect_error("new_housepi_model")
  assert_household_data(data.frame(), NA_character_) |>
    expect_error("empty")
  assert_household_data(NULL, NA_character_) |>
    expect_error("Missing argument")
  assert_household_data(NULL, "test") |>
    expect_s3_class("data.frame")
})

test_that("Model generator inconsistent calls are caught", {
  # household model object
  hhgen <- HousepiGenerator$new("test")

  # helper suggests correect template name
  HousepiGenerator$new("case ascertaned") |>
    expect_error("you mean.*case ascertained")

  # model name imposed redundancy
  name1 <- HousepiGenerator$new("case ascertained")$name
  name2 <- HousepiGenerator$new("case_ascertained")$name
  name3 <- HousepiGenerator$new("case-ascertained")$name
  unique(c(name1, name2, name3)) |>
    expect_equal("case ascertained")

  # inference without definitions
  run_inference(hhgen) |>
    expect_error("HousepiModel")
})

test_that("Model choice helper handles malformed input", {
  # NA falls back to 'test' with a warning
  expect_warning(hhgen <- HousepiGenerator$new(NA_character_),
                 regexp = "test")
  expect_equal(hhgen$name, "test")
  expect_true(hhgen$is_test_model())

  # non-character input is rejected
  HousepiGenerator$new(42) |>
    expect_error("non-empty")

  # length != 1 character input is rejected
  HousepiGenerator$new(c("test", "snake_eater")) |>
    expect_error("non-empty")

  # '.' as a model choice is rejected with a clear message
  HousepiGenerator$new(".") |>
    expect_error("different folder")

  # unknown template suggests closest match (none for unrelated input)
  HousepiGenerator$new("fox.die") |>
    expect_error("not a folder nor a template")
})

test_that("Compile flags are validated against model source", {
  # unknown user flag → stop
  HousepiGenerator$new("test", flags = "_HAIDARAAA") |>
    expect_error("No visible option")
  # discrete-time flag is only allowed when the model declares it
  HousepiGenerator$new("test", discr = TRUE) |>
    expect_error("discrete time")
})

test_that("Active fields are read-only and expose state", {
  hhgen <- HousepiGenerator$new("test")

  # all active fields error on assignment
  hhgen |>
    helper_test_read_only(
      "name", "raw_name", "setup", "mcmc_ptr", "source_dir",
      "data_required", "last_flags",
      "cols_data_wide", "info_data_wide", "infec_status_set"
    )

  # field values match expected internals
  expect_equal(hhgen$raw_name, "test")
  expect_equal(hhgen$name, "test")
  expect_true(dir.exists(hhgen$source_dir))
  expect_identical(hhgen$last_flags, character())
  expect_identical(hhgen$cols_data_wide, character())
  expect_null(hhgen$info_data_wide)
  expect_identical(hhgen$infec_status_set, integer())
  expect_false(is.null(hhgen$mcmc_ptr))
  expect_true(is.list(hhgen$setup))
})

test_that("get_reference points to the model source", {
  hhgen <- HousepiGenerator$new("test")
  ref <- hhgen$get_reference()
  expect_true(file.exists(ref))
  expect_equal(basename(ref), "model.cpp")
})

test_that("data_required reflects the underlying template", {
  test_gen <- HousepiGenerator$new("test")
  test_gen$data_required |>
    expect_length(0L)

  ait_gen <- HousepiGenerator$new("case ascertained")
  expect_true("start_date" %in% ait_gen$data_required)
  expect_true("end_date" %in% ait_gen$data_required)
})

test_that("factorize_data_wide behaves as expected", {
  hhgen <- HousepiGenerator$new("test")

  env <- new.env(parent = emptyenv())
  utils::data("dummy_household", package = "housepi", envir = env)
  hhgen$factorize_data_wide(env$dummy_household) |>
    expect_silent()
  hhgen$cols_data_wide |>
    expect_setequal(colnames(env$dummy_household))
  hhgen$infec_status_set |>
    expect_equal(1L)
  is.list(hhgen$info_data_wide) |>
    expect_true()

  # round-trip refactor keeps the original character columns intact
  reloaded <- hhgen$refactor_data_wide(env$dummy_household)
  reloaded$id |>
    expect_equal(env$dummy_household$id)
  reloaded$house_id |>
    expect_equal(env$dummy_household$house_id)

  # data long format
  hhgen$import_data_long(list(a = 1)) |>
    expect_error("data.frame")
})

test_that("user paths, flags and libs reach the compile command", {
  log <- suppressMessages(capture.output(
    hhgen <- HousepiGenerator$new("case ascertained",
                                  discr = TRUE, force = TRUE,
                                  paths = "/tmp/diamond_dogs",
                                  libs  = "/tmp/mother_base"),
    type = "output"
  )) |>
    paste(collapse = "\n")

  # match across compilation commands
  expect_match(log, "-I/tmp/diamond_dogs", fixed = TRUE)
  expect_match(log,
               paste0("-I", system.file("include", package = "housepi")),
               fixed = TRUE)
  expect_match(log, "-D_HOUSEPI_DISCRETE_TIME", fixed = TRUE)
  expect_match(log, "-L/tmp/mother_base", fixed = TRUE)
})

test_that("as_template is refused for a built-in template", {
  expect_warning(hhgen <- HousepiGenerator$new("test", as_template = TRUE),
                 regexp = "browse_housepi_templates")
  # construction still succeeds
  expect_false(is.null(hhgen$mcmc_ptr))
})

test_that("as_template is dropped on the 'test' fallback", {
  # a missing model falls back to 'test': only fallback warning should fire
  warns <- character()
  withCallingHandlers(
    hhgen <- HousepiGenerator$new(as_template = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("No model specified", warns)))
  expect_false(any(grepl("Cannot store template", warns)))
  expect_false(is.null(hhgen$mcmc_ptr))
})

test_that("a session-loaded template warns when compile flags change", {
  path <- file.path(tempdir(), "shadow_moses")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  browse_housepi_templates("case_ascertained", path, replace = TRUE)

  # build + load a discrete template, recording its flags in a sidecar
  HousepiGenerator$new(path, discr = TRUE, as_template = TRUE, quiet = TRUE)
  expect_identical(readLines(dll_flags_file(template_dll_file(path))),
                   "_HOUSEPI_DISCRETE_TIME")

  # the loaded discrete dll cannot be replaced by the continuous variant:
  # the warning must fire in this direction too (was previously silent)
  HousepiGenerator$new(path, as_template = TRUE, quiet = TRUE) |>
    expect_warning("replace existing template")
})

test_that("custom model directory without model.cpp is rejected", {
  path <- file.path(tempdir(), "outer_heaven")
  dir.create(path, showWarnings = FALSE)

  HousepiGenerator$new(path) |>
    expect_error("outer_heaven.*was not found")
})

test_that("an error is raised when tracked referen is lost", {
  dest <- tempfile("solid_snake_")
  on.exit(unlink(dest, recursive = TRUE), add = TRUE)

  browse_housepi_templates("final_size", path = dest)
  hhgen <- HousepiGenerator$new(dest, quiet = TRUE)
  unlink(hhgen$get_reference())

  hhgen$data_required |>
    expect_error("[Ii]mpossible to find reference file")
})


test_that("init_pointers is a no-op when pointers already exist", {
  hhgen <- HousepiGenerator$new("test")
  ptr_before <- hhgen$mcmc_ptr
  expect_false(is.null(ptr_before))
  hhgen$init_pointers()
  expect_identical(hhgen$mcmc_ptr, ptr_before)
})

test_that("Model has no leaks", {
  mem_used_mb <- function() gc()["Vcells", "(Mb)"]

  hhgen <- HousepiGenerator$new("case ascertained")
  rm(hhgen)
  invisible(gc())
  mem_before <- mem_used_mb()

  # heavy allocation loop
  for (i in cli::cli_progress_along(seq_len(200))) {
    hhgen <- HousepiGenerator$new("case ascertained")
    rm(hhgen)
    if (i %% 20 == 0) gc()
  }
  invisible(gc())
  expect_lt(mem_used_mb() - mem_before, 5)
})
