library(testthat)

# Random Metal Gear-flavoured strings used in name-validation cases.
mgs_codename <- "snake_eater"
mgs_glitch <- "fox.die"

test_that("Compilation flags are ok", {
  expect_false(system.file("include", package = "BH") == "")
})

test_that("Model generator inconsistent calls", {
  # blank household model object
  hhgen <- HousepiGenerator$new("blank")

  # helper suggests correect template name
  expect_error(HousepiGenerator$new("augment example"),
               regexp = "augmented example")

  # model name imposed redundancy
  hhgen1 <- HousepiGenerator$new("augmented example")
  hhgen2 <- HousepiGenerator$new("augmented_example")
  hhgen3 <- HousepiGenerator$new("augmented-example")
  expect_equal(unique(c(hhgen1$name, hhgen2$name, hhgen3$name)),
               "augmented example")

  # inference without definitions
  expect_error(hhgen$run_inference())
})

test_that("Model choice helper handles malformed input", {
  # NA falls back to 'blank' with a warning
  expect_warning(hhgen <- HousepiGenerator$new(NA_character_),
                 regexp = "blank")
  expect_equal(hhgen$name, "blank")
  expect_true(hhgen$is_blank())

  # non-character input is rejected
  expect_error(HousepiGenerator$new(42), regexp = "non-empty string")

  # length != 1 character input is rejected
  expect_error(HousepiGenerator$new(c("blank", mgs_codename)),
               regexp = "non-empty string")

  # '.' as a model choice is rejected with a clear message
  expect_error(HousepiGenerator$new("."), regexp = "different folder")

  # unknown template suggests closest match (none for unrelated input)
  expect_error(HousepiGenerator$new(mgs_glitch),
               regexp = "not a folder nor a template")
})

test_that("Compile flags are validated against model source", {
  # unknown user flag → stop
  expect_error(
    HousepiGenerator$new("blank", flags = "_HOUSEPI_NOT_A_REAL_FLAG"),
    regexp = "No visible option"
  )
  # discrete-time flag is only allowed when the model declares it
  expect_error(HousepiGenerator$new("blank", discr = TRUE),
               regexp = "discrete time")
})

test_that("Active fields are read-only and expose state", {
  hhgen <- HousepiGenerator$new("blank")

  # all active fields error on assignment
  ro_fields <- c("name", "raw_name", "setup", "mcmc_ptr", "source_dir",
                 "data_required", "last_flags",
                 "cols_data_wide", "info_data_wide", "infec_status_set")
  for (f in ro_fields) {
    expect_error(eval(parse(text = sprintf("hhgen$%s <- NULL", f))),
                 regexp = "read only")
  }

  # field values match expected internals
  expect_equal(hhgen$raw_name, "blank")
  expect_equal(hhgen$name, "blank")
  expect_true(dir.exists(hhgen$source_dir))
  expect_identical(hhgen$last_flags, character())
  expect_identical(hhgen$cols_data_wide, character())
  expect_null(hhgen$info_data_wide)
  expect_identical(hhgen$infec_status_set, integer())
  expect_false(is.null(hhgen$mcmc_ptr))
  expect_true(is.list(hhgen$setup))
})

test_that("get_reference points to the model source", {
  hhgen <- HousepiGenerator$new("blank")
  ref <- hhgen$get_reference()
  expect_true(file.exists(ref))
  expect_equal(basename(ref), "model.cpp")
})

test_that("data_required reflects the underlying template", {
  blank_gen <- HousepiGenerator$new("blank")
  expect_length(blank_gen$data_required, 0)

  ait_gen <- HousepiGenerator$new("augmented infection time")
  expect_true("end_followup" %in% ait_gen$data_required)
})

test_that("import_data_wide validates inputs and populates state", {
  hhgen <- HousepiGenerator$new("blank")

  expect_error(hhgen$import_data_wide("not a frame"),
               regexp = "data.frame")
  expect_error(hhgen$import_data_wide(data.frame()),
               regexp = "empty")
  expect_error(hhgen$import_data_long(list(a = 1)),
               regexp = "data.frame")

  env <- new.env(parent = emptyenv())
  utils::data("blank_house", package = "housepi", envir = env)
  hhgen$import_data_wide(env$blank_house)
  expect_setequal(hhgen$cols_data_wide, colnames(env$blank_house))
  expect_equal(hhgen$infec_status_set, 1L)
  expect_true(is.list(hhgen$info_data_wide))

  # round-trip refactor keeps the original character columns intact
  reloaded <- hhgen$refactor_data_wide(env$blank_house)
  expect_equal(reloaded$id, env$blank_house$id)
  expect_equal(reloaded$house_id, env$blank_house$house_id)
})

test_that("init_pointers is a no-op when pointers already exist", {
  hhgen <- HousepiGenerator$new("blank")
  ptr_before <- hhgen$mcmc_ptr
  expect_false(is.null(ptr_before))
  hhgen$init_pointers()
  expect_identical(hhgen$mcmc_ptr, ptr_before)
})

test_that("Model has no leaks", {
  mem_used_mb <- function() gc()["Vcells", "(Mb)"]

  hhgen <- HousepiGenerator$new("augmented example")
  rm(hhgen)
  invisible(gc())
  mem_before <- mem_used_mb()

  # heavy allocation loop
  for (i in cli::cli_progress_along(seq_len(200))) {
    hhgen <- HousepiGenerator$new("augmented example")
    rm(hhgen)
    if (i %% 20 == 0) gc()
  }
  invisible(gc())
  expect_lt(mem_used_mb() - mem_before, 5)
})
