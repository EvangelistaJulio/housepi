library(testthat)

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
