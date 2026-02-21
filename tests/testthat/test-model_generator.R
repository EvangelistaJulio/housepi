library(testthat)

test_that("Compilation flags are ok", {
  expect_false(system.file("include", package = "BH") == "")
})

test_that("Model generator inconsistent calls", {
  # blank household model object
  hhgen <- housepi:::HousepiGenerator$new("blank")

  # helper suggests correect template name
  expect_error(housepi:::HousepiGenerator$new("augment example"),
               regexp = "augmented example")

  # model name imposed redundancy
  hhgen1 <- housepi:::HousepiGenerator$new("augmented example")
  hhgen2 <- housepi:::HousepiGenerator$new("augmented_example")
  hhgen3 <- housepi:::HousepiGenerator$new("augmented-example")
  expect_equal(unique(c(hhgen1$name, hhgen2$name, hhgen3$name)),
               "augmented example")

  # inference without definitions
  expect_error(hhgen$run_inference())
})

test_that("Model has no leaks", {
  mem_used_mb <- function() gc()["Vcells", "(Mb)"]

  # warmup
  hhgen <- housepi:::HousepiGenerator$new("augmented example")
  rm(hhgen)
  invisible(gc())
  mem_before <- mem_used_mb()

  # heavy allocation loop
  for (i in cli::cli_progress_along(seq_len(200))) {
    hhgen <- housepi:::HousepiGenerator$new("augmented example")
    rm(hhgen)
    if (i %% 20 == 0) gc()
  }
  invisible(gc())
  expect_lt(mem_used_mb() - mem_before, 5)
})

test_that("Modules do not get overwritten", {
  # NOTE: use profile as model status flag
  hhgen_1 <- housepi:::HousepiGenerator$new("augmented_infection_time")
  out_1 <- hhgen_1$.__enclos_env__$private$.rcpp$export_infec_profile()
  expect_gt(length(out_1), 0)
  
  hhgen_2 <- housepi:::HousepiGenerator$new("chain_binomial")
  out_2 <- hhgen_2$.__enclos_env__$private$.rcpp$export_infec_profile()
  expect_length(out_2, 0)
  
  # reload first model base Model again
  hhgen_3 <- housepi:::HousepiGenerator$new("augmented_infection_time")
  out_3 <- hhgen_3$.__enclos_env__$private$.rcpp$export_infec_profile()
  expect_identical(out_1, out_3)
})
