library(housepi)
library(testthat)

test_that("Model has no leaks", {
  
  mem_used_mb <- function() gc()["Vcells", "(Mb)"]
  
  # warmup
  hhmodel <- new_housepi_model("augmented_example")
  rm(hhmodel)
  
  invisible(gc())
  mem_before <- mem_used_mb()
  
  for (i in cli::cli_progress_along(seq_len(100))) {
    hhmodel <- new_housepi_model("augmented_example")
    rm(hhmodel)
    if (i %% 20 == 0) gc()
  }
  invisible(gc())
  
  expect_lt(mem_used_mb() - mem_before, 10)
})

# test_that("Repeated model creation does not leak memory", {
#   
#   mem_used_mb <- function() gc()["Vcells", "(Mb)"]
#   
#   # warmup
#   hhmodel <- new_housepi_model("augmented_example")
#   hhmodel$add_data_wide(housepi::survey)
#   rm(hhmodel)
#   
#   invisible(gc())
#   mem_before <- mem_used_mb()
#   
#   for (i in cli::cli_progress_along(seq_len(100))) {
#     hhmodel <- new_housepi_model("augmented_example")
#     hhmodel$add_data_wide(housepi::survey)
#     rm(hhmodel)
#     if (i %% 20 == 0) gc()
#   }
#   invisible(gc())
#   
#   expect_lt(mem_used_mb() - mem_before, 10)
# })
