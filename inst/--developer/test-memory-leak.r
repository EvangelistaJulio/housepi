library(housepi)
library(testthat)

test_that("Repeated model creation does not leak memory", {
  
  skip_on_cran()  
  
  mem_used_mb <- function() {
    gc(full = TRUE)["Vcells", "used"] * 8 / 1024^2
  }
  
  invisible(gc(full = TRUE))
  mem_before <- mem_used_mb()
  
  for (i in seq_len(50)) {
    hhmodel <- new_housepi_model("augmented_example") |> 
      load_parameter_template() |> 
      load_augmented_template() |> 
      add_data_wide(housepi::survey)
    rm(hhmodel)
  }
  invisible(gc(full = TRUE))
  
  expect_lt(mem_used_mb() - mem_before, 5)
})
