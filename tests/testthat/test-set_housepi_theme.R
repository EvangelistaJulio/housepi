library(testthat)

####  Definitions  ####

test_that("The housepi ggplot theme is set", {
  skip_if_not_installed("ggplot2")
  set_housepi_theme()
  ggplot2::theme_get()$panel.background$fill |> expect_identical("#FEF9EC")
  ggplot2::theme_get()$strip.background$fill |> expect_identical("white")
  set_housepi_theme(colour_strips = TRUE)
  ggplot2::theme_get()$strip.background$fill |> expect_identical("#FBE09B")
})
