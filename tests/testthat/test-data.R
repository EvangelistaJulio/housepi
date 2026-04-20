library(testthat)

test_that("Flu 2014 Matsumoto is correctly pre-processed", {
  hhmodel <- flu_2014_matsumoto |>
    new_housepi_model("chain_binomial_grouped")

  cases <- summarise_cases(hhmodel)
  expect_true(is.data.frame(cases))
  # expect same household count as in Endo, et al. ('19)
  # NOTE: in Table 1 missing sex was discarded, but I can't replicate that
  expect_equal(
    subset(cases, type == "household" & comment == "total")$counts,
    10486
  )
  # expect same household count as in Endo, et al. ('19)
  # expect_equal(
  #   subset(cases, type == "household" & comment == "distinct")$counts,
  #   465
  # )

  expect_true(is.data.frame(summarise_house_sizes(hhmodel)))
})
