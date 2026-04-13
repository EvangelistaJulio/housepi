library(testthat)

aux_add <- function(self, by = NA_character_,
                    reference = NULL, prefix = "cr",
                    logic_table = NULL, params = NULL,
                    prior = NA_character_, ...,
                    min = NA_real_, max = NA_real_,
                    proposal = NA_character_,
                    proposal_scale = NA_real_,
                    .warn_over = TRUE) {
  self$add_quo(by = by, prefix = prefix,
               reference = reference,
               logic_table = logic_table,
               params = params,
               prior = prior, rlang::enquos(...),
               min = rlang::enquo(min),
               max = rlang::enquo(max),
               proposal = proposal,
               proposal_scale = proposal_scale,
               .warn_over = .warn_over)
  invisible(self)
}

# a dummy dataset and relative information for factorization
dummy_df <- data.frame(age_group = c("child", "child", "adult", "child"),
                       sex = c("M", "F", "M", "M"),
                       stringsAsFactors = FALSE)
dummy_factor_info <- list(age_group = list(levels = c("child", "adult")),
                          sex = list(levels = c("M", "F")))

test_that("Contact rate custom mode stores correctly", {
  # logic table
  param_name <- c("cr_a2c", "cr_c2a", "cr_c2c")
  lt <- data.frame(from_age = c(rep("adult", 2), rep("child", 2)),
                   to_age = rep(c("adult", "child"), 2),
                   param_name = c(".reference", param_name))
  # custom fine grained parameter override
  overrides <- list(
    list(value = 1, prior = "uniform", min = 1, max = 5,
         proposal = "normal", proposal_scale = 0.5),
    2.34
  )
  fancy_param_1 <- param_name[1] # param to override as variable
  fancy_param_2 <- param_name[2] # param to override as constant
  default_param <- param_name[3] # default param definitions
  names(overrides) <- c(fancy_param_1, fancy_param_2)

  cr <- HousepiContactRate$new()

  expect_error(aux_add(cr, by = "age", logic_table = lt),
               regexp = "Impossible")

  aux_add(cr, logic_table = lt, prior = "lognormal", sdlog = 1.8,
          params = overrides)

  # definitions stored with mode sentinel
  expect_equal(cr$defs$mode, "custom")

  # logic table stored
  expect_true(is.data.frame(cr$logic_table))
  expect_equal(nrow(cr$logic_table), 4)

  # parameters generated for non-reference
  expect_length(cr$parameters, 3)
  expect_true(all(param_name %in% names(cr$parameters)))

  # default parameter definitions
  expect_equal(cr$parameters[[default_param]]$prior, "lognormal")
  expect_equal(cr$parameters[[default_param]]$sdlog, "1.8|N")
  # overrides as variable
  expect_equal(cr$parameters[[fancy_param_1]]$prior, "uniform")
  expect_equal(cr$parameters[[fancy_param_1]]$min, 1)
  # overrides as constant
  expect_equal(cr$parameters[[fancy_param_2]], 2.34)

  # implement should not clear the logic table
  cr$implement(dummy_df, dummy_factor_info)
  expect_true(is.data.frame(cr$logic_table))
  expect_equal(cr$logic_table$param_name, lt$param_name)

  # tidy table
  tidy <- cr$tidy_table
  expect_true(is.data.frame(tidy))
  expect_equal(nrow(tidy), 4)
})
