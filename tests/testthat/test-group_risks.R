library(testthat)

# a dummy dataset for the plain logic-table helpers
dummy_df <- data.frame(house_id = 1,
                       age = c(rep("child", 2), rep("adult", 2)),
                       sex = c("M", "F", "M", "M"))
# dummy group risk logic table
dummy_strata <- data.frame(isolation = c(FALSE, TRUE),
                           param_name = c(".reference", "inf_isol"))

# Model over the dummy features. Feature columns stay character on purpose:
# levels follow order of first appearance, which fixes the reference stratum.
mk_risk_model <- function(age = c("child", "child", "adult", "adult")) {
  data.frame(house_id = 1, status = 0,
             age = age,
             sex = c("M", "F", "M", "M"),
             isolation = c(FALSE, FALSE, TRUE, TRUE)) |>
    new_housepi_model("final_size")
}

####################################
## GROUP RISKS : HELPER FUNCTIONS ##
####################################

test_that("Group risks logic table computes combinations and counts", {
  tbl <- group_risk_logic_table(dummy_df, c("age", "sex"))

  # counts should sum to number of rows and include expected frequencies
  expect_equal(sum(tbl$n), nrow(dummy_df))
  expect_true(all(sort(tbl$n) == sort(c(2L, 1L, 1L))))
  # ensure columns are present
  expect_true(all(c("age", "sex", "n") %in% colnames(tbl)))
})

test_that("Associate param names decorates table and marks reference", {
  tbl <- group_risk_logic_table(dummy_df, c("age", "sex"))
  res <- associate_group_param_names(tbl, prefix = "mgs", reference = ".auto")

  expect_true("param_name" %in% colnames(res))
  # exactly one .reference row (largest group)
  expect_equal(sum(res$param_name == ".reference"), 1L)
  # other param_names should start with provided prefix
  others <- res$param_name[res$param_name != ".reference"]
  expect_true(all(grepl("^mgs_", others)))

  # absolute by default: every stratum keeps its own parameter
  abs_res <- associate_group_param_names(tbl, prefix = "mgs", reference = NULL)
  expect_false(".reference" %in% abs_res$param_name)
  expect_true(all(grepl("^mgs_", abs_res$param_name)))
})

test_that("Custom strata logic table is validated", {
  # must be a data.frame
  mk_risk_model() |>
    add_group_risk("susceptibility", strata = list(a = 1)) |>
    expect_error("data\\.frame")
  # must have param_name column
  mk_risk_model() |>
    add_group_risk("susceptibility", strata = data.frame(x = 1)) |>
    expect_error("param_name")
  # must have at least one row
  mk_risk_model() |>
    add_group_risk("susceptibility",
                   strata = data.frame(param_name = character(0))) |>
    expect_error("at least one row")
  # reserved feature column names rejected
  mk_risk_model() |>
    add_group_risk("susceptibility",
                   strata = data.frame(param_name = "a", id = 1)) |>
    expect_error("[Rr]eserved")
})

##################################
## GROUP RISKS : STANDARD CALLS ##
##################################

test_that("Group risks integrate with data", {
  hhmodel <- mk_risk_model() |>
    add_group_risk("community", by = c("age", "sex"), joint = TRUE) |>
    expect_silent()

  tidy <- get_group_risks(hhmodel)
  expect_true(is.data.frame(tidy))
  expect_true(all(tidy$group_risk == "community"))
  # joint expansion: one rule group over both features
  expect_equal(length(unique(tidy$rule_group)), 1L)
  expect_true(all(c("age", "sex") %in% colnames(tidy)))
  # absolute by default: no reference stratum
  expect_equal(sum(tidy$param_name == ".reference"), 0L)

  # the baseline the risk scales collapses to a constant
  params <- get_parameters(hhmodel)
  expect_equal(params$value[params$name == "alpha"], 1)
  expect_true(params$const[params$name == "alpha"])

  # mcmc setup is the absolute default
  params <- params[!params$const, ]
  expect_true(all(params$prior == "uniform"))
  expect_true(all(is.na(params$min)))
  expect_true(all(is.na(params$max)))
  expect_true(all(params$proposal == "rw_log"))
})

test_that("Reference turns group risks relative with its own mcmc setup", {
  hhmodel <- mk_risk_model() |>
    add_group_risk("community", by = c("age", "sex"), joint = TRUE,
                   reference = ".auto") |>
    expect_silent()

  tidy <- get_group_risks(hhmodel)
  expect_equal(sum(tidy$param_name == ".reference"), 1L)
  # the reference is the largest stratum
  expect_equal(tidy$n[tidy$param_name == ".reference"], max(tidy$n))

  # mcmc setup is default (Cauchemez et al., 2009)
  params <- get_parameters(hhmodel)
  expect_true(all(params$prior == "lognormal"))
  expect_equal(params$args[[1]], c(meanlog = 0, sdlog = 3))
  expect_true(all(params$min == 0))
  expect_true(all(params$max == 10))
  expect_true(all(params$proposal == "rw_log"))
})

test_that("Group risks hold without data and re-derive on add_data_wide", {
  # risks ingest data: a data-less model holds them, no risk params yet
  metalgear <- new_housepi_model(model = "final_size") |>
    add_parameter(import = "template") |>
    # relative: the template's baseline stays free alongside the risk
    add_group_risk("susceptibility", by = "member", reference = ".auto")
  expect_silent(metalgear)
  expect_false(any(grepl("^sus_member", metalgear$parameters$name)))

  # attaching data materialises the risk parameters and runs
  diamond <- add_data_wide(metalgear, flu_2014_matsumoto)
  full_params <- diamond$parameters$name
  expect_true(any(grepl("^sus_member", full_params)))
  expect_true(is.finite(evaluate_log_lik(diamond)))

  # replacing data re-derives: a single-level subset drops the stale risk params
  student <- flu_2014_matsumoto[flu_2014_matsumoto$member == "student", ]
  mother_base <- add_data_wide(diamond, student)
  expect_false(any(grepl("^sus_member", mother_base$parameters$name)))
  expect_true(is.finite(evaluate_log_lik(mother_base)))

  # re-attaching the full data restores the exact original parameter set
  restored <- add_data_wide(mother_base, flu_2014_matsumoto)
  expect_setequal(restored$parameters$name, full_params)
})

test_that("Group risks expand on their own covariates only", {
  df  <- sarscov2_2021_icofs |>
    dplyr::select(-start_date) |>
    dplyr::filter(house_id %in% utils::head(unique(house_id), 50))

  # a data column the risk does not read leaves the expansion untouched
  risks <- df |>
    new_housepi_model("case_ascertained") |>
    add_group_risk("susceptibility", by = "age") |>
    get_group_risks()
  expect_setequal(risks$group_risk, "susceptibility")
  expect_true("sus_age_child" %in% risks$param_name)

  # the covariate itself must be a genuine data column
  df |>
    dplyr::select(-age) |>
    new_housepi_model("case_ascertained") |>
    add_group_risk("susceptibility", by = "age") |>
    expect_error("[Uu]nknown feature.*age")

  # without data the read surface warns and keeps the table shape
  expect_warning(
    risks <- new_housepi_model(model = "final_size") |>
      add_group_risk("susceptibility", by = "member") |>
      get_group_risks(),
    regexp = "without data"
  )
  expect_equal(nrow(risks), 0L)
  expect_named(risks, c("group_risk", "param_name", "n", "rule_group"))
})

test_that("Group risks sentinels work correctly", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("final_size")

  # empty risk definition
  hhmodel |>
    add_group_risk("susceptibility", by = NA) |>
    expect_error("[Ii]nvalid risk definition")

  # wrong type for covariate column choice
  hhmodel |>
    add_group_risk("susceptibility", by = 4) |>
    expect_error("[Aa]rgument .by. must be of type character")

  # not implemented group
  hhmodel |>
    add_group_risk("cobras", by = "sex") |>
    expect_error("[Uu]nused group risk.*cobras")

  # incomplete custom prior definition
  hhmodel |>
    add_group_risk("susceptibility", by = "sex", max = 3) |>
    expect_error("[Ss]pare prior arguments")

  hhmodel |>
    add_group_risk("susceptibility", by = "test_date") |>
    expect_error("[Mm]ust be a scalar")
})

test_that("Group risks reference support behaves as expected", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("final_size")

  # an untagged custom table is absolute
  strata <- data.frame(sex = c("F", "M"),
                       param_name = c("ladies", "gentlemen"))
  hhmodel |>
    add_group_risk("susceptibility", strata = strata) |>
    expect_no_error()
  # reference has no meaning against a custom table
  hhmodel |>
    add_group_risk("susceptibility", reference = ".auto", strata = strata) |>
    expect_error("[Nn]ot used with a custom")

  # a tagged row makes the same table relative
  ocelot <- strata
  ocelot$param_name[1] <- ".reference"
  hhmodel |>
    add_group_risk("susceptibility", strata = ocelot) |>
    get_group_risks() |>
    subset(param_name == ".reference") |>
    nrow() |>
    expect_equal(1L)
})

test_that("Custom strata carry at most one reference per rule group", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("final_size")

  # two references within the same rule group are two competing scales
  liquid <- data.frame(sex = c("F", "M"),
                       param_name = c(".reference", ".reference"))
  hhmodel |>
    add_group_risk("susceptibility", strata = liquid) |>
    expect_error("[Rr]ows.*one rule group")

  # distinct rule groups each carry their own reference
  solid <- data.frame(sex = c("F", "M", NA, NA),
                      age = c(NA, NA, "child", "adult"),
                      param_name = c(".reference", "gentlemen",
                                     ".reference", "grown_ups"))
  hhmodel |>
    add_group_risk("susceptibility", strata = solid) |>
    expect_no_error()
})

test_that("Parameter options work correctly", {
  # the collapsed baseline is not one of the risk parameters
  sampled_params <- function(x) {
    params <- get_parameters(x)
    params[!params$const, ]
  }

  params <- mk_risk_model() |>
    add_group_risk("susceptibility", by = c("age", "sex"),
                   prior = "uniform", min = 0, max = 10, joint = TRUE) |>
    sampled_params()
  expect_true(all(params$prior == "uniform"))
  expect_true(all(params$min == 0))
  expect_true(all(params$max == 10))
  # a sampled risk parameter defaults to a random walk
  expect_true(all(params$proposal == "rw"))

  params <- mk_risk_model() |>
    add_group_risk("infectivity", by = c("age", "sex"),
                   prior = "gamma", shape = 2, scale = 1, joint = TRUE) |>
    sampled_params()
  expect_true(all(params$prior == "gamma"))
  expect_equal(params$args[[1]], c(shape = 2, scale = 1))
  expect_true(all(is.na(params$min)))
  expect_true(all(is.na(params$max)))
  expect_true(all(params$proposal == "rw"))
})

test_that("Load list works as expected", {
  hhmodel <- mk_risk_model() |>
    add_group_risk(
      community = list(by = c("age", "sex"), prefix = "ct",
                       prior = "uniform", min = 0, max = 10),
      susceptibility = list(by = c("age", "sex"),
                            prior = "lognormal", sdlog = 3)
    ) |>
    expect_silent()

  tidy <- get_group_risks(hhmodel)
  expect_setequal(unique(tidy$group_risk), c("community", "susceptibility"))
  expect_true(any(grepl("^ct_", tidy$param_name)))
  expect_true(any(grepl("^sus_", tidy$param_name)))
})

#############################################
## GROUP RISKS : CUSTOM STRATA DEFINITIONS ##
#############################################

test_that("Group risks deal with custom logic table correctly", {
  hhmodel <- mk_risk_model() |>
    add_group_risk(name = "infectivity", strata = dummy_strata,
                   prior = "lognormal", meanlog = 0, sdlog = 2)

  tidy <- get_group_risks(hhmodel)
  expect_true(is.data.frame(tidy))
  expect_true(all(tidy$group_risk == "infectivity"))
  # the custom strata drives the feature column and the rows
  expect_true("isolation" %in% colnames(tidy))
  expect_setequal(tidy$param_name, c(".reference", "inf_isol"))

  # parameters generated for non-reference only
  params <- get_parameters(hhmodel)
  expect_equal(params$name, "inf_isol")
  expect_equal(params$prior, "lognormal")
  expect_equal(params$args[[1]], c(meanlog = 0, sdlog = 2))
})

test_that("Group risk parameter overrides work", {
  strata <- data.frame(age = c("child", "adult", "elderly"),
                       param_name = c(".reference", "sus_old", "sus_elderly"))
  overrides <- list(
    sus_old = 1.73,
    sus_elderly = list(value = 0.5, prior = "lognormal", sdlog = 1.2,
                       min = NA_real_, max = 10,
                       proposal = "rw_log", proposal_scale = 1)
  )
  hhmodel <- mk_risk_model(c("child", "adult", "elderly", "adult")) |>
    add_group_risk(name = "susceptibility", strata = strata,
                   params = overrides, prior = "gamma", shape = 2, scale = 1)

  params <- get_parameters(hhmodel)
  expect_setequal(params$name, c("sus_old", "sus_elderly"))
  # overrides used instead of raw params
  old <- subset(params, name == "sus_old")
  expect_equal(old$value, 1.73)
  expect_true(old$const)
  elderly <- subset(params, name == "sus_elderly")
  expect_equal(elderly$prior, "lognormal")
  expect_equal(elderly$args[[1]], c(sdlog = 1.2))

  hhmodel |>
    add_group_risk(name = "susceptibility",
                   strata = data.frame(
                     sex = c("M", "F"), param_name = c("smtg", ".reference")
                   ),
                   params = list(unknown = 2.4)) |>
    expect_error("[Uu]nknown")
})

test_that("Group risk implementation preserves strata entries", {
  hhmodel <- mk_risk_model() |>
    # add a standard auto-mode risk
    add_group_risk("community", by = "age", prefix = "com") |>
    # add a custom-mode risk
    add_group_risk(name = "infectivity",
                   strata = data.frame(
                     age = c("adult", "child"),
                     param_name = c("inf_old", ".reference")
                   ))

  tidy <- get_group_risks(hhmodel)
  expect_setequal(unique(tidy$group_risk), c("community", "infectivity"))
  # strata entry is unchanged
  subset(tidy, group_risk == "infectivity")$param_name |>
    expect_setequal(c("inf_old", ".reference"))
  expect_true(any(grepl("^com_", tidy$param_name)))

  add_group_risk(hhmodel, "infectivity", by = "age") |>
    expect_warning("[Oo]verwriting.*infectivity")
  add_group_risk(hhmodel, "community", strata = dummy_strata) |>
    expect_warning("[Oo]verwriting.*community")
})

test_that("Group risk mixed auto and custom modes", {
  hhmodel <- mk_risk_model() |>
    add_group_risk("community", by = "age", prefix = "com") |>
    add_group_risk(name = "infectivity", strata = dummy_strata)

  tidy <- get_group_risks(hhmodel)
  expect_setequal(unique(tidy$group_risk), c("community", "infectivity"))
  # both modes materialise their parameters
  param_names <- get_parameters(hhmodel)$name
  expect_true("inf_isol" %in% param_names)
  expect_true(any(grepl("^com_", param_names)))

  # a feature the data does not carry fails as soon as data is attached
  data.frame(house_id = 1, status = 0, age = c("child", "adult")) |>
    new_housepi_model("final_size") |>
    add_group_risk(name = "infectivity", strata = dummy_strata) |>
    expect_error("[Uu]nknown feature.*isolation")
})

test_that("Absoluteness follows the reference across every definition mode", {
  # the default prior tracks whether a risk carries absolute values
  risk_prior <- function(...) {
    params <- get_parameters(add_group_risk(mk_risk_model(), ...))
    unique(params$prior[!params$const])
  }
  tagged <- data.frame(sex = c("F", "M"),
                       param_name = c(".reference", "gentlemen"))
  untagged <- data.frame(sex = c("F", "M"),
                         param_name = c("ladies", "gentlemen"))

  # auto mode: absolute unless a reference is resolved
  expect_equal(risk_prior("susceptibility", by = "sex"), "uniform")
  expect_equal(risk_prior("susceptibility", by = "sex",
                          reference = ".auto"), "lognormal")
  expect_equal(risk_prior("susceptibility", by = "sex",
                          reference = c(sex = "M")), "lognormal")
  # custom mode: the tag alone decides
  expect_equal(risk_prior("susceptibility", strata = untagged), "uniform")
  expect_equal(risk_prior("susceptibility", strata = tagged), "lognormal")
})

##############################
## GROUP RISKS : END-TO-END ##
##############################

test_that("An absolute risk collapses the baseline it scales", {
  # alpha is never declared: the community risk carries the whole rate
  bigboss <- flu_2014_matsumoto |>
    new_housepi_model("final_size") |>
    add_group_risk("community", by = "member") |>
    add_parameter(beta = 0.1, gamma = 1, hhsize_ref = 4)

  expect_true(is.finite(evaluate_log_lik(bigboss)))
  params <- get_parameters(bigboss)
  expect_true(params$const[params$name == "alpha"])
  expect_equal(params$value[params$name == "alpha"], 1)

  # the collapsed baseline is carried into the posterior
  draws <- bigboss |>
    run_inference(steps = 2000, thin = 2, adapt = 500, chains = 2,
                  progress = FALSE) |>
    summarise_draws()
  expect_false("alpha" %in% draws$variable)
  expect_true(any(grepl("^com_member", draws$variable)))
})

test_that("Two absolute risks on one baseline are not identifiable", {
  flu_2014_matsumoto |>
    new_housepi_model("final_size") |>
    add_group_risk("susceptibility", by = "member") |>
    add_group_risk("infectivity", by = "member") |>
    add_parameter(alpha = 0.01, gamma = 1, hhsize_ref = 4) |>
    evaluate_log_lik() |>
    expect_error("non-identifiability.*susceptibility.*infectivity")
})

#########################################
## GROUP RISKS : OBSERVATION COLLISION ##
#########################################

# Data carrying a numeric feature a risk can stratify on.
mk_obs_risk_data <- function() {
  data.frame(house_id = c(1, 1), status = c(1, 0),
             severity = c(2, 3), test_date = c(2, 6))
}

test_that("A group risk keyed on an observation is rejected", {
  # the risk feeds the transmission model the observation is generated from
  mk_obs_risk_data() |>
    new_housepi_model("final_size") |>
    add_group_risk("susceptibility", by = "severity") |>
    add_observation("severity", model = "uniform", min = 0, max = 5) |>
    expect_error("keyed on observation.*severity")

  # the reverse declaration order fails just the same
  mk_obs_risk_data() |>
    new_housepi_model("final_size") |>
    add_observation("severity", model = "uniform", min = 0, max = 5) |>
    add_group_risk("susceptibility", by = "severity") |>
    expect_error("keyed on observation.*severity")
})

test_that("A deterministic household observation cannot key a group risk", {
  # severity stands in for its data column: derived, an observation all the same
  mk_obs_risk_data() |>
    new_housepi_model("final_size") |>
    add_observation(severity = min(test_date@D), level = "house") |>
    add_group_risk("susceptibility", by = "severity") |>
    expect_error("keyed on observation.*severity")
})

test_that("A custom strata table cannot key a group risk on an observation", {
  strata <- data.frame(severity = c(2, 3),
                       param_name = c(".reference", "sus_high"))
  mk_obs_risk_data() |>
    new_housepi_model("final_size") |>
    add_observation("severity", model = "uniform", min = 0, max = 5) |>
    add_group_risk(name = "susceptibility", strata = strata) |>
    expect_error("keyed on observation.*severity")
})
