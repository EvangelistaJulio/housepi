library(testthat)

# a dummy dataset for the plain logic-table helpers
dummy_df <- data.frame(house_id = 1,
                       age = c(rep("child", 2), rep("adult", 2)),
                       sex = c("M", "F", "M", "M"))
# dummy group risk logic table
dummy_strata <- data.frame(isolation = c(FALSE, TRUE),
                           param_name = c(".reference", "inf_isol"))

# Model over the dummy covariates. Covariate columns stay character on purpose:
# levels follow order of first appearance, which fixes the reference stratum.
mk_risk_model <- function(age = c("child", "child", "adult", "adult")) {
  data.frame(house_id = 1, status = 0,
             age = age,
             sex = c("M", "F", "M", "M"),
             isolation = c(FALSE, FALSE, TRUE, TRUE)) |>
    new_housepi_model("final_size")
}

# a small "member"-stratified dataset standing in for flu_2014_matsumoto:
# five households, "student" is the only role present in every household.
flu_small <- data.frame(
  house_id = c(1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5),
  status = c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1),
  member = c("student", "other", "student", "sibling",
             "student", "mother", "father",
             "student", "other",
             "student", "mother", "father")
)

# shared base model for icofs-backed validation tests: add_group_risk()
# copies before mutating, so reusing this fixture across tests is safe.
icofs_model <- sarscov2_2021_icofs |> new_housepi_model("final_size")

# a simulation table giving every parameter of the model a value
pin_all_params <- function(hhmodel, value = 0.5) {
  names <- get_parameters(hhmodel)$name
  as.data.frame(stats::setNames(as.list(rep(value, length(names))), names))
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
  # reserved covariate column names rejected
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
  # joint expansion: one rule group over both covariates
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
  diamond <- add_data_wide(metalgear, flu_small)
  full_params <- diamond$parameters$name
  expect_true(any(grepl("^sus_member", full_params)))
  expect_true(is.finite(evaluate_log_lik(diamond)))

  # replacing data re-derives: a single-level subset drops the stale risk params
  student <- flu_small[flu_small$member == "student", ]
  mother_base <- add_data_wide(diamond, student)
  expect_false(any(grepl("^sus_member", mother_base$parameters$name)))
  expect_true(is.finite(evaluate_log_lik(mother_base)))

  # re-attaching the full data restores the exact original parameter set
  restored <- add_data_wide(mother_base, flu_small)
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
    expect_error("[Uu]nknown covariate.*age")

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
  hhmodel <- icofs_model

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
    expect_error("[Mm]ust be a unidimensional type")
})

test_that("Group risks reference support behaves as expected", {
  hhmodel <- icofs_model

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
  hhmodel <- icofs_model

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
  # the custom strata drives the covariate column and the rows
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

  # a covariate the data does not carry fails as soon as data is attached
  data.frame(house_id = 1, status = 0, age = c("child", "adult")) |>
    new_housepi_model("final_size") |>
    add_group_risk(name = "infectivity", strata = dummy_strata) |>
    expect_error("[Uu]nknown covariate.*isolation")
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

###############################
## GROUP RISKS : STRATA PINS ##
###############################

# Custom strata over `flu_small`, pinning one stratum and fixing another.
pin_strata <- function(value = c(NA, 1.2, 0.5), const = c(NA, NA, TRUE)) {
  data.frame(member = c("student", "mother", "father"),
             param_name = c(".reference", "sus_mum", "sus_dad"),
             value = value, const = const)
}

test_that("A strata value pins the stratum and keeps its prior", {
  hhmodel <- flu_small |>
    new_housepi_model("final_size") |>
    add_group_risk("susceptibility", strata = pin_strata(),
                   prior = "lognormal", sdlog = 3)

  params <- get_parameters(hhmodel)
  mum <- subset(params, name == "sus_mum")
  expect_equal(mum$value, 1.2)
  expect_false(mum$const)
  expect_equal(mum$prior, "lognormal")

  # a fixed stratum is a constant parameter: it carries no prior
  dad <- subset(params, name == "sus_dad")
  expect_equal(dad$value, 0.5)
  expect_true(dad$const)
  expect_true(is.na(dad$prior))
})

test_that("Pin columns surface on the read table only when a pin exists", {
  hhmodel <- flu_small |>
    new_housepi_model("final_size") |>
    add_group_risk("susceptibility", strata = pin_strata(),
                   prior = "lognormal", sdlog = 3)

  tidy <- get_group_risks(hhmodel)
  expect_true(all(c("value", "const") %in% colnames(tidy)))
  # padding keeps the columns typed across sub-tables
  expect_type(tidy$value, "double")
  expect_type(tidy$const, "logical")
  expect_equal(tidy$value[tidy$param_name == "sus_mum"], 1.2)
  expect_true(is.na(tidy$value[tidy$param_name == ".reference"]))

  # an unpinned risk keeps the table it had
  flu_small |>
    new_housepi_model("final_size") |>
    add_group_risk("susceptibility", by = "member") |>
    get_group_risks() |>
    colnames() |>
    intersect(c("value", "const")) |>
    expect_length(0)
})

test_that("Invalid strata pins are rejected", {
  mk <- function(strata) {
    flu_small |>
      new_housepi_model("final_size") |>
      add_group_risk("susceptibility", strata = strata,
                     prior = "lognormal", sdlog = 3)
  }

  # a fixed stratum needs the value it is fixed at
  mk(pin_strata(value = c(NA, 1.2, NA), const = c(NA, NA, TRUE))) |>
    expect_error("sus_dad")
  # the reference is 1 by definition
  mk(pin_strata(value = c(0.9, 1.2, 0.5))) |>
    expect_error("[Rr]eference")
  mk(pin_strata(const = c(TRUE, NA, TRUE))) |>
    expect_error("[Rr]eference")
  # the columns carry a fixed type
  mk(transform(pin_strata(), value = as.character(value))) |>
    expect_error("numeric")
  mk(transform(pin_strata(), const = as.character(const))) |>
    expect_error("logical")

  # a stratum spread over two rules carries a single pin
  data.frame(member = c("student", "mother", "mother"),
             param_name = c(".reference", "sus_mum", "sus_mum"),
             value = c(NA, 1.2, 0.4)) |>
    mk() |>
    get_parameters() |>
    expect_error("[Cc]onflicting")
})

test_that("A pinned stratum takes no entry in params", {
  flu_small |>
    new_housepi_model("final_size") |>
    add_group_risk("susceptibility", strata = pin_strata(),
                   params = list(sus_mum = 2), prior = "lognormal",
                   sdlog = 3) |>
    expect_error("sus_mum")
})

#################################
## GROUP RISKS : STATUS STRATA ##
#################################

# two observed infected statuses, plus a data column duplicating them on the
# infected members only: the covariate a status stratification replaces.
status <- c("escaped", "symptomatic", "asymptomatic", "symptomatic",
            "escaped", "asymptomatic", "symptomatic", "escaped",
            "asymptomatic", "symptomatic", "asymptomatic", "escaped")
snake_data <- data.frame(
  house_id = rep(1:4, each = 3),
  status = status,
  shadow = ifelse(status == "escaped", NA, status)
)

test_that("Status is accepted as a stratification axis", {
  # both gates let 'status' through for a group risk
  snake_model <- snake_data |> new_housepi_model("final_size")
  snake_model |>
    add_group_risk("infectivity", by = "status") |>
    expect_silent()
  snake_model |>
    add_group_risk("infectivity",
                   strata = data.frame(status = "symptomatic",
                                       param_name = "inf_ocelot")) |>
    expect_silent()

  # ... and neither exposes it as a required data column
  module <- snake_model |>
    add_group_risk(infectivity = list(
      by = "status", reference = list(status = "symptomatic")
    )) |>
    (\(m) m$.__enclos_env__$private$group_risks_)()
  expect_false("status" %in% module$data_required)
  expect_false("status" %in% module$covariate_names)
  expect_true(module$stratifies_on_status)

  # one parameter per non-reference stratum, named after the label
  expect_equal(unname(module$param_required), "inf_status_asymptomatic")
})

test_that("Only the risks the model declares may be status stratified", {
  snake_model <- snake_data |> new_housepi_model("final_size")
  # infectivity is read after the member's status is settled, the others are
  # read on members that carry no infected status
  snake_model |>
    add_group_risk("susceptibility", by = "status") |>
    expect_error("cannot be stratified by the infected status")
  snake_model |>
    add_group_risk("community", by = "status") |>
    expect_error("cannot be stratified by the infected status")

  # contact risks keep no status axis at all
  snake_model |>
    add_contact_risk(by = "status") |>
    expect_error("reserved")

  # a model that declares no status level cannot stratify on it. The risk
  # waits for the declaration, so decorator order is free, and the wall reports
  # that it never came
  mk_risk_model() |>
    add_group_risk("infectivity", by = "status") |>
    expect_silent()
  mk_risk_model() |>
    add_group_risk("infectivity", by = "status") |>
    add_parameter(alpha = 0.01, beta = 0.2, gamma = 1, hhsize_ref = 3) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("model declares no level")
})

test_that("Status strata are declared, not read off the data", {
  # the observed labels are the model dictionary: a stratum outside it is
  # unknown, warned about, and yields no parameter
  liquid <- snake_data
  liquid$status[liquid$status == "asymptomatic"] <- "escaped"
  liquid |>
    new_housepi_model("final_size") |>
    add_group_risk(infectivity = list(
      strata = data.frame(status = c("symptomatic", "asymptomatic"),
                          param_name = c(".reference", "inf_liquid"))
    )) |>
    expect_warning("[Uu]nknown definition.*asymptomatic")

  # with no status observed at all, no sample size can pick a reference
  snake_data |>
    dplyr::mutate(status = NA) |>
    new_housepi_model("final_size") |>
    add_latent(infec_time = list(model = "uniform", min = 0, max = 5)) |>
    add_status(symptomatic = 0.3, asymptomatic = 0.7) |>
    add_group_risk(infectivity = list(by = "status", reference = ".auto")) |>
    expect_error("reference = list\\(status")
})

test_that("A status stratum reads as a covariate duplicating it", {
  # the workaround a status axis replaces: a data column carrying the status
  # of the infected members, NA elsewhere, so escapees match no rule either
  pin <- function(hhmodel, name, value) {
    hhmodel <- add_parameter(hhmodel, alpha = 0.01, beta = 0.2,
                             gamma = 1, hhsize_ref = 3)
    args <- list(hhmodel, 0.4)
    names(args) <- c("", name)
    suppressWarnings(do.call(add_parameter, args))
  }
  by_status <- snake_data |>
    new_housepi_model("final_size") |>
    add_group_risk(infectivity = list(
      by = "status", reference = list(status = "symptomatic")
    )) |>
    pin("inf_status_asymptomatic")
  by_shadow <- snake_data |>
    new_housepi_model("final_size") |>
    add_group_risk(infectivity = list(
      by = "shadow", reference = list(shadow = "symptomatic")
    )) |>
    pin("inf_shadow_asymptomatic")

  # every parameter is pinned, so the likelihood is a deterministic probe
  expect_false(anyNA(get_parameters(by_status)$value))
  expect_identical(evaluate_log_lik(by_status), evaluate_log_lik(by_shadow))
})

test_that("A status stratum drives simulated transmission", {
  # the drawn status, not a frozen covariate, selects the infectivity: with the
  # asymptomatic stratum near zero, a run that draws only asymptomatic cases
  # must transmit far less than one that draws only symptomatic ones
  outer_heaven <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 150)) |>
    dplyr::mutate(status = NA, start_date = 0, end_date = 21)

  base <- outer_heaven |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_infectivity_profile(import = "template") |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    add_status(symptomatic = 0.5, asymptomatic = 0.5) |>
    add_group_risk(infectivity = list(
      by = "status", reference = list(status = "symptomatic")
    )) |>
    add_parameter(alpha = 1e-12, beta = 20, gamma = 1, hhsize_ref = 4,
                  inf_status_asymptomatic = 1e-3) |>
    suppressWarnings()

  # the status weights must reach the categorical as literals
  mostly_sym <- base |>
    add_status(symptomatic = 0.9, asymptomatic = 0.1) |>
    suppressWarnings()
  mostly_asy <- base |>
    add_status(symptomatic = 0.1, asymptomatic = 0.9) |>
    suppressWarnings()
  attack_rate <- function(hhmodel) {
    simulated <- suppressWarnings(run_simulation(hhmodel))
    mean(as.character(simulated$status) != "escaped")
  }

  # alpha is negligible, so every case beyond the introductions is transmitted
  expect_gt(attack_rate(mostly_sym), 1.5 * attack_rate(mostly_asy))
})

##############################
## GROUP RISKS : END-TO-END ##
##############################

test_that("An absolute risk collapses the baseline it scales", {
  # alpha is never declared: the community risk carries the whole rate
  hhmodel <- flu_small |>
    new_housepi_model("final_size") |>
    add_group_risk("community", by = "member") |>
    add_parameter(beta = 0.1, gamma = 1, hhsize_ref = 4)

  expect_true(is.finite(evaluate_log_lik(hhmodel)))
  params <- get_parameters(hhmodel)
  expect_true(params$const[params$name == "alpha"])
  expect_equal(params$value[params$name == "alpha"], 1)

  # the collapsed baseline is carried into the posterior
  draws <- hhmodel |>
    run_inference(steps = 5L, thin = 1L, adapt = 0L, chains = 1L) |>
    suppressWarnings() |>
    summarise_draws()
  expect_false("alpha" %in% draws$variable)
  expect_true(any(grepl("^com_member", draws$variable)))
})

test_that("Two absolute risks on one baseline are not identifiable", {
  flu_small |>
    new_housepi_model("final_size") |>
    add_group_risk("susceptibility", by = "member") |>
    add_group_risk("infectivity", by = "member") |>
    add_parameter(alpha = 0.01, gamma = 1, hhsize_ref = 4) |>
    evaluate_log_lik() |>
    expect_error("non-identifiability.*susceptibility.*infectivity")
})

test_that("A simulation table pins the scales it provides", {
  hhmodel <- flu_small |>
    new_housepi_model("final_size") |>
    add_group_risk("susceptibility", by = "member") |>
    add_group_risk("infectivity", by = "member") |>
    add_parameter(alpha = 0.01, gamma = 1, hhsize_ref = 4)

  # every scale reads a value: no product of free scales is left to resolve
  expect_silent(run_simulation(hhmodel, param_table = pin_all_params(hhmodel)))
})

test_that("A single fixed stratum anchors the scale of its group", {
  hhmodel <- flu_small |>
    new_housepi_model("final_size") |>
    add_parameter(alpha = list(prior = "uniform", min = 0.2, max = 0.9),
                  beta = 0.1, gamma = 1, hhsize_ref = 4) |>
    add_group_risk("community", by = "member")

  # one stratum read from the table fixes the level of the whole group
  expect_silent(
    run_simulation(hhmodel, param_table = data.frame(com_member_mother = 0.5))
  )
  # the baseline keeps the prior it was given
  params <- get_parameters(hhmodel)
  expect_false(params$const[params$name == "alpha"])

  # a stratum held constant anchors the same way
  hhmodel |>
    add_parameter(com_member_student = 1.7) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_silent()
})

test_that("A pinned baseline is not collapsed onto the reference", {
  hhmodel <- flu_small |>
    new_housepi_model("final_size") |>
    add_parameter(alpha = list(prior = "uniform", min = 0, max = 1),
                  beta = 0.1, gamma = 1, hhsize_ref = 4) |>
    add_group_risk("community", by = "member")

  expect_silent(run_simulation(hhmodel, param_table = pin_all_params(hhmodel)))

  # with nothing pinning them, the baseline still carries the whole rate
  hhmodel |>
    evaluate_log_lik() |>
    expect_warning("[Bb]aseline.*constant reference")
})

test_that("A strata pin reaches simulation and stratifies on data only", {
  hhmodel <- flu_small |>
    new_housepi_model("final_size") |>
    add_parameter(alpha = 0.01, beta = 0.1, gamma = 1, hhsize_ref = 4) |>
    add_group_risk("susceptibility", strata = pin_strata(),
                   prior = "lognormal", sdlog = 3)

  # every scale reads a value, and the pin columns are no stratification axis
  expect_silent(run_simulation(hhmodel))
})

test_that("A fixed stratum holds across draws, a pinned one is resampled", {
  mk <- function(strata) {
    flu_small |>
      new_housepi_model("final_size") |>
      add_parameter(alpha = 0.01, beta = 0.1, gamma = 1, hhsize_ref = 4) |>
      add_group_risk("susceptibility", strata = strata,
                     prior = "lognormal", sdlog = 3)
  }
  twice <- function(hhmodel) {
    c(evaluate_log_lik(hhmodel), evaluate_log_lik(hhmodel))
  }

  # every stratum fixed: inference reads the values, nothing is drawn
  fixed <- mk(pin_strata(value = c(NA, 1.2, 0.5), const = c(NA, TRUE, TRUE)))
  expect_silent(out <- twice(fixed))
  expect_equal(out[1], out[2])

  # a value= pin is a simulation-only override: the prior still drives inference
  drawn <- mk(pin_strata(value = c(NA, 1.2, 0.5), const = c(NA, NA, NA)))
  expect_false(isTRUE(all.equal(twice(drawn)[1], twice(drawn)[2])))
})

test_that("A fixed stratum anchors the scale of its rule group", {
  # the group is anchored by the constant, so no baseline is collapsed onto it
  flu_small |>
    new_housepi_model("final_size") |>
    add_parameter(alpha = list(prior = "uniform", min = 0.2, max = 0.9),
                  beta = 0.1, gamma = 1, hhsize_ref = 4) |>
    add_group_risk("community", strata = pin_strata(
      value = c(NA, 1.2, 0.5), const = c(NA, TRUE, TRUE)
    ), prior = "lognormal", sdlog = 3) |>
    evaluate_log_lik() |>
    expect_silent()
})

test_that("Strata counts come from the engine cross-tabulation", {
  # a cell the data does not carry forms no stratum of its own, so the counts
  # add up to the members that carry a value, not to every member
  shadow_moses <- data.frame(
    house_id = c("h1", "h1", "h2", "h2", "h3"),
    status = c(1, 0, 1, 1, 0),
    role = c("snake", "otacon", "snake", NA, "otacon"),
    veteran = c(TRUE, FALSE, TRUE, FALSE, NA)
  )
  tbl <- shadow_moses |>
    new_housepi_model("final_size", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_group_risk(susceptibility = list(by = "role", reference = ".auto")) |>
    add_group_risk(infectivity = list(by = "veteran", reference = ".auto")) |>
    get_group_risks()

  role <- tbl[tbl$group_risk == "susceptibility", ]
  expect_setequal(as.character(role$role), c("snake", "otacon"))
  expect_equal(sum(role$n), 4L)

  # a logical column keeps both of its levels, and drops the missing one
  veteran <- tbl[tbl$group_risk == "infectivity", ]
  expect_setequal(as.character(veteran$veteran), c("TRUE", "FALSE"))
  expect_equal(sum(veteran$n), 4L)
})

test_that("A wildcard stratum survives a missing covariate cell", {
  # NA reads as missing in a data cell and as the wildcard in a rule cell: a
  # column that happens to carry an empty cell must not collapse the two
  shadow_moses <- data.frame(
    id = c("snake", "otacon", "meryl"), house_id = "h1",
    status = c(1, 0, 0), role = c("spy", "engineer", NA)
  )
  strata <- data.frame(role = c("spy", NA),
                       param_name = c(".reference", "sus_rest"))
  snap <- shadow_moses |>
    new_housepi_model("final_size", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_group_risk("susceptibility", strata = strata,
                   params = list(sus_rest = 3)) |>
    (\(x) x$export_house_snapshot("h1"))()
  risk <- stats::setNames(snap$group_risks$susceptibility, snap$group_risks$id)

  # the specific stratum sets the unit scale
  expect_equal(unname(risk["snake"]), 1)
  # the wildcard covers every other member, the ones carrying a value...
  expect_equal(unname(risk["otacon"]), 3)
  # ...and the ones that do not
  expect_equal(unname(risk["meryl"]), 3)
})

test_that("Status strata count the infected members only", {
  # the status is read off the member, not off a data column: escaped members
  # carry no infected label, so they land in no stratum
  tbl <- sarscov2_2021_icofs |>
    new_housepi_model("final_size", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_group_risk(infectivity = list(by = "status", reference = ".auto")) |>
    get_group_risks()

  expect_equal(sum(tbl$n), sum(sarscov2_2021_icofs$status > 0))
  expect_setequal(as.character(tbl$status), c("1", "2"))
})

#########################################
## GROUP RISKS : OBSERVATION COLLISION ##
#########################################

# Data carrying a numeric covariate a risk can stratify on.
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
