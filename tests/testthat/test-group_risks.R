library(testthat)

# Model over two covariates. Covariate columns stay character on purpose:
# levels follow order of first appearance, which fixes the reference stratum.
mk_risk <- local({
  hhdata <- data.frame(house_id = 1, status = 0,
                       age = c("child", "child", "adult", "adult"),
                       sex = c("M", "F", "M", "M"),
                       isolation = c(FALSE, FALSE, TRUE, TRUE))
  strata <- data.frame(isolation = c(FALSE, TRUE),
                       param_name = c(".reference", "inf_isol"))
  list(data = hhdata, strata = strata,
       model = new_housepi_model(hhdata, "final-size"))
})

# five households over a "member" covariate, "student" being the only role
# present in every household
mk_flu <- local({
  hhdata <- data.frame(
    house_id = c(1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5),
    status = c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    member = c("student", "other", "student", "sibling",
               "student", "mother", "father",
               "student", "other",
               "student", "mother", "father")
  )
  hhmodel <- new_housepi_model(hhdata, "final-size")
  list(data = hhdata, model = hhmodel,
       pinned = add_parameter(hhmodel, alpha = 0.01, beta = 0.1, gamma = 1,
                              hhsize_ref = 4))
})

# real dataset for validation: add_group_risk() copies before mutating
mk_icofs <- local({
  list(model = new_housepi_model(sarscov2_2021_icofs, "final-size"))
})

# two observed infected statuses, plus a data column duplicating them on the
# infected members only: the covariate a status stratification replaces
mk_status <- local({
  status <- c("escaped", "symptomatic", "asymptomatic", "symptomatic",
              "escaped", "asymptomatic", "symptomatic", "escaped",
              "asymptomatic", "symptomatic", "asymptomatic", "escaped")
  hhdata <- data.frame(house_id = rep(1:4, each = 3), status = status,
                       datum = ifelse(status == "escaped", NA, status))
  list(data = hhdata, model = new_housepi_model(hhdata, "final-size"))
})

# numeric covariates a risk can stratify on, and an observation can read
mk_observed <- local({
  hhmodel <- data.frame(house_id = c(1, 1), status = c(1, 0),
                        datum = c(2, 3), datum_2 = c(2, 6)) |>
    new_housepi_model("final-size")
  list(model = hhmodel)
})

# a simulation table giving every parameter of the model a value
pin_all_params <- function(hhmodel, value = 0.5) {
  names <- get_parameters(hhmodel)$name
  as.data.frame(stats::setNames(as.list(rep(value, length(names))), names))
}

# custom strata over `mk_flu`, pinning one stratum and fixing another
pin_strata <- function(value = c(NA, 1.2, 0.5), const = c(NA, NA, TRUE)) {
  data.frame(member = c("student", "mother", "father"),
             param_name = c(".reference", "sus_mum", "sus_dad"),
             value = value, const = const)
}

########################################
###  GROUP RISKS : HELPER FUNCTIONS  ###
########################################

test_that("Group risks logic table computes combinations and counts", {
  tbl <- group_risk_logic_table(mk_risk$data, c("age", "sex"))

  # counts sum to the number of rows, one count per observed combination
  expect_equal(sum(tbl$n), nrow(mk_risk$data))
  expect_equal(sort(tbl$n), c(1L, 1L, 2L))
  expect_true(all(c("age", "sex", "n") %in% colnames(tbl)))
})

test_that("Associate param names decorates table and marks reference", {
  tbl <- group_risk_logic_table(mk_risk$data, c("age", "sex"))
  res <- associate_group_param_names(tbl, prefix = "sus", reference = ".auto")

  # exactly one .reference row (largest group)
  expect_equal(sum(res$param_name == ".reference"), 1L)
  # other param_names start with the provided prefix
  others <- res$param_name[res$param_name != ".reference"]
  expect_true(all(grepl("^sus_", others)))

  # absolute by default: every stratum keeps its own parameter
  abs_res <- associate_group_param_names(tbl, prefix = "sus", reference = NULL)
  expect_false(".reference" %in% abs_res$param_name)
  expect_true(all(grepl("^sus_", abs_res$param_name)))
})

test_that("An empty logic table passes through grouping and ordering", {
  tbl <- compute_rule_groups(
    data.frame(age = character(), param_name = character(), n = integer())
  )
  expect_equal(nrow(tbl), 0L)
  expect_true("rule_group" %in% colnames(tbl))
  expect_identical(reorder_risks_logic_table(tbl), tbl)
})

test_that("Custom strata logic table is validated", {
  # must be a data.frame
  mk_risk$model |>
    add_group_risk("susceptibility", strata = list(a = 1)) |>
    expect_error("data\\.frame")
  # must have param_name column
  mk_risk$model |>
    add_group_risk("susceptibility", strata = data.frame(x = 1)) |>
    expect_error("param_name")
  # must have at least one row
  mk_risk$model |>
    add_group_risk("susceptibility",
                   strata = data.frame(param_name = character(0))) |>
    expect_error("at least one row")
  # reserved covariate column names rejected
  mk_risk$model |>
    add_group_risk("susceptibility",
                   strata = data.frame(param_name = "a", id = 1)) |>
    expect_error("[Rr]eserved")
})

######################################
###  GROUP RISKS : STANDARD CALLS  ###
######################################

test_that("Group risks integrate with data", {
  hhmodel <- mk_risk$model |>
    add_group_risk("community", by = c("age", "sex"), joint = TRUE)

  tidy <- get_group_risks(hhmodel)
  expect_true(all(tidy$group_risk == "community"))
  # joint expansion: one rule group over both covariates
  expect_equal(length(unique(tidy$rule_group)), 1L)
  expect_true(all(c("age", "sex") %in% colnames(tidy)))
  # absolute by default: no reference stratum
  expect_equal(sum(tidy$param_name == ".reference"), 0L)

  # mcmc setup is the absolute default
  params <- get_parameters(hhmodel)
  params <- params[!params$const, ]
  expect_true(all(params$min == 0))
  expect_true(all(params$max == 1))
  expect_true(all(params$proposal == "rw_log"))
})

test_that("Reference turns group risks relative with its own mcmc setup", {
  hhmodel <- mk_risk$model |>
    add_group_risk("community", by = c("age", "sex"), joint = TRUE,
                   reference = ".auto")

  tidy <- get_group_risks(hhmodel)
  expect_equal(sum(tidy$param_name == ".reference"), 1L)
  # the reference is the largest stratum
  expect_equal(tidy$n[tidy$param_name == ".reference"], max(tidy$n))

  # mcmc setup is default (Cauchemez et al., 2009)
  params <- get_parameters(hhmodel)
  expect_equal(params$args[[1]], c(meanlog = 0, sdlog = 3))
  expect_true(all(params$min == 0))
  expect_true(all(params$max == 10))
  expect_true(all(params$proposal == "rw_log"))
})

test_that("Group risks hold without data and re-derive on add_data_wide", {
  # risks ingest data: a data-less model holds them, no risk params yet
  hhmodel <- new_housepi_model(model = "final-size") |>
    add_parameter(import = "template") |>
    # relative: the template's baseline stays free alongside the risk
    add_group_risk("susceptibility", by = "member", reference = ".auto")
  expect_false(any(grepl("^sus_member", hhmodel$parameters$name)))

  # attaching data materialises the risk parameters and runs
  hhmodel_2 <- add_data_wide(hhmodel, mk_flu$data)
  full_params <- hhmodel_2$parameters$name
  expect_true(any(grepl("^sus_member", full_params)))
  expect_true(is.finite(evaluate_log_lik(hhmodel_2)))

  # replacing data re-derives: a single-level subset drops the stale risk params
  hhmodel_3 <- add_data_wide(hhmodel_2,
                             mk_flu$data[mk_flu$data$member == "student", ])
  expect_false(any(grepl("^sus_member", hhmodel_3$parameters$name)))
  expect_true(is.finite(evaluate_log_lik(hhmodel_3)))

  # re-attaching the full data restores the exact original parameter set
  add_data_wide(hhmodel_3, mk_flu$data)$parameters$name |>
    expect_setequal(full_params)
})

test_that("Group risks expand on their own covariates only", {
  hhdata <- sarscov2_2021_icofs |>
    dplyr::select(-start_date) |>
    dplyr::filter(house_id %in% utils::head(unique(house_id), 50))

  # a data column the risk does not read leaves the expansion untouched
  risks <- hhdata |>
    new_housepi_model("case-ascertained") |>
    add_group_risk("susceptibility", by = "age") |>
    get_group_risks()
  expect_setequal(risks$group_risk, "susceptibility")
  expect_true("sus_age_child" %in% risks$param_name)

  # without data the read surface warns and keeps the table shape
  expect_warning(
    risks <- new_housepi_model(model = "final-size") |>
      add_group_risk("susceptibility", by = "member") |>
      get_group_risks(),
    regexp = "without data"
  )
  expect_equal(nrow(risks), 0L)
  expect_named(risks, c("group_risk", "param_name", "n", "rule_group"))
})

test_that("Group risks sentinels work correctly", {
  # empty risk definition
  mk_icofs$model |>
    add_group_risk("susceptibility", by = NA) |>
    expect_error("[Ii]nvalid risk definition")

  # wrong type for covariate column choice
  mk_icofs$model |>
    add_group_risk("susceptibility", by = 4) |>
    expect_error("[Aa]rgument .by. must be of type character")

  # not implemented group
  mk_icofs$model |>
    add_group_risk("unknown", by = "sex") |>
    expect_error("[Uu]nused group risk.*unknown")

  # incomplete custom prior definition
  mk_icofs$model |>
    add_group_risk("susceptibility", by = "sex", max = 3) |>
    expect_error("[Ss]pare prior arguments")

  mk_icofs$model |>
    add_group_risk("susceptibility", by = "test_date") |>
    expect_error("[Mm]ust be a unidimensional type")
})

test_that("An automatic reference names an observed stratum of its axes", {
  mk_risk$model |>
    add_group_risk("susceptibility", by = "age",
                   reference = list(age = "elderly")) |>
    expect_error("individual with age = .elderly. is never observed")
  mk_risk$model |>
    add_group_risk("susceptibility", by = c("age", "sex"),
                   reference = list(age = "adult")) |>
    expect_error("by. and .reference. should have equal dimensions")
  mk_risk$model |>
    add_group_risk("susceptibility", by = "age",
                   reference = list(sex = "M")) |>
    expect_error("[Uu]nrecognized names in reference list")
})

test_that("Group risks reference support behaves as expected", {
  # an untagged custom table is absolute
  strata <- data.frame(sex = c("F", "M"),
                       param_name = c("ladies", "gentlemen"))
  mk_icofs$model |>
    add_group_risk("susceptibility", strata = strata) |>
    expect_no_error()
  # reference has no meaning against a custom table
  mk_icofs$model |>
    add_group_risk("susceptibility", reference = ".auto", strata = strata) |>
    expect_error("[Nn]ot used with a custom")

  # a tagged row makes the same table relative
  strata$param_name[1] <- ".reference"
  mk_icofs$model |>
    add_group_risk("susceptibility", strata = strata) |>
    get_group_risks() |>
    subset(param_name == ".reference") |>
    nrow() |>
    expect_equal(1L)
})

test_that("Custom strata carry at most one reference per rule group", {
  # two references within the same rule group are two competing scales
  mk_icofs$model |>
    add_group_risk("susceptibility", strata = data.frame(
      sex = c("F", "M"), param_name = c(".reference", ".reference")
    )) |>
    expect_error("[Rr]ows.*one rule group")

  # distinct rule groups each carry their own reference
  mk_icofs$model |>
    add_group_risk("susceptibility", strata = data.frame(
      sex = c("F", "M", NA, NA),
      age = c(NA, NA, "child", "adult"),
      param_name = c(".reference", "gentlemen", ".reference", "grown_ups")
    )) |>
    expect_no_error()
})

test_that("Parameter options work correctly", {
  # the collapsed baseline is not one of the risk parameters
  sampled_params <- function(hhmodel) {
    params <- get_parameters(hhmodel)
    params[!params$const, ]
  }

  params <- mk_risk$model |>
    add_group_risk("susceptibility", by = c("age", "sex"),
                   prior = "uniform", min = 0, max = 10, joint = TRUE) |>
    sampled_params()
  expect_true(all(params$prior == "uniform"))
  expect_true(all(params$min == 0))
  expect_true(all(params$max == 10))
  # a sampled risk parameter defaults to a random walk
  expect_true(all(params$proposal == "rw"))

  params <- mk_risk$model |>
    add_group_risk("infectivity", by = c("age", "sex"),
                   prior = "gamma", shape = 2, scale = 1, joint = TRUE) |>
    sampled_params()
  expect_true(all(params$prior == "gamma"))
  expect_equal(params$args[[1]], c(shape = 2, scale = 1))
  expect_true(all(params$min == 0))
  expect_true(all(is.na(params$max)))
  expect_true(all(params$proposal == "rw"))
})

test_that("Load list works as expected", {
  tidy <- mk_risk$model |>
    add_group_risk(
      community = list(by = c("age", "sex"), prefix = "ct",
                       prior = "uniform", min = 0, max = 10),
      susceptibility = list(by = c("age", "sex"),
                            prior = "lognormal", sdlog = 3)
    ) |>
    get_group_risks()

  expect_setequal(unique(tidy$group_risk), c("community", "susceptibility"))
  expect_true(any(grepl("^ct_", tidy$param_name)))
  expect_true(any(grepl("^sus_", tidy$param_name)))
})

#################################################
###  GROUP RISKS : CUSTOM STRATA DEFINITIONS  ###
#################################################

test_that("Group risks deal with custom logic table correctly", {
  hhmodel <- mk_risk$model |>
    add_group_risk("infectivity", strata = mk_risk$strata,
                   prior = "lognormal", meanlog = 0, sdlog = 2)

  tidy <- get_group_risks(hhmodel)
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

test_that("A stratum naming an unobserved level is skipped", {
  # the data carries no elderly: the stratum names a level that would encode
  # to NA, the wildcard, so it is dropped rather than widened
  strata <- data.frame(age = c("child", "elderly"),
                       param_name = c(".reference", "sus_elderly"))
  expect_warning(
    hhmodel <- mk_risk$model |>
      add_group_risk("susceptibility", strata = strata,
                     prior = "lognormal", meanlog = 0, sdlog = 2),
    "[Uu]nknown definition.*elderly"
  )

  # only the reference survives, and it asks for no parameter
  expect_setequal(get_group_risks(hhmodel)$param_name, ".reference")
  expect_false("sus_elderly" %in% get_parameters(hhmodel)$name)
})

test_that("A strata table emptied by unknown levels stratifies nothing", {
  # neither level is in the data, so every stratum goes: the risk is left
  # with the reference alone, which every member reads
  strata <- data.frame(age = c("elderly", "infant"),
                       param_name = c(".reference", "sus_infant"))
  expect_warning(
    hhmodel <- mk_risk$model |>
      add_group_risk("susceptibility", strata = strata,
                     prior = "lognormal", meanlog = 0, sdlog = 2),
    "[Uu]nknown definitions.*elderly"
  )

  expect_equal(get_group_risks(hhmodel)$param_name, ".reference")
  expect_false("sus_infant" %in% get_parameters(hhmodel)$name)
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
  hhmodel <- mk_risk$data |>
    transform(age = c("child", "adult", "elderly", "adult")) |>
    new_housepi_model("final-size") |>
    add_group_risk("susceptibility", strata = strata,
                   params = overrides, prior = "gamma", shape = 2, scale = 1)

  params <- get_parameters(hhmodel)
  expect_setequal(params$name, c("sus_old", "sus_elderly"))
  # overrides used instead of raw params
  old <- get_parameters(hhmodel, "sus_old")
  expect_equal(old$value, 1.73)
  expect_true(old$const)
  elderly <- get_parameters(hhmodel, "sus_elderly")
  expect_equal(elderly$prior, "lognormal")
  expect_equal(elderly$args[[1]], c(sdlog = 1.2, meanlog = 0))

  hhmodel |>
    add_group_risk("susceptibility",
                   strata = data.frame(
                     sex = c("M", "F"), param_name = c("sus_male", ".reference")
                   ),
                   params = list(unknown = 2.4)) |>
    expect_error("[Uu]nknown")
})

test_that("Group risk mixed auto and custom modes preserve strata entries", {
  hhmodel <- mk_risk$model |>
    add_group_risk(
      community = list(by = "age", prefix = "com"),
      infectivity = list(strata = mk_risk$strata)
    )

  tidy <- get_group_risks(hhmodel)
  expect_setequal(unique(tidy$group_risk), c("community", "infectivity"))
  # strata entry is unchanged
  subset(tidy, group_risk == "infectivity")$param_name |>
    expect_setequal(c(".reference", "inf_isol"))
  # both modes materialise their parameters
  param_names <- get_parameters(hhmodel)$name
  expect_true("inf_isol" %in% param_names)
  expect_true(any(grepl("^com_", param_names)))

  add_group_risk(hhmodel, "infectivity", by = "age") |>
    expect_warning("[Oo]verwriting.*infectivity")
  add_group_risk(hhmodel, "community", strata = mk_risk$strata) |>
    expect_warning("[Oo]verwriting.*community")

  # a covariate the data does not carry fails as soon as data is attached
  data.frame(house_id = 1, status = 0, age = c("child", "adult")) |>
    new_housepi_model("final-size") |>
    add_group_risk("infectivity", strata = mk_risk$strata) |>
    expect_error("[Uu]nknown covariate.*isolation")
})

test_that("Absoluteness follows the reference across every definition mode", {
  # the default prior tracks whether a risk carries absolute values
  risk_prior <- function(...) {
    params <- get_parameters(add_group_risk(mk_risk$model, ...))
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

###################################
###  GROUP RISKS : STRATA PINS  ###
###################################

test_that("A strata value pins the stratum and keeps its prior", {
  hhmodel <- mk_flu$model |>
    add_group_risk("susceptibility", strata = pin_strata(),
                   prior = "lognormal", sdlog = 3)
  mum <- get_parameters(hhmodel, "sus_mum")
  expect_equal(mum$value, 1.2)
  expect_false(mum$const)
  expect_equal(mum$prior, "lognormal")

  # a fixed stratum is a constant parameter: it carries no prior
  dad <- get_parameters(hhmodel, "sus_dad")
  expect_equal(dad$value, 0.5)
  expect_true(dad$const)
  expect_true(is.na(dad$prior))

  # pin columns surface on the read table only when a pin exists
  tidy <- get_group_risks(hhmodel)
  expect_true(all(c("value", "const") %in% colnames(tidy)))
  # padding keeps the columns typed across sub-tables
  expect_type(tidy$value, "double")
  expect_type(tidy$const, "logical")
  expect_equal(tidy$value[tidy$param_name == "sus_mum"], 1.2)
  expect_true(is.na(tidy$value[tidy$param_name == ".reference"]))

  # an unpinned risk keeps the table it had
  mk_flu$model |>
    add_group_risk("susceptibility", by = "member") |>
    get_group_risks() |>
    colnames() |>
    intersect(c("value", "const")) |>
    expect_length(0)
})

test_that("Invalid strata pins are rejected", {
  mk <- function(strata) {
    add_group_risk(mk_flu$model, "susceptibility", strata = strata,
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

  # a stratum spread over two rules carries a single pin, even when the data
  # carries no grandparent to keep the second rule
  data.frame(member = c("student", "mother", "grandparent"),
             param_name = c(".reference", "sus_mum", "sus_mum"),
             value = c(NA, 1.2, 0.4)) |>
    mk() |>
    expect_error("[Cc]onflicting.*value.*sus_mum.*rows 2 and 3")
})

test_that("A pinned stratum takes no entry in params", {
  mk_flu$model |>
    add_group_risk("susceptibility", strata = pin_strata(),
                   params = list(sus_mum = 2), prior = "lognormal",
                   sdlog = 3) |>
    expect_error("sus_mum")
})

#####################################
###  GROUP RISKS : STATUS STRATA  ###
#####################################

test_that("Status is accepted as a stratification axis", {
  # both gates let 'status' through for a group risk
  mk_status$model |>
    add_group_risk("infectivity", by = "status") |>
    get_group_risks() |>
    getElement("status") |>
    as.character() |>
    expect_setequal(c("symptomatic", "asymptomatic"))
  mk_status$model |>
    add_group_risk("infectivity",
                   strata = data.frame(status = "symptomatic",
                                       param_name = "inf_symptomatic")) |>
    get_parameters() |>
    subset(!const) |>
    getElement("name") |>
    expect_identical("inf_symptomatic")

  # one parameter per non-reference stratum, named after the label
  mk_status$model |>
    add_group_risk(infectivity = list(
      by = "status", reference = list(status = "symptomatic")
    )) |>
    get_parameters() |>
    getElement("name") |>
    expect_identical("inf_status_asymptomatic")
})

test_that("Only the risks the model declares may be status stratified", {
  # infectivity is read after the member's status is settled, the others are
  # read on members that carry no infected status
  mk_status$model |>
    add_group_risk("susceptibility", by = "status") |>
    expect_error("cannot be stratified by the infected status")
  mk_status$model |>
    add_group_risk("community", by = "status") |>
    expect_error("cannot be stratified by the infected status")

  # contact risks keep no status axis at all
  mk_status$model |>
    add_contact_risk(by = "status") |>
    expect_error("reserved")

  # a model that declares no status level cannot stratify on it. The risk
  # waits for the declaration, so decorator order is free, and the wall reports
  # that it never came
  expect_no_error(
    hhmodel <- add_group_risk(mk_risk$model, "infectivity", by = "status")
  )
  hhmodel |>
    add_parameter(alpha = 0.01, beta = 0.2, gamma = 1, hhsize_ref = 3) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("model declares no level")
})

test_that("Status strata are declared, not read off the data", {
  # the observed labels are the model dictionary: a stratum outside it is
  # unknown, warned about, and yields no parameter
  hhdata <- mk_status$data
  hhdata$status[hhdata$status == "asymptomatic"] <- "escaped"
  hhdata |>
    new_housepi_model("final-size") |>
    add_group_risk(infectivity = list(
      strata = data.frame(status = c("symptomatic", "asymptomatic"),
                          param_name = c(".reference", "inf_asymptomatic"))
    )) |>
    expect_warning("[Uu]nknown definition.*asymptomatic")

  # with no status observed at all, no sample size can pick a reference
  hhmodel <- mk_status$data |>
    dplyr::mutate(status = NA) |>
    new_housepi_model("final-size") |>
    add_latent(infec_time = list(model = "uniform", min = 0, max = 5)) |>
    add_status(symptomatic = 0.3, asymptomatic = 0.7)
  hhmodel |>
    add_group_risk(infectivity = list(by = "status", reference = ".auto")) |>
    expect_error("reference = list\\(status")

  # a custom table stands all the same, its sample sizes unknown
  hhmodel |>
    add_group_risk(infectivity = list(
      strata = data.frame(status = c("symptomatic", "asymptomatic"),
                          param_name = c(".reference", "inf_asymptomatic"))
    )) |>
    get_group_risks() |>
    getElement("n") |>
    expect_equal(c(NA_integer_, NA_integer_))
})

test_that("A status stratum reads as a covariate duplicating it", {
  # the workaround a status axis replaces: a data column carrying the status
  # of the infected members, NA elsewhere, so escapees match no rule either
  pin <- function(hhmodel, name) {
    hhmodel <- add_parameter(hhmodel, alpha = 0.01, beta = 0.2,
                             gamma = 1, hhsize_ref = 3)
    suppressWarnings(add_parameter(hhmodel, import = stats::setNames(
      list(0.4), name
    )))
  }
  by_status <- mk_status$model |>
    add_group_risk(infectivity = list(
      by = "status", reference = list(status = "symptomatic")
    )) |>
    pin("inf_status_asymptomatic")
  by_datum <- mk_status$model |>
    add_group_risk(infectivity = list(
      by = "datum", reference = list(datum = "symptomatic")
    )) |>
    pin("inf_datum_asymptomatic")

  # every parameter is pinned, so the likelihood is a deterministic probe
  expect_false(anyNA(get_parameters(by_status)$value))
  expect_identical(evaluate_log_lik(by_status), evaluate_log_lik(by_datum))
})

test_that("A status stratum drives simulated transmission", {
  # the drawn status, not a frozen covariate, selects the infectivity: with the
  # asymptomatic stratum near zero, a run that draws only asymptomatic cases
  # must transmit far less than one that draws only symptomatic ones
  hhmodel <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 150)) |>
    dplyr::mutate(status = NA, start_date = 0, end_date = 21) |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
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
  attack_rate <- function(symptomatic) {
    simulated <- hhmodel |>
      add_status(symptomatic = !!symptomatic,
                 asymptomatic = !!(1 - symptomatic)) |>
      run_simulation() |>
      suppressWarnings()
    mean(as.character(simulated$status) != "escaped")
  }

  # alpha is negligible, so every case beyond the introductions is transmitted
  expect_gt(attack_rate(0.9), 1.5 * attack_rate(0.1))
})

test_that("Status strata count the infected members only", {
  # the status is read off the member, not off a data column: escaped members
  # carry no infected label, so they land in no stratum
  tbl <- mk_icofs$model |>
    add_parameter(import = "template") |>
    add_group_risk(infectivity = list(by = "status", reference = ".auto")) |>
    get_group_risks()

  expect_equal(sum(tbl$n), sum(sarscov2_2021_icofs$status > 0))
  expect_setequal(as.character(tbl$status), c("1", "2"))
})

##############################
## GROUP RISKS : END-TO-END ##
##############################

test_that("An absolute risk collapses the baseline it scales", {
  # alpha is never declared: the community risk carries the whole rate
  hhmodel <- mk_flu$model |>
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

test_that("A baseline pinned by an absolute risk follows its redeclaration", {
  hhmodel <- mk_flu$model |>
    add_group_risk("community", by = "member", reference = "student") |>
    add_parameter(beta = 0.1, gamma = 1, hhsize_ref = 4)
  # the same risk declared absolute first, then relative
  hhmodel_2 <- mk_flu$model |>
    add_group_risk("community", by = "member") |>
    add_parameter(beta = 0.1, gamma = 1, hhsize_ref = 4)
  expect_warning(
    hhmodel_2 <- hhmodel_2 |>
      add_group_risk("community", by = "member", reference = "student"),
    "[Oo]verwriting.*community"
  )

  by_name <- function(x) x[order(x$name), , drop = FALSE]
  by_name(get_parameters(hhmodel_2)) |>
    expect_equal(by_name(get_parameters(hhmodel)), ignore_attr = TRUE)
  expect_false("alpha" %in% get_parameters(hhmodel_2)$name)
  hhmodel_2 |>
    evaluate_log_lik() |>
    expect_error("[Mm]issing parameter definition.*alpha")
})

test_that("Two absolute risks on one baseline are identifiable once pinned", {
  hhmodel <- mk_flu$model |>
    add_group_risk("susceptibility", by = "member") |>
    add_group_risk("infectivity", by = "member") |>
    add_parameter(alpha = 0.01, gamma = 1, hhsize_ref = 4)
  hhmodel |>
    evaluate_log_lik() |>
    expect_error("[Nn]on-identifiability.*susceptibility.*infectivity")

  # a simulation table gives every scale a value: no product of free scales
  # is left to resolve
  run_simulation(hhmodel, param_table = pin_all_params(hhmodel)) |>
    expect_s3_class("data.frame")
})

test_that("A single fixed stratum anchors the scale of its group", {
  hhmodel <- mk_flu$model |>
    add_parameter(alpha = list(prior = "uniform", min = 0.2, max = 0.9),
                  beta = 0.1, gamma = 1, hhsize_ref = 4) |>
    add_group_risk("community", by = "member")

  # one stratum read from the table fixes the level of the whole group
  hhmodel |>
    run_simulation(param_table = data.frame(com_member_mother = 0.5)) |>
    expect_s3_class("data.frame")
  # the baseline keeps the prior it was given
  params <- get_parameters(hhmodel)
  expect_false(params$const[params$name == "alpha"])

  # a stratum held constant anchors the same way
  hhmodel |>
    add_parameter(com_member_student = 1.7) |>
    suppressWarnings() |> # suppress overwrite warning
    evaluate_log_lik() |>
    is.finite() |>
    expect_true()
})

test_that("A pinned baseline is not collapsed onto the reference", {
  hhmodel <- mk_flu$model |>
    add_parameter(alpha = list(prior = "uniform", min = 0, max = 1),
                  beta = 0.1, gamma = 1, hhsize_ref = 4) |>
    add_group_risk("community", by = "member")

  run_simulation(hhmodel, param_table = pin_all_params(hhmodel)) |>
    expect_s3_class("data.frame")

  # with nothing pinning them, the baseline still carries the whole rate
  hhmodel |>
    evaluate_log_lik() |>
    expect_warning("[Bb]aseline.*constant reference")
})

test_that("A strata pin reaches simulation and stratifies on data only", {
  # every scale reads a value, and the pin columns are no stratification axis
  mk_flu$pinned |>
    add_group_risk("susceptibility", strata = pin_strata(),
                   prior = "lognormal", sdlog = 3) |>
    run_simulation() |>
    expect_s3_class("data.frame")
})

test_that("A fixed stratum holds across draws, a pinned one is resampled", {
  mk <- function(strata) {
    mk_flu$pinned |>
      add_group_risk("susceptibility", strata = strata,
                     prior = "lognormal", sdlog = 3)
  }
  twice <- function(hhmodel) {
    c(evaluate_log_lik(hhmodel), evaluate_log_lik(hhmodel))
  }

  # every stratum fixed: inference reads the values, nothing is drawn
  log_lik <- twice(mk(pin_strata(value = c(NA, 1.2, 0.5),
                                 const = c(NA, TRUE, TRUE))))
  log_lik[1] |> expect_equal(log_lik[2])

  # a value= pin is a simulation-only override: the prior still drives inference
  log_lik <- twice(mk(pin_strata(value = c(NA, 1.2, 0.5),
                                 const = c(NA, NA, NA))))
  expect_false(isTRUE(all.equal(log_lik[1], log_lik[2])))
})

test_that("A fixed stratum anchors the scale of its rule group", {
  # the group is anchored by the constant, so no baseline is collapsed onto it
  log_lik <- mk_flu$model |>
    add_parameter(alpha = list(prior = "uniform", min = 0.2, max = 0.9),
                  beta = 0.1, gamma = 1, hhsize_ref = 4) |>
    add_group_risk("community",
      strata = pin_strata(value = c(NA, 1.2, 0.5), const = c(NA, TRUE, TRUE)),
      prior = "lognormal", sdlog = 3
    ) |>
    evaluate_log_lik()
  expect_true(is.finite(log_lik))
})

test_that("Strata counts come from the engine cross-tabulation", {
  # a cell the data does not carry forms no stratum of its own, so the counts
  # add up to the members that carry a value, not to every member
  hhmodel <- data.frame(house_id = c("h1", "h1", "h2", "h2", "h3"),
                        status = c(1, 0, 1, 1, 0),
                        datum = c("a", "b", "a", NA, "b"),
                        datum_2 = c(TRUE, FALSE, TRUE, FALSE, NA)) |>
    new_housepi_model("final-size", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_group_risk(susceptibility = list(by = "datum", reference = ".auto"),
                   infectivity = list(by = "datum_2", reference = ".auto"))

  datum <- hhmodel |> get_group_risks("susceptibility")
  expect_setequal(as.character(datum$datum), c("a", "b"))
  sum(datum$n) |> expect_equal(4L)

  # a logical column keeps both of its levels, and drops the missing one
  datum_2 <- hhmodel |> get_group_risks("infectivity")
  expect_setequal(as.character(datum_2$datum_2), c("TRUE", "FALSE"))
  sum(datum_2$n) |> expect_equal(4L)
})

test_that("A wildcard stratum survives a missing covariate cell", {
  # NA reads as missing in a data cell and as the wildcard in a rule cell: a
  # column that happens to carry an empty cell must not collapse the two
  strata <- data.frame(datum = c("a", NA),
                       param_name = c(".reference", "sus_rest"))
  snapshot <- data.frame(id = c("m1", "m2", "m3"), house_id = "h1",
                         status = c(1, 0, 0), datum = c("a", "b", NA)) |>
    new_housepi_model("final-size", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_group_risk("susceptibility", strata = strata,
                   params = list(sus_rest = 3)) |>
    inspect_house_state("h1")
  risk <- stats::setNames(snapshot$group_risks$susceptibility,
                          snapshot$group_risks$id)

  # the specific stratum sets the unit scale
  expect_equal(unname(risk["m1"]), 1)
  # the wildcard covers every other member, the ones carrying a value...
  expect_equal(unname(risk["m2"]), 3)
  # ...and the ones that do not
  expect_equal(unname(risk["m3"]), 3)
})

test_that("A wildcard stratum counts the members it resolves to", {
  # the specific stratum takes `m1` off the wildcard
  strata <- data.frame(datum = c("a", NA),
                       param_name = c(".reference", "sus_rest"))
  hhmodel <- data.frame(id = c("m1", "m2", "m3"), house_id = "h1",
                        status = c(1, 0, 0), datum = c("a", "b", NA)) |>
    new_housepi_model("final-size", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_group_risk("susceptibility", strata = strata,
                   params = list(sus_rest = 3))
  risks <- get_group_risks(hhmodel)

  risk <- inspect_house_state(hhmodel, "h1")$group_risks$susceptibility
  risks$n[risks$param_name == ".reference"] |> expect_equal(sum(risk == 1))
  risks$n[risks$param_name == "sus_rest"] |> expect_equal(sum(risk == 3))
})

#########################################
## GROUP RISKS : OBSERVATION COLLISION ##
#########################################

test_that("A group risk keyed on an observation is rejected", {
  # the risk feeds the transmission model the observation is generated from
  mk_observed$model |>
    add_group_risk("susceptibility", by = "datum") |>
    add_observation("datum", model = "uniform", min = 0, max = 5) |>
    expect_error("keyed on observation.*datum")

  # the reverse declaration order fails just the same
  mk_observed$model |>
    add_observation("datum", model = "uniform", min = 0, max = 5) |>
    add_group_risk("susceptibility", by = "datum") |>
    expect_error("keyed on observation.*datum")

  # a custom strata table takes no loophole
  mk_observed$model |>
    add_observation("datum", model = "uniform", min = 0, max = 5) |>
    add_group_risk("susceptibility", strata = data.frame(
      datum = c(2, 3), param_name = c(".reference", "sus_high")
    )) |>
    expect_error("keyed on observation.*datum")

  # a deterministic household observation standing in for its data column is
  # derived, an observation all the same
  mk_observed$model |>
    add_observation(datum = min(datum_2@D), level = "house") |>
    add_group_risk("susceptibility", by = "datum") |>
    expect_error("keyed on observation.*datum")
})
