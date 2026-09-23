library(testthat)

# model with no data of its own beyond the 'test' default
mk_wrapper <- local({
  list(model = new_housepi_model(model = "test"))
})

# households over a member covariate the risk wrappers stratify on
mk_member <- local({
  hhmodel <- data.frame(
    hhid = c(1, 1, 2, 2, 2),
    status = c(1, 0, 1, 0, 0),
    member = c("student", "other", "student", "mother", "father")
  ) |>
    new_housepi_model("final-size")
  list(model = hhmodel)
})

###################################
## EXPOSE CLASS : DISPATCH       ##
###################################

test_that("Wrappers reject objects that are not a housepi model", {
  obj <- "not a model"
  err <- "must be of class"
  add_data_long(obj, data.frame()) |> expect_error(err)
  add_data_wide(obj, data.frame()) |> expect_error(err)
  get_data_wide(obj) |> expect_error(err)
  get_data_long(obj) |> expect_error(err)
  get_parameters(obj) |> expect_error(err)
  get_latent(obj) |> expect_error(err)
  get_group_risks(obj) |> expect_error(err)
  get_infectivity_profiles(obj) |> expect_error(err)
  add_contact_risk(obj) |> expect_error(err)
  add_parameter(obj) |> expect_error(err)
  add_latent(obj) |> expect_error(err)
  add_infectivity_profile(obj, "gamma") |> expect_error(err)
  evaluate_log_lik(obj) |> expect_error(err)
  run_inference(obj) |> expect_error(err)
  run_simulation(obj) |> expect_error(err)
  add_group_risk(obj) |> expect_error(err)
  run_custom_predictive_check(obj, NULL) |> expect_error(err)
  run_epi_predictive_check(obj, NULL) |> expect_error(err)
})

test_that("Data provided at construction is exposed via get_data_wide", {
  hhdata <- get_data_wide(new_housepi_model(flu_1918_baltimore, "test"))
  expect_true(all(c("hhid", "status") %in% colnames(hhdata)))

  # id automatically added
  expect_identical(hhdata$id, as.character(seq_len(nrow(hhdata))))
  # hhid cast to character
  expect_identical(hhdata$hhid, as.character(flu_1918_baltimore$hhid))
  # statuses are integer
  expect_identical(hhdata$status, flu_1918_baltimore$status)
})

test_that("Variadic labels are not swallowed by formal matching", {
  # single-letter categorical weights c/d/p would partial-match the
  # const/discr/prior formals of add_parameter; they reach the prior call
  param <- mk_wrapper$model |>
    add_parameter("rho", prior = "categorical", a = 0.2, b = 0.1, c = 0.7) |>
    get_parameters()
  expect_equal(param$prior, "categorical")
  expect_equal(param$args[[1]], c(a = 0.2, b = 0.1, c = 0.7))

  # bulk form
  param <- mk_wrapper$model |>
    add_parameter(
      rho = list(prior = "categorical", a = 0.2, b = 0.1, c = 0.7),
      rho_2 = list(prior = "categorical", d = 0.4, p = 0.6)
    ) |>
    get_parameters()
  expect_equal(param$name, c("rho", "rho_2"))
  expect_equal(param$prior, c("categorical", "categorical"))
  expect_equal(param$args[[1]], c(a = 0.2, b = 0.1, c = 0.7))
  expect_equal(param$args[[2]], c(d = 0.4, p = 0.6))
})

test_that("Bulk import mirrors inline name = value pairs", {
  defs <- list(rho = 1, rho_2 = list(prior = "uniform"))

  # an external list passed via `import` equals the inline dotted form
  mk_wrapper$model |>
    add_parameter(import = defs) |>
    get_parameters() |>
    expect_identical(
      get_parameters(add_parameter(mk_wrapper$model, rho = 1,
                                   rho_2 = list(prior = "uniform")))
    )

  # `import` cannot be combined with a single name or with inline pairs
  mk_wrapper$model |>
    add_parameter("rho_3", import = defs) |>
    expect_error("[Mm]ixing")
  mk_wrapper$model |>
    add_parameter(rho_3 = 1, import = defs) |>
    expect_error("[Mm]ixing")

  # the template sentinel follows the same exclusion rules
  mk_wrapper$model |>
    add_parameter("rho_3", import = "template") |>
    expect_error("[Mm]ixing")
  mk_wrapper$model |>
    add_latent(hidden = 1, import = "template") |>
    expect_error("[Mm]ixing")

  # a character import naming no known source is a typo, not a definition
  mk_wrapper$model |>
    add_parameter(import = "templte") |>
    expect_error("[Uu]nknown import source")
})

test_that("Bulk import resolves symbols in the calling frame", {
  # `datum_local` exists only inside the caller
  param <- local({
    datum_local <- 3
    mk_wrapper$model |>
      add_parameter(import = list(rho = quote(datum_local * 2))) |>
      get_parameters()
  })
  expect_equal(param$value, 6)
})

test_that("add_observation routes every input branch on `level`", {
  hhmodel <- data.frame(hhid = "h1", status = 1, record = 3, datum = 2) |>
    new_housepi_model("test")

  # default level is member
  hhmodel |>
    add_observation("record", model = "normal") |>
    get_observations() |>
    getElement("level") |>
    expect_identical("member")

  # single, bulk and imported definitions
  single <- hhmodel |>
    add_observation("record_2", value = min(datum@D), level = "house")
  bulk <- hhmodel |>
    add_observation(record_2 = min(datum@D), level = "house")
  imported <- hhmodel |>
    add_observation(import = list(record_2 = quote(min(datum@D))),
                    level = "house")

  expect_identical(get_observations(bulk), get_observations(single))
  expect_identical(get_observations(imported), get_observations(single))
  expect_identical(get_observations(single)[["level"]], "house")

  # a member-level function reaches the definition through `func`
  single <- hhmodel |>
    add_observation("record_2", func = min(datum@D), status = 1)
  bulk <- hhmodel |>
    add_observation(record_2 = list(func = min(datum@D), status = 1))
  imported <- hhmodel |>
    add_observation(
      import = list(record_2 = list(func = quote(min(datum@D)), status = 1))
    )

  expect_identical(get_observations(bulk), get_observations(single))
  expect_identical(get_observations(imported), get_observations(single))
  expect_identical(get_observations(single)[["func"]], "min")

  # an unknown level is rejected before anything is stored
  hhmodel |>
    add_observation(record_2 = min(datum@D), level = "unknown") |>
    expect_error("should be one of")
})

test_that("Advanced setup options route through `...`", {
  hhmodel <- mk_wrapper$model |>
    add_parameter(rho = list(prior = "gamma", shape = 2, scale = 1))

  # a name that is no advanced option is rejected, closest match suggested
  hhmodel |>
    run_inference(adapt_windows = 2L) |>
    expect_error("adapt_windows.*adapt_window")

  # the adaptation window is carried to the engine and read back
  hhmodel |>
    run_inference(steps = 4L, thin = 1L, adapt = 4L, chains = 1L,
                  adapt_window = 2L) |>
    suppressWarnings()
  expect_equal(hhmodel$setup$adapt_window, 2)
})

test_that("Progress bar renders, refreshes and closes", {
  # progressbar is time sensitive, run has to last well beyond 2 secs
  # I cannot run a long test in the CI, so I set a time limit and check that
  on.exit(setTimeLimit())
  bar <- capture.output(type = "output",
    notice <- capture.output(type = "message", { # R time-limit notice
      outcome <- tryCatch({
        setTimeLimit(elapsed = 3.5, transient = TRUE)
        hhmodel <- mk_wrapper$model |>
          add_parameter("rho", prior = "uniform") |>
          run_inference(steps = 1e10, thin = 1e4, progress = 3)
      }, interrupt = function(e) e)
    })
  )
  expect_s3_class(outcome, "interrupt")

  # an opening "0s<?s" line, then refreshes every `progress` seconds,
  # the last reporting the time left, "3s<6h17m"
  secs <- as.integer(sub(".*% ([0-9]+)s<.*", "\\1", bar))
  expect_identical(secs, c(0L, 3L))
  expect_match(bar[1], "0s<\\?s")
  expect_match(bar[length(bar)], "[0-9]+s<[0-9]+[dhms]")
})

######################################
###  EXPOSE CLASS : RISK WRAPPERS  ###
######################################

test_that("Group risk wrapper input branches", {
  hhmodel <- add_group_risk(mk_member$model, "community", by = "member")
  groups <- get_group_risks(hhmodel)
  expect_gt(nrow(groups), 0)

  # the bulk form overwrites the single one with the same definition
  expect_warning(
    hhmodel_2 <- add_group_risk(hhmodel, community = "member"),
    "[Oo]verwriting.*community"
  )
  expect_identical(get_group_risks(hhmodel_2), groups)

  # load custom definitions
  strata <- tibble::tribble(~member,   ~param_name,
                            "student", "sus_student",
                            NA,        ".reference")
  hhmodel <- hhmodel |>
    add_group_risk(
      susceptibility = list(strata = strata, params = list(sus_student = 1.23))
    )
  get_parameters(hhmodel, "sus_student")$value |> expect_equal(1.23)

  # check for non-existing custom covariates
  strata <- tibble::tribble(~member,   ~param_name,
                            "unknown", "sus_student",
                            NA,        ".reference")
  hhmodel |>
    add_group_risk("susceptibility", strata = strata, append = FALSE) |>
    expect_warning("[Uu]nknown.*member.*unknown")

  # append = FALSE drops the existing group risks
  expect_true("susceptibility" %in% get_group_risks(hhmodel)[["group_risk"]])
  hhmodel |>
    add_group_risk(community = list(by = "member"), append = FALSE) |>
    get_group_risks() |>
    getElement("group_risk") |>
    expect_setequal("community")
})

test_that("Contact risk wrapper input branches", {
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "student",    "student",  "ctr_student",
    NA,           NA,         ".reference"
  )

  # 'by' and 'strata' together is ambiguous
  mk_member$model |>
    add_contact_risk(by = "member", strata = strata) |>
    expect_error("[Ii]mpossible")

  # neither 'by' nor 'strata' is empty
  mk_member$model |>
    add_contact_risk() |>
    expect_error("[Ii]nvalid")

  # the wrapper returns a new model, leaving its input untouched
  hhmodel <- mk_member$model |>
    add_contact_risk(strata = strata, params = list(ctr_student = 1.21))
  expect_equal(nrow(mk_member$model$contact_risks), 0L)
  expect_gt(nrow(get_contact_risks(hhmodel)), 0L)
  expect_equal(get_contact_risks(hhmodel), hhmodel$contact_risks)
  # the rebuilt parameter is exposed on the returned model
  expect_true("ctr_student" %in% get_parameters(hhmodel)$name)

  # second call warns about overwrite
  hhmodel |>
    add_contact_risk(strata = strata, params = list(ctr_student = 1.21)) |>
    expect_warning("[Oo]verwriting")
})
