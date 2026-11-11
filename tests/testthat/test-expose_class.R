library(testthat)

####  Helpers  ####

# model with no data of its own beyond the 'test' default
mk_wrapper <- local({
  list(model = new_housepi_model(model = "test"))
})

# households over a member covariate the rate wrappers stratify on
mk_member <- local({
  hh_model <- data.frame(
    hh_id = c(1, 1, 2, 2, 2),
    status = c(1, 0, 1, 0, 0),
    member = c("student", "other", "student", "mother", "father")
  ) |>
    new_housepi_model("final-size")
  list(model = hh_model)
})

####  Dispatch  ####

test_that("Wrappers reject objects that are not a housepi model", {
  obj <- "not a model"
  err <- "must be of class"
  add_long_dataset(obj, data.frame()) |> expect_error(err)
  add_dataset(obj, data.frame()) |> expect_error(err)
  get_dataset(obj) |> expect_error(err)
  get_long_dataset(obj) |> expect_error(err)
  get_parameters(obj) |> expect_error(err)
  get_latent(obj) |> expect_error(err)
  get_group_rates(obj) |> expect_error(err)
  get_infectivity_profiles(obj) |> expect_error(err)
  get_reporting_process(obj) |> expect_error(err)
  add_contact_rate(obj) |> expect_error(err)
  add_parameter(obj) |> expect_error(err)
  add_latent(obj) |> expect_error(err)
  add_infectivity_profile(obj, "gamma") |> expect_error(err)
  add_reporting_process(obj, dataset = identity) |> expect_error(err)
  evaluate_log_lik(obj) |> expect_error(err)
  run_inference(obj) |> expect_error(err)
  run_simulation(obj) |> expect_error(err)
  add_group_rate(obj) |> expect_error(err)
  check_predictive(obj, NULL) |> expect_error(err)
})

test_that("Data provided at construction is exposed via get_dataset", {
  hh_data <- get_dataset(new_housepi_model(flu_1918_baltimore, "test"))
  colnames(hh_data) |> expect_contains(c("hh_id", "status"))

  # id automatically added
  hh_data$id |> expect_identical(as.character(seq_len(nrow(hh_data))))
  # hh_id cast to character
  hh_data$hh_id |> expect_identical(as.character(flu_1918_baltimore$hh_id))
  # statuses are integer
  hh_data$status |> expect_identical(flu_1918_baltimore$status)
})

test_that("A getter refusing a name points to the close one, or to all", {
  hh_model <- mk_wrapper$model |> add_parameter(rho = 1)
  get_parameters(hh_model, "rhoo") |>
    expect_error("[Nn]o definition named .rhoo.*[Dd]id you mean .rho")
  get_parameters(hh_model, "unknown") |>
    expect_error("[Nn]o definition named .unknown.*[Oo]mit .name.")
})

test_that("Variadic labels are not swallowed by formal matching", {
  # single-letter categorical weights c/d/p would partial-match the
  # const/discr/prior formals of add_parameter; they reach the prior call
  param <- mk_wrapper$model |>
    add_parameter("rho", prior = "categorical", a = 0.2, b = 0.1, c = 0.7) |>
    get_parameters()
  param$prior |> expect_identical("categorical")
  param$args[[1]] |> expect_equal(c(a = 0.2, b = 0.1, c = 0.7))

  # bulk form
  param <- mk_wrapper$model |>
    add_parameter(
      rho = list(prior = "categorical", a = 0.2, b = 0.1, c = 0.7),
      rho_2 = list(prior = "categorical", d = 0.4, p = 0.6)
    ) |>
    get_parameters()
  param$name |> expect_identical(c("rho", "rho_2"))
  param$prior |> expect_identical(c("categorical", "categorical"))
  param$args[[1]] |> expect_equal(c(a = 0.2, b = 0.1, c = 0.7))
  param$args[[2]] |> expect_equal(c(d = 0.4, p = 0.6))
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

  # so are the arguments only single mode reads
  mk_wrapper$model |>
    add_parameter(import = defs, max = 2) |>
    expect_error("[Mm]ixing")

  # the template sentinel follows the same exclusion rules
  mk_wrapper$model |>
    add_parameter("rho_3", import = "template") |>
    expect_error("[Mm]ixing")
  mk_wrapper$model |>
    add_latent(hidden = 1, import = "template") |>
    expect_error("[Mm]ixing")
  mk_wrapper$model |>
    add_latent(import = "template", discr = TRUE) |>
    expect_error("[Mm]ixing")

  # a character import naming no known source is a typo, not a definition
  mk_wrapper$model |>
    add_parameter(import = "templte") |>
    expect_error("[Uu]nknown import source")
})

test_that("Status probabilities import as a list or from the template", {
  # its member escaped: the status process labels the statuses
  hh_model <- new_housepi_model(data.frame(hh_id = "h1", status = 0), "test")
  # an external list equals the inline probabilities, labelled or not
  hh_model |>
    add_status_process(import = list(symptomatic = 0.7, asymptomatic = 0.3)) |>
    get_status_process() |>
    expect_identical(get_status_process(add_status_process(
      hh_model, symptomatic = 0.7, asymptomatic = 0.3
    )))
  hh_model |>
    add_status_process(import = list(0, 1)) |>
    get_status_process() |>
    expect_identical(get_status_process(add_status_process(hh_model, 0, 1)))

  # an imported table equals the `strata` argument
  strata <- data.frame(status = c("symptomatic", "asymptomatic"),
                       weight = c(0.7, 0.3))
  hh_model |>
    add_status_process(import = list(strata = strata)) |>
    get_status_process() |>
    expect_identical(get_status_process(add_status_process(
      hh_model, strata = strata
    )))
  hh_model |>
    add_status_process(strata = strata, import = list(strata = strata)) |>
    expect_error("[Mm]ixing imported and inline status probabilities")

  # `import` cannot be combined with inline probabilities
  mk_wrapper$model |>
    add_status_process(symptomatic = 1, import = list(asymptomatic = 1)) |>
    expect_error("[Mm]ixing imported and inline status probabilities")
  mk_wrapper$model |>
    add_status_process(symptomatic = 1, import = "template") |>
    expect_error("[Mm]ixing imported and inline status probabilities")
  mk_wrapper$model |>
    add_status_process(import = "templte") |>
    expect_error("[Uu]nknown import source")

  # a model shipping no status template
  new_housepi_model(model = "final-size") |>
    add_status_process(import = "template") |>
    expect_error("[Ii]mpossible to use templated status process")
  expect_warning(
    hh_model <- add_status_process(mk_wrapper$model, import = "template"),
    "[Nn]o template status process definitions.*test model"
  )
  get_status_process(hh_model) |> expect_nrow(0L)
})

test_that("A definition named like a single-mode argument is not dropped", {
  hh_model <- data.frame(hh_id = "h1", status = 1, datum = 2) |>
    new_housepi_model("test")
  err <- "[Dd]efinition can.t be named"

  mk_wrapper$model |>
    add_parameter(rho = 1, min = 2, rho_2 = 3) |>
    expect_error(paste0(err, ".*min"))
  hh_model |>
    add_latent(hidden = min(datum), const = 1) |>
    expect_error(paste0(err, ".*const"))
  hh_model |>
    add_observation(record = min(datum), func = max(datum),
                    level = "house") |>
    expect_error(paste0(err, ".*func"))
  mk_member$model |>
    add_group_rate(community = "member", reference = ".auto") |>
    expect_error(paste0(err, ".*reference"))
  data.frame(hh_id = "h1", status = 1, end_date = 20) |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_infectivity_profile(profile = list(func = "gamma"), discr = TRUE) |>
    expect_error(paste0(err, ".*discr"))

  # the arguments multiple mode reads apply to every definition
  hh_model |>
    add_latent(hidden = min(datum), hidden_2 = max(datum), status = 1) |>
    get_latent() |>
    getElement("status") |>
    expect_identical(c("1", "1"))
})

test_that("Single-mode arguments without the key name the missing key", {
  err <- "[Mm]issing .%s. of the single %s definition"
  mk_wrapper$model |>
    add_parameter(prior = "uniform", max = 3) |>
    expect_error(sprintf(err, "name", "parameter"))
  hh_model <- data.frame(hh_id = "h1", status = 1, end_date = 20) |>
    new_housepi_model("case-ascertained", quiet = TRUE)
  hh_model |>
    add_infectivity_profile(discr = TRUE) |>
    expect_error(sprintf(err, "func", "infectivity profile"))

  # a key set to `NA` is missing
  mk_wrapper$model |>
    add_parameter(name = NA, prior = "uniform", max = 3) |>
    expect_error(sprintf(err, "name", "parameter"))
  hh_model |>
    add_infectivity_profile(func = NA, discr = FALSE, a = 1) |>
    expect_error(sprintf(err, "func", "infectivity profile"))
})

test_that("A bulk definition built outside the call binds like the literal", {
  param <- mk_wrapper$model |>
    add_parameter(rho = list(prior = "uniform", max = 3)) |>
    get_parameters()
  def <- list(prior = "uniform", max = 3)
  mk_wrapper$model |>
    add_parameter(rho = def) |>
    get_parameters() |>
    expect_identical(param)
  mk_wrapper$model |>
    add_parameter(rho = c(list(prior = "uniform"), list(max = 3))) |>
    get_parameters() |>
    expect_identical(param)

  # a name bound to no list stays a reference
  mk_wrapper$model |>
    add_parameter(rho = 1, rho_2 = rho) |>
    get_parameters("rho_2") |>
    getElement("func") |>
    expect_identical("rho|P")

  def_rate <- list(by = "member")
  mk_member$model |>
    add_group_rate(community = def_rate) |>
    get_group_rates() |>
    expect_identical(
      get_group_rates(add_group_rate(mk_member$model, community = "member"))
    )
})

test_that("Bulk import resolves symbols in the calling frame", {
  # `datum_local` exists only inside the caller
  param <- local({
    datum_local <- 3
    mk_wrapper$model |>
      add_parameter(import = list(rho = quote(datum_local * 2))) |>
      get_parameters()
  })
  param$const |> expect_equal(6)
})

test_that("add_observation routes every input branch on `level`", {
  hh_model <- data.frame(hh_id = "h1", status = 1, record = 3, datum = 2) |>
    new_housepi_model("test")

  # default level is member
  hh_model |>
    add_observation("record", model = "normal") |>
    get_observations() |>
    getElement("level") |>
    expect_identical("member")

  # single, bulk and imported definitions
  single <- hh_model |>
    add_observation("record_2", func = min(datum), level = "house")
  bulk <- hh_model |>
    add_observation(record_2 = min(datum), level = "house")
  imported <- hh_model |>
    add_observation(import = list(record_2 = quote(min(datum))),
                    level = "house")

  get_observations(bulk) |> expect_identical(get_observations(single))
  get_observations(imported) |> expect_identical(get_observations(single))
  get_observations(single)[["level"]] |> expect_identical("house")

  # a member-level function reaches the definition through `func`
  single <- hh_model |>
    add_observation("record_2", func = min(datum), status = 1)
  bulk <- hh_model |>
    add_observation(record_2 = list(func = min(datum), status = 1))
  imported <- hh_model |>
    add_observation(
      import = list(record_2 = list(func = quote(min(datum)), status = 1))
    )

  get_observations(bulk) |> expect_identical(get_observations(single))
  get_observations(imported) |> expect_identical(get_observations(single))
  get_observations(single)[["func"]] |> expect_identical("min")

  # an unknown level is rejected before anything is stored
  hh_model |>
    add_observation(record_2 = min(datum), level = "unknown") |>
    expect_error("[Uu]nknown .level. .unknown")
  # a prefix of a level is not read as that level
  hh_model |>
    add_observation(record_2 = min(datum), level = "hous") |>
    expect_error("[Uu]nknown .level.*[Dd]id you mean .house")
  hh_model |>
    add_observation(record_2 = min(datum), level = "h") |>
    expect_error("[Uu]nknown .level.*add_observation")
})

test_that("Reporting process wrapper replaces only the slots passed", {
  # a model starts without one
  get_reporting_process(mk_wrapper$model) |>
    expect_identical(list(dataset = NULL, long_dataset = NULL))

  hh_model <- mk_wrapper$model |>
    add_reporting_process(dataset = identity, long_dataset = nrow)
  # an omitted slot is kept
  hh_model_2 <- add_reporting_process(hh_model, long_dataset = ncol)
  get_reporting_process(hh_model_2) |>
    expect_identical(list(dataset = identity, long_dataset = ncol))
  # NULL removes a slot
  add_reporting_process(hh_model_2, dataset = NULL) |>
    get_reporting_process() |>
    expect_identical(list(dataset = NULL, long_dataset = ncol))
  # the input model is left untouched
  get_reporting_process(hh_model) |>
    expect_identical(list(dataset = identity, long_dataset = nrow))
})

test_that("Reporting process wrapper rejects anything but a function", {
  mk_wrapper$model |>
    add_reporting_process() |>
    expect_error("[Nn]o reporting process given")
  mk_wrapper$model |>
    add_reporting_process(dataset = "not a function") |>
    expect_error("[Rr]eporting process .dataset. must be a function")
})

test_that("Advanced settings route through `...`", {
  hh_model <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")

  # a name that is no advanced option is rejected, closest match suggested
  hh_model |>
    run_simulation(time_steps = 0.5) |>
    expect_error("time_steps.*time_step")

  # the integration step is carried to the engine and read back
  simulated <- run_simulation(hh_model, time_step = 0.5)
  simulated |> expect_s3_class("data.frame")
  hh_model$settings$time_step |> expect_equal(0.5)
  # threads, a formal, is carried the same way
  run_simulation(hh_model, threads = 2L) |> expect_s3_class("data.frame")
  hh_model$settings$threads |> expect_equal(2)
})

test_that("Progress bar prints every `progress` seconds until interrupted", {
  # the run outlasts the time limit, which stops it: an interrupt at the
  # engine's check, an error while a print runs
  on.exit(setTimeLimit())
  # the bar prints as messages
  printed <- capture.output(type = "message", {
    outcome <- tryCatch({
      setTimeLimit(elapsed = 3.5, transient = TRUE)
      mk_wrapper$model |>
        add_parameter("rho", prior = "uniform") |>
        run_inference(steps = 1e10, thin = 1e4, progress = 1)
    }, interrupt = \(e) e, error = \(e) e)
  })
  outcome |> expect_s3_class(c("interrupt", "error"))

  # every line reports a percentage: the interrupted bar prints no closing line
  # interactive prints overwrite one line, each led by a carriage return
  bar <- grep("^sample [|]", unlist(strsplit(printed, "\r")), value = TRUE)
  length(bar) |> expect_gte(2L)
  bar |> expect_match("%")
  # elapsed whole seconds before "<", one second apart
  secs <- as.integer(sub(".*% ([0-9]+)s<.*", "\\1", bar))
  diff(secs) |> expect_all_equal(1L)
})

####  Rate wrappers  ####

test_that("Group rate wrapper input branches", {
  hh_model <- add_group_rate(mk_member$model, "community", by = "member")
  groups <- get_group_rates(hh_model)
  nrow(groups) |> expect_gt(0)

  # the bulk form overwrites the single one with the same definition
  expect_warning(
    hh_model_2 <- add_group_rate(hh_model, community = "member"),
    "[Oo]verwriting.*community"
  )
  get_group_rates(hh_model_2) |> expect_identical(groups)

  # load custom definitions
  strata <- tibble::tribble(~member,   ~param_name,
                            "student", "sus_student",
                            NA,        ".reference")
  hh_model <- hh_model |>
    add_group_rate(
      susceptibility = list(strata = strata, params = list(sus_student = 1.23))
    )
  get_parameters(hh_model, "sus_student")$const |> expect_equal(1.23)

  # a covariate value the data never carries resolves no member
  strata <- tibble::tribble(~member,   ~param_name,
                            "unknown", "sus_student",
                            NA,        ".reference")
  hh_model |>
    add_group_rate("susceptibility", strata = strata, append = FALSE) |>
    expect_error("1 rule no member resolves to.*member = .unknown.")

  # append = FALSE drops the existing group rates
  get_group_rates(hh_model)[["group_rate"]] |> expect_contains("susceptibility")
  hh_model |>
    add_group_rate(community = list(by = "member"), append = FALSE) |>
    get_group_rates() |>
    getElement("group_rate") |>
    expect_setequal("community")
})

test_that("Contact rate wrapper input branches", {
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "student",    "mother",   "ctr_student",
    NA,           NA,         ".reference"
  )

  # 'by' and 'strata' together is ambiguous
  mk_member$model |>
    add_contact_rate(by = "member", strata = strata) |>
    expect_error("[Ii]mpossible")

  # neither 'by' nor 'strata' is empty
  mk_member$model |>
    add_contact_rate() |>
    expect_error("[Ii]nvalid")

  # the wrapper returns a new model, leaving its input untouched
  hh_model <- mk_member$model |>
    add_contact_rate(strata = strata, params = list(ctr_student = 1.21))
  get_contact_rates(mk_member$model) |> expect_nrow(0L)
  nrow(get_contact_rates(hh_model)) |> expect_gt(0L)
  # the rebuilt parameter is exposed on the returned model
  get_parameters(hh_model)$name |> expect_contains("ctr_student")

  # second call warns about overwrite
  hh_model |>
    add_contact_rate(strata = strata, params = list(ctr_student = 1.21)) |>
    expect_warning("[Oo]verwriting")
})
