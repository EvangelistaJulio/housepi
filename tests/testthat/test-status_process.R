library(testthat)

####  Helpers  ####

# model carrying the data columns the definition blocks refer to,
# its member escaped: the status process labels the statuses
mk_latent <- local({
  hh_model <- data.frame(hh_id = "h1", status = 0, datum = 1,
                         record = 3, datum_2 = 2) |>
    new_housepi_model("test")
  list(model = hh_model)
})

# saturated outbreak: every household infects all of its members,
# the status process draws the status, so the data carries none
mk_outbreak <- local({
  hh_data <- sarscov2_2021_icofs |>
    dplyr::filter(hh_id %in% utils::head(unique(hh_id), 50)) |>
    dplyr::mutate(end_date = 21) |>
    subset(select = -status)
  hh_model <- hh_data |>
    new_housepi_model("case-ascertained") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(alpha = 1e-12, beta = 100, gamma = 1, hh_size_ref = 4)
  list(data = hh_data, model = hh_model)
})

# two infected statuses under the template definitions
mk_status <- local({
  hh_data <- data.frame(id = c("m1", "m2", "m3", "m4"), hh_id = "h1",
                        status = c("symptomatic", "asymptomatic",
                                   "escaped", "escaped"),
                        end_date = 21)
  hh_model <- hh_data |>
    new_housepi_model("case-ascertained") |>
    add_parameter(import = "template") |>
    add_latent(import = "template")
  list(data = hh_data, model = hh_model)
})

# every case_ascertained parameter pinned
mk_pinned <- local({
  hh_model <- new_housepi_model(model = "case-ascertained") |>
    add_parameter(import = helper_onset_house$param) |>
    add_latent(infec_time = infec_date) |>
    add_infectivity_profile(func = "gamma", shape = shape, scale = scale)
  list(model = hh_model)
})

####  Definitions  ####

test_that("A status process needs probabilities, within the status storage", {
  mk_latent$model |>
    add_status_process() |>
    expect_error("[Nn]o status probability")
  mk_latent$model |>
    add_status_process(symptomatic = 0.5, symptomatic = 0.5) |>
    expect_error("[Rr]epeats rule at row 2")

  # the engine stores at most 254 infected statuses
  weights <- function(k) {
    stats::setNames(as.list(rep(1 / k, k)), paste0("unit_", seq_len(k)))
  }
  do.call(add_status_process, c(list(mk_latent$model), weights(255))) |>
    expect_error("255 statuses.*maximum of 254")
  do.call(add_status_process, c(list(mk_latent$model), unname(weights(255)))) |>
    expect_error("255 statuses.*maximum of 254")
})

test_that("Status labels name every probability or none", {
  mk_latent$model |>
    add_status_process(symptomatic = 0.7, 0.3) |>
    expect_error("[Ss]tatus labels were given to some probabilities only")
})

test_that("The status process reads back one row per rule", {
  # inline probabilities: one rule per status, the names label them
  status <- mk_latent$model |>
    add_status_process(symptomatic = 0.7, asymptomatic = 0.3) |>
    get_status_process()
  status$status |> expect_identical(c("symptomatic", "asymptomatic"))
  status$weight |> expect_equal(c(0.7, 0.3))
  status$n |> expect_identical(c(1L, 1L))
  status$rule_group |> expect_identical(c(1L, 1L))

  # unnamed probabilities draw integer statuses
  status <- mk_latent$model |>
    add_status_process(0, 1, 0) |>
    get_status_process()
  status$status |> expect_identical(1:3)
  status$weight |> expect_equal(c(0, 1, 0))

  # a table reads back its covariates
  strata <- data.frame(datum = c(1, NA),
                       status = c("symptomatic", "asymptomatic"),
                       weight = c(2, 1))
  status <- mk_latent$model |>
    add_status_process(strata = strata) |>
    get_status_process()
  colnames(status) |>
    expect_identical(c("status", "weight", "n", "rule_group", "datum"))
  status$datum |> expect_identical(c(1, NA))
  # a read-back table declares the same status process
  mk_latent$model |>
    add_status_process(strata = status) |>
    get_status_process() |>
    expect_identical(status)

  # no status process declared
  status <- get_status_process(mk_latent$model)
  status |> expect_nrow(0L)
  colnames(status) |> expect_identical(c("status", "weight", "n", "rule_group"))
})

test_that("A status table keys every rule group on status", {
  mk_latent$model |>
    add_status_process(strata = data.frame(status = "symptomatic")) |>
    expect_error("[Tt]able has no column.*weight")
  mk_latent$model |>
    add_status_process(strata = data.frame(status = TRUE, weight = 1)) |>
    expect_error("[Cc]olumn.*status.*must hold status labels or whole numbers")
  mk_latent$model |>
    add_status_process(strata = data.frame(status = 1.5, weight = 1)) |>
    expect_error("[Cc]olumn.*status.*must hold status labels or whole numbers")

  # the datum rules scale both statuses alike
  strata <- data.frame(datum = c(NA, NA, 1), status = c("a", "b", NA),
                       weight = c(1, 2, 3), datum_2 = c(NA, NA, 2))
  mk_latent$model |>
    add_status_process(strata = strata) |>
    expect_error("[Rr]ule group of the status process does not key on.*status")

  mk_latent$model |>
    add_status_process(symptomatic = 1, strata = data.frame(status = "a",
                                                            weight = 1)) |>
    expect_error("[Mm]ixing a status table with inline status probabilities")
})

test_that("A status table is checked against the data", {
  # the single member carries datum 1
  mk_latent$model |>
    add_status_process(strata = data.frame(datum = 2, status = "a",
                                           weight = 1)) |>
    expect_error("[Ll]eaves 1 member unmatched")
  mk_latent$model |>
    add_status_process(strata = data.frame(datum = c(1, 3), status = "a",
                                           weight = 1)) |>
    expect_error("1 rule no member resolves to.*datum.*3")

  # at datum 1 and status a, two rules of one group key one column each
  strata <- data.frame(datum = c(1, NA, 1), status = c(NA, "a", "b"),
                       weight = 1)
  mk_latent$model |>
    add_status_process(strata = strata) |>
    expect_error("[Tt]ied rules.*datum.*status")

  # a factor status reads as its labels
  strata <- data.frame(datum = 1, status = factor(c("a", "b")), weight = 0)
  mk_latent$model |>
    add_status_process(strata = strata) |>
    expect_error("zero weight for 1 member")
})

test_that("Statuses are integer codes or labels throughout", {
  err <- "[Bb]oth as integer codes and as labels"
  hh_data <- data.frame(hh_id = "h1", status = c(1, 2, 0), datum = 1)
  # integer data attached last, beside status process labels
  new_housepi_model(model = "final-size") |>
    add_status_process(symptomatic = 0.5, asymptomatic = 0.5) |>
    add_dataset(hh_data) |>
    expect_error(err)
  # integer data beside a latent label
  hh_data |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform", min = 0, max = 1,
                             status = "symptomatic")) |>
    expect_error(err)
  # integer data beside a group-rate reference label
  hh_data |>
    new_housepi_model("final-size") |>
    add_group_rate(infectivity = list(
      by = "status", reference = list(status = "symptomatic")
    )) |>
    expect_error(err)

  # labelled data beside integer statuses
  hh_data_2 <- transform(hh_data, status = c("symptomatic", "asymptomatic",
                                             "escaped"))
  hh_data_2 |>
    new_housepi_model("final-size") |>
    add_status_process(0.5, 0.5) |>
    expect_error(err)
  hh_data_2 |>
    new_housepi_model("final-size") |>
    add_group_rate("infectivity", strata = data.frame(
      status = c(1, 2), param_name = c("inf_1", "inf_2")
    )) |>
    expect_error(err)

  # definitions mixing integer statuses and labels, the data carrying none
  data.frame(hh_id = "h1", status = 0, datum = 1) |>
    new_housepi_model("test") |>
    add_status_process(symptomatic = 0.5, asymptomatic = 0.5) |>
    add_latent(hidden = list(model = "uniform", min = 0, max = 1,
                             status = 1)) |>
    expect_error(err)
})

test_that("Status process labels are labels of the data", {
  hh_data <- data.frame(hh_id = "h1", status = c("symptomatic", "escaped"),
                        datum = 1)
  # each close data label is suggested for the label it corrects
  hh_data |>
    new_housepi_model("final-size") |>
    add_status_process(symptomatc = 0.5, asymptomatic = 0.5) |>
    expect_error(paste0("[Ss]tatus labels .symptomatc. and .asymptomatic. ",
                        "are not in the data.*[Dd]id you mean .symptomatic. ",
                        "for .symptomatc."))

  # a factor declares its levels, carried by a member or not
  transform(hh_data, status = factor(
    status, levels = c("escaped", "symptomatic", "asymptomatic")
  )) |>
    new_housepi_model("final-size") |>
    add_status_process(symptomatic = 0.5, asymptomatic = 0.5) |>
    get_status_process() |>
    getElement("status") |>
    expect_identical(c("symptomatic", "asymptomatic"))
})

test_that("Status probabilities are non-negative numbers", {
  # a parameter, a latent or a data column is not a number
  mk_status$model |>
    add_status_process(symptomatic = alpha, asymptomatic = 3) |>
    expect_error("[Ss]tatus probability .symptomatic. is not a number")
  mk_status$model |>
    add_status_process(symptomatic = infec_time, asymptomatic = 3) |>
    expect_error("[Ss]tatus probability .symptomatic. is not a number")
  mk_latent$model |>
    add_status_process(symptomatic = datum, asymptomatic = 3) |>
    expect_error("[Ss]tatus probability .symptomatic. is not a number")

  mk_latent$model |>
    add_status_process(symptomatic = -1, asymptomatic = 3) |>
    expect_error("[Ss]tatus weights must be non-negative")
  # some status must be drawable
  mk_latent$model |>
    add_status_process(symptomatic = 0, asymptomatic = 0) |>
    expect_error("zero weight for 1 member")
})

test_that("Status is a constant in inference", {
  # refined statuses already in the data; infec_time is read from data so the
  # density is deterministic and the only difference is the status process.
  hh_model <- helper_onset_house$data |>
    transform(status = c("symptomatic", "asymptomatic", "symptomatic",
                         "escaped", "escaped", "escaped", "escaped"),
              infec_date = infec_time) |>
    add_dataset(hh_model = mk_pinned$model)

  # the status process drops out of the density (status is fixed data)
  hh_model |>
    add_status_process(symptomatic = 0.7, asymptomatic = 0.3) |>
    evaluate_log_lik() |>
    expect_equal(evaluate_log_lik(hh_model))
})

####  End-to-end  ####

test_that("Status process generates a labelled status factor in simulation", {
  simulated <- mk_outbreak$model |>
    # categorical drawn first; labels (declaration order) name the statuses
    add_status_process(symptomatic = 0, asymptomatic = 1) |>
    run_simulation()

  # status is a factor: code 0 -> escaped, 1 -> symptomatic, 2 -> asymptomatic
  simulated$status |> expect_s3_class("factor")
  levels(simulated$status) |>
    expect_identical(c("escaped", "symptomatic", "asymptomatic"))
  # the second label is certain: a swapped label order or a draw stuck on the
  # first code labels some infected member symptomatic
  as.character(simulated$status[simulated$status != "escaped"]) |>
    expect_all_equal("asymptomatic")

  # feed the factor-status output back for inference: factorization maps the
  # labels to the same codes (escaped -> 0), so the density is finite
  hh_model <- simulated |>
    dplyr::select(id, hh_id, status, infec_date = infec_time, end_date) |>
    add_dataset(hh_model = mk_pinned$model)
  evaluate_log_lik(hh_model) |> expect_all_finite()

  # the data-derived dictionary matches the simulated infected labels
  levels(get_dataset(hh_model)$status) |>
    expect_identical(c("escaped", "symptomatic", "asymptomatic"))
})

test_that("A factor status of the status process orders the drawn labels", {
  strata <- data.frame(status = factor(c("asymptomatic", "symptomatic"),
                                       levels = c("symptomatic",
                                                  "asymptomatic")),
                       weight = c(1, 0))
  hh_model <- add_status_process(mk_outbreak$model, strata = strata)
  # the rule table reads back as given
  get_status_process(hh_model)$status |> expect_identical(strata$status)

  simulated <- run_simulation(hh_model)
  levels(simulated$status) |>
    expect_identical(c("escaped", "symptomatic", "asymptomatic"))
  as.character(simulated$status[simulated$status != "escaped"]) |>
    expect_all_equal("asymptomatic")
})

test_that("Unnamed status probabilities draw integer status codes", {
  simulated <- mk_outbreak$model |>
    add_status_process(0, 1, 0, 0) |>
    run_simulation()

  # the second of four statuses is certain, and none is labelled
  infected <- simulated$status[simulated$status != 0]
  infected |> expect_all_equal(2)
  simulated$status |> expect_type("integer")
})

test_that("A status table draws the status of each stratum", {
  # adults are symptomatic, children asymptomatic: the catch-all weighs 0
  strata <- data.frame(age = c("adult", "child", NA),
                       status = c("symptomatic", "asymptomatic", NA),
                       weight = c(1, 1, 0))
  simulated <- mk_outbreak$model |>
    add_status_process(strata = strata) |>
    run_simulation()
  status <- as.character(simulated$status)
  status[simulated$age == "adult"] |> expect_all_equal("symptomatic")
  status[simulated$age == "child"] |> expect_all_equal("asymptomatic")
})

test_that("Rule groups of a status table multiply", {
  # neither covariate alone settles the status, their product does
  strata <- data.frame(
    age = rep(c("adult", "child", NA, NA), each = 3),
    sex = rep(c(NA, NA, "F", "M"), each = 3),
    status = rep(c("symptomatic", "asymptomatic", "mild"), 4),
    weight = c(1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0)
  )
  simulated <- mk_outbreak$model |>
    add_status_process(strata = strata) |>
    run_simulation()
  status <- as.character(simulated$status)
  cell <- paste(simulated$age, simulated$sex)
  status[cell == "adult F"] |> expect_all_equal("symptomatic")
  status[cell == "adult M"] |> expect_all_equal("asymptomatic")
  status[cell == "child F"] |> expect_all_equal("mild")
  status[cell == "child M"] |> expect_all_equal("asymptomatic")
})

test_that("A status process draws statuses without any latent definition", {
  simulated <- mk_outbreak$data |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = 100, beta = 100, gamma = 1, hh_size_ref = 4) |>
    add_status_process(symptomatic = 0.5, asymptomatic = 0.5) |>
    run_simulation()

  levels(simulated$status) |>
    expect_identical(c("escaped", "symptomatic", "asymptomatic"))
  # saturated community: every member is infected, at either drawn status
  as.character(simulated$status) |>
    expect_in(c("symptomatic", "asymptomatic"))
  unique(as.character(simulated$status)) |> expect_length(2L)
})
