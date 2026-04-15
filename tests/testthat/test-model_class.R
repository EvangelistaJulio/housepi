library(testthat)

test_that("Model class inconsistent calls", {
  # Default blank
  expect_warning(new_housepi_model(), regexp = "Falling back to")

  # blank household model object
  hhmodel <- new_housepi_model("blank")

  # inference without definitions
  expect_error(hhmodel$run_inference(), regexp = "Nothing to sample")
})

test_that("Model class blank calls", {
  # blank household model object
  hhmodel <- new_housepi_model("blank")
  snapshot <- hhmodel$export_house_snapshot("house")
  expect_equal(length(snapshot$param), 0)
  expect_equal(length(snapshot$augm), 0)
  expect_equal(length(snapshot$group_rate), 0)

  # blank default data wide format
  expect_equal(hhmodel$data_wide, helper_test_get_data("blank_house"))
  expect_equal(snapshot$data_wide, hhmodel$data_wide)

  # correct data long format
  data_long <- data.frame(id = snapshot$data_wide$id[1],
                          time = seq(3), measure = seq(3)**3)
  hhmodel$add_data_long(data_long)
  expect_equal(hhmodel$data_long, data_long)

  # factorization
  complex_data <- helper_test_get_data("sarscov2_2021_icofs")
  complex_data$dates <- as.Date(runif(dim(complex_data)[1], max = 100))
  hhmodel$add_data_wide(complex_data)
  expect_equal(hhmodel$data_wide, complex_data)

  # templates
  expect_warning(hhmodel$load_parameter_template(), regexp = "blank")
  expect_warning(hhmodel$load_augmented_template(), regexp = "blank")

  # zero log-likelihood
  expect_warning(x <- hhmodel$evaluate_log_lik(), regexp = "blank")
  expect_equal(x, 0)

  # blank simulation
  hhmodel <- new_housepi_model("blank")
  hhmodel$add_augmented(name = "dummy", value = 3)
  expect_warning(out <- hhmodel$run_simulation(), regexp = "blank")
  expect_true(
    all(sapply(out$dummy[out$status == 0], function(x) identical(x, NA_real_)))
  )
})

test_that("Rlang calls work as expected", {
  hhmodel <- new_housepi_model("blank")

  val <- 1
  hhmodel$load_parameter_list(
    var1 = list(prior = "normal", mean = indiv_shape + 1),
    var2 = list(prior = "normal", mean = indiv_shape + !!val),
    var3 = list(prior = "normal", mean = ~ indiv_shape + 1),
    var4 = list(prior = "normal", mean = ~ indiv_shape + !!val)
  )
  aux <- hhmodel$parameters
  aux$name <- NULL
  expect_equal(nrow(unique(aux)), 1)

  val <- 3
  hhmodel$load_augmented_list(
    var1 = list(model = "normal", mean = indiv_shape + 3),
    var2 = list(model = "normal", mean = indiv_shape + !!val),
    var3 = list(model = "normal", mean = ~ indiv_shape + 3),
    var4 = list(model = "normal", mean = ~ indiv_shape + !!val)
  )
  aux <- hhmodel$augmented
  aux$name <- NULL
  expect_equal(nrow(unique(aux)), 1)
})

test_that("Discretization is enabled as expected", {
  # blank household model object
  hhmodel <- new_housepi_model("blank")
  hhmodel$load_parameter_list(
    sample = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE),
    uniform = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                   proposal = "uniform"),
    normal = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                  proposal = "normal"),
    lognormal = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                     proposal = "lognormal")
  )

  results <- hhmodel$run_inference()
  y <- posterior::as_draws_matrix(get_param_draws(results))
  expect_true(all(y == floor(y)))
})

test_that("Augmented aliasing works as expected", {
  # check that model arguments are correctly integrated
  hhmodel <- new_housepi_model("blank") |>
    add_data_wide(data.frame(id = seq(1, 1000), house_id = 1, status = 1)) |>
    load_augmented_list(
      hour = list(model = "uniform", min = day, max = day + 1),
      day = list(model = "gamma", shape = 2, scale = 1)
    )
  simulated <- run_simulation(hhmodel)
  expect_true(all(simulated$day <= simulated$hour &
                    simulated$hour <= simulated$day + 1))
})

test_that("Nested definitions works as expected", {
  hhmodel <- new_housepi_model("blank")

  # define Box-Muller transformation
  hhmodel$load_augmented_list(
    norm = list(model = "normal"),
    u1 = list(model = "uniform"),
    u2 = list(model = "uniform"),
    z0 = sqrt(- 2 * log(u1)) * cos(2 * pi * u2)
  )

  results <- hhmodel$run_inference(steps = 1e5, progress = 0)
  sample_norm <- results$augm_draws$tommy_norm
  sample_boxmul <- results$augm_draws$tommy_z0

  expect_equal(mean(sample_norm), 0, tolerance = 0.05)
  expect_equal(sd(sample_norm), 1, tolerance = 0.05)
  probs <- c(0.025, 0.5, 0.975)
  q <- setNames(stats::quantile(sample_norm, probs), NULL)
  expect_equal(q, qnorm(probs), tolerance = 0.05)

  expect_equal(mean(sample_boxmul), 0, tolerance = 0.05)
  expect_equal(sd(sample_boxmul), 1, tolerance = 0.05)
  q <- setNames(stats::quantile(sample_boxmul, probs), NULL)
  expect_equal(q, qnorm(probs), tolerance = 0.05)
})

test_that("Load augmented list works fine for external lists", {
  hhmodel <- new_housepi_model("blank")

  external_list <- list(z0 = ~ date@D, z1 = ~ sqrt(log(z0)))
  # raise error for on injection
  expect_error(hhmodel$load_augmented_list(external_list),
               regexp = "injection operator")

  tmp <- hhmodel$load_augmented_list(!!!external_list)$augmented
  hhmodel$load_augmented_list(z0 = date@D, z1 = sqrt(log(z0)))
  expect_equal(tmp, hhmodel$augmented)
})

test_that("Augmented traces export only relevant variables", {
  hhmodel <- new_housepi_model("blank")
  infec_df <- subset(hhmodel$data_wide, status > 0)
  hhmodel$load_parameter_list(alpha = 1, beta = list(prior = "uniform"))
  hhmodel$load_augmented_list(
    # constants
    blocked = 42,
    # private variables
    .hidden = list(model = "uniform"),
    .secret = list(model = "normal"),
    # public non constant function
    layer = .secret - .hidden,
    # equal to parameters
    alpha = alpha@P,
    alpha_beta = alpha@P + beta@P,
    # equal to other augmented objects
    dummy = blocked,
    alias = layer,
    # implicit function
    impl = list(func = "gamma", func_args = c(shape = 3)),
    # constant function
    invar = list(func = "gamma", func_args = c(x = dummy, shape = 3)),
    # function
    expl = list(func = "gamma", func_args = c(x = alias, shape = 3))
  )
  result <- hhmodel$run_inference(id_draws = infec_df$id)
  expect_equal(
    posterior::variables(result$augm_draws),
    paste(infec_df$id, c("layer", "alias", "alpha_beta", "expl"), sep = "_")
  )
})

test_that("Custom models shares function registries", {
  hhmodel_0 <- new_housepi_model("blank")
  expect_error(expose_function("growth"), regexp = "Unknown function")

  custom <- system.file(file.path("examples", "custom"), package = "housepi")
  hhmodel_1 <- new_housepi_model(custom, quiet = TRUE)
  expect_no_error(expose_function("growth")$func(seq(3)))
})

test_that("Modules do not get overwritten", {
  # NOTE: I use profile as model status flag
  hhmodel_1 <- new_housepi_model("augmented_infection_time")
  out_1 <- hhmodel_1$infec_profile
  expect_gt(length(out_1), 0)

  hhmodel_2 <- new_housepi_model("chain_binomial")
  out_2 <- hhmodel_2$infec_profile
  expect_length(out_2, 0)

  # reload first model base Model again
  hhmodel_3 <- new_housepi_model("augmented_infection_time")
  out_3 <- hhmodel_3$infec_profile
  expect_identical(out_1, out_3)
})

test_that("Print behavior works as expected", {
  # with data and parameters
  hhmodel <- new_housepi_model("chain_binomial")
  hhmodel$add_data_wide(helper_test_get_data("flu_1918_frederick"))
  hhmodel$load_parameter_template()
  expect_snapshot(print(hhmodel))

  # without parameters
  hhmodel2 <- new_housepi_model("chain_binomial")
  hhmodel2$add_data_wide(helper_test_get_data("flu_1918_frederick"))
  expect_snapshot(print(hhmodel2))

  # blank model
  expect_snapshot(print(new_housepi_model("blank")))
})

test_that("Active fields are read-only", {
  hhmodel <- new_housepi_model("chain_binomial")
  hhmodel$add_data_wide(helper_test_get_data("flu_1918_frederick"))
  expect_error(hhmodel$data_wide <- "La", regexp = "read only")
  expect_error(hhmodel$data_long <- "Li", regexp = "read only")
  expect_error(hhmodel$parameters <- "Lu", regexp = "read only")
  expect_error(hhmodel$augmented <- "Le", regexp = "read only")
  expect_error(hhmodel$infec_profile <- "Lo", regexp = "read only")
  expect_error(hhmodel$group_rates <- "The", regexp = "read only")
  expect_error(hhmodel$name <- "Patriots", regexp = "read only")
  expect_error(hhmodel$setup <- "Cipher", regexp = "read only")
})

###########################################################################
# Contact rate snapshot tests
###########################################################################

test_that("Snapshot without contact rates returns empty matrix", {
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  data_wide <- helper_test_get_data("sarscov2_2021_icofs")
  hhmodel$add_data_wide(data_wide)
  hhmodel$load_parameter_template()

  hid <- unique(data_wide$house_id)[1]
  snapshot <- hhmodel$export_house_snapshot(hid)

  # no contact rate defs -> empty matrix
  cr_mtx <- snapshot$contact_rate
  expect_true(is.matrix(cr_mtx))
  expect_equal(nrow(cr_mtx), 0)
  expect_equal(ncol(cr_mtx), 0)
})

test_that("Snapshot respects constant contact rate definitions", {
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  data_wide <- helper_test_get_data("sarscov2_2021_icofs")
  hhmodel$add_data_wide(data_wide)
  hhmodel$load_parameter_template()

  # define contact rates by sex with reference M->M (symmetric by default)
  hhmodel$add_contact_rate(
    by = "sex",
    reference = c(from_sex = "M", to_sex = "M"),
    prefix = "cr",
    prior = "lognormal", sdlog = 1
  )
  # override as constants for deterministic checking
  # symmetric: M->F and F->M share param cr_sex_F_to_sex_M
  hhmodel$load_parameter_list(
    cr_sex_F_to_sex_M = 2.0,
    cr_sex_F_to_sex_F = 1.5,
    append = TRUE
  )

  # pick a household with mixed sexes
  house_ids <- unique(data_wide$house_id)
  mixed <- NULL
  for (hid in house_ids) {
    sub <- data_wide[data_wide$house_id == hid, ]
    sexes <- unique(sub$sex)
    if (length(sexes) > 1 && nrow(sub) >= 2) {
      mixed <- hid
      break
    }
  }
  expect_false(is.null(mixed))

  snapshot <- hhmodel$export_house_snapshot(mixed)
  cr_mtx <- snapshot$contact_rate
  snap_wide <- snapshot$data_wide
  n <- nrow(snap_wide)

  expect_true(is.matrix(cr_mtx))
  expect_equal(nrow(cr_mtx), n)
  expect_equal(ncol(cr_mtx), n)

  # diagonal should be reference (1.0)
  expect_equal(unname(diag(cr_mtx)), rep(1, n))

  # M->M off-diagonal pairs: reference (1.0)
  m_ids <- snap_wide$id[snap_wide$sex == "M"]
  for (i in m_ids) for (j in m_ids) {
    expect_equal(cr_mtx[i, j], 1.0)
  }

  # symmetric: M->F and F->M should both be 2.0
  f_ids <- snap_wide$id[snap_wide$sex == "F"]
  if (length(m_ids) > 0 && length(f_ids) > 0) {
    expect_equal(cr_mtx[m_ids[1], f_ids[1]], 2.0)
    expect_equal(cr_mtx[f_ids[1], m_ids[1]], 2.0)
  }
  # F->F should be 1.5
  if (length(f_ids) > 1) {
    expect_equal(
      cr_mtx[f_ids[1], f_ids[2]], 1.5
    )
  }
})

test_that("Snapshot with multi-dim contact rates (joint)", {
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  data_wide <- helper_test_get_data("sarscov2_2021_icofs")
  # discretize age for a small number of levels
  data_wide$age_class <- ifelse(data_wide$age < 18, "child", "adult")
  hhmodel$add_data_wide(data_wide)
  hhmodel$load_parameter_template()

  # two features jointly: age_class x sex -> 4 from x 4 to = 16
  hhmodel$add_contact_rate(
    by = c("age_class", "sex"),
    reference = c(
      from_age_class = "adult", from_sex = "M",
      to_age_class = "adult", to_sex = "M"
    ),
    prefix = "cr",
    prior = "lognormal", sdlog = 1
  )

  # verify logic table dimensions
  ct <- hhmodel$contact_rates
  expect_equal(nrow(ct), 16)
  expect_equal(sum(ct$param_name == ".reference"), 1)

  # set all non-reference params to constant 2.0
  pnames <- setdiff(ct$param_name, ".reference")
  plist <- stats::setNames(
    as.list(rep(2.0, length(pnames))), pnames
  )
  expect_warning(
    do.call(hhmodel$load_parameter_list,
            c(plist, list(append = TRUE))),
    regexp = "Overwriting"
  )

  # pick a household with >= 2 members
  house_ids <- unique(data_wide$house_id)
  hid <- NULL
  for (h in house_ids) {
    if (sum(data_wide$house_id == h) >= 3) {
      hid <- h
      break
    }
  }
  expect_false(is.null(hid))

  snapshot <- hhmodel$export_house_snapshot(hid)
  cr_mtx <- snapshot$contact_rate
  snap_wide <- snapshot$data_wide

  expect_true(is.matrix(cr_mtx))
  expect_equal(nrow(cr_mtx), nrow(snap_wide))

  # diagonal is 1.0 (self-contact)
  expect_equal(unname(diag(cr_mtx)), rep(1, nrow(snap_wide)))

  # adult-M -> adult-M off-diag should be 1.0 (reference)
  am_ids <- snap_wide$id[
    snap_wide$age_class == "adult" & snap_wide$sex == "M"
  ]
  for (i in am_ids) for (j in am_ids) {
    expect_equal(cr_mtx[i, j], 1.0)
  }

  # any non-reference pair should be 2.0
  other_ids <- setdiff(snap_wide$id, am_ids)
  if (length(am_ids) > 0 && length(other_ids) > 0) {
    expect_equal(
      cr_mtx[am_ids[1], other_ids[1]], 2.0
    )
  }
})

test_that("Snapshot with non-joint contact rates (multiplicative)", {
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  data_wide <- helper_test_get_data("sarscov2_2021_icofs")
  data_wide$age_class <- ifelse(data_wide$age < 18, "child", "adult")
  hhmodel$add_data_wide(data_wide)
  hhmodel$load_parameter_template()

  # two features independently: 2 sub-tables of 4 rows each
  hhmodel$add_contact_rate(
    by = c("age_class", "sex"),
    joint = FALSE,
    reference = c(age_class = "adult", sex = "M"),
    prefix = "cr",
    prior = "lognormal", sdlog = 1
  )

  ct <- hhmodel$contact_rates
  expect_equal(nrow(ct), 8) # 4 + 4

  # set age non-reference params to 2.0 and sex to 3.0
  age_params <- setdiff(
    ct$param_name[!is.na(ct$from_age_class)], ".reference"
  )
  sex_params <- setdiff(
    ct$param_name[!is.na(ct$from_sex)], ".reference"
  )
  plist <- c(
    stats::setNames(as.list(rep(2.0, length(age_params))),
                    age_params),
    stats::setNames(as.list(rep(3.0, length(sex_params))),
                    sex_params)
  )
  expect_warning(
    do.call(hhmodel$load_parameter_list,
            c(plist, list(append = TRUE))),
    regexp = "Overwriting"
  )

  # pick a household with >= 3 members
  house_ids <- unique(data_wide$house_id)
  hid <- NULL
  for (h in house_ids) {
    if (sum(data_wide$house_id == h) >= 3) {
      hid <- h
      break
    }
  }
  expect_false(is.null(hid))

  snapshot <- hhmodel$export_house_snapshot(hid)
  cr_mtx <- snapshot$contact_rate
  snap_wide <- snapshot$data_wide

  expect_true(is.matrix(cr_mtx))
  expect_equal(unname(diag(cr_mtx)), rep(1, nrow(snap_wide)))

  # adult-M -> adult-M: both reference -> 1.0
  am_ids <- snap_wide$id[
    snap_wide$age_class == "adult" & snap_wide$sex == "M"
  ]
  for (i in am_ids) for (j in am_ids) {
    expect_equal(cr_mtx[i, j], 1.0)
  }

  # adult-M -> child-M: age=2.0, sex=ref -> 2.0
  cm_ids <- snap_wide$id[
    snap_wide$age_class == "child" & snap_wide$sex == "M"
  ]
  if (length(am_ids) > 0 && length(cm_ids) > 0) {
    expect_equal(cr_mtx[am_ids[1], cm_ids[1]], 2.0)
  }

  # adult-M -> adult-F: age=ref, sex=3.0 -> 3.0
  af_ids <- snap_wide$id[
    snap_wide$age_class == "adult" & snap_wide$sex == "F"
  ]
  if (length(am_ids) > 0 && length(af_ids) > 0) {
    expect_equal(cr_mtx[am_ids[1], af_ids[1]], 3.0)
  }

  # adult-M -> child-F: age=2.0 * sex=3.0 -> 6.0
  cf_ids <- snap_wide$id[
    snap_wide$age_class == "child" & snap_wide$sex == "F"
  ]
  if (length(am_ids) > 0 && length(cf_ids) > 0) {
    expect_equal(cr_mtx[am_ids[1], cf_ids[1]], 6.0)
  }
})

test_that("Add data type errors are conserved", {
  hhmodel <- new_housepi_model("blank")
  expect_error(hhmodel$add_data_wide(data.frame()), regexp = "empty")
  expect_error(hhmodel$add_data_wide("ground"), regexp = "data.frame")
  expect_error(hhmodel$add_data_long("zero"), regexp = "data.frame")
})
