library(testthat)

# Two houses over three risk classes, with group and contact risks fixed at
# values the R reference below reads back: `_1` holds a single infected member,
# `_3` three of them. Class `a` is the reference class.
mk_house <- local({
  groups <- list(community = list("class", reference = "a"),
                 susceptibility = list("class", reference = "a"),
                 infectivity = list("class", reference = "a"))
  contacts <- list(by = "class", symmetric = FALSE,
                   reference = c(from_class = "a", to_class = "b"))

  grisks <- list(a = list(com = 1., sus = 1., inf = 1.),
                 b = list(com = 2., sus = 3., inf = 0.2),
                 c = list(com = 0.5, sus = 0.8, inf = 1.5))
  crisks_1 <- list(class_b_to_class_a = 2.3)
  crisks_3 <- c(crisks_1,
                list(class_a_to_class_c = 0.7,
                     class_c_to_class_a = 0.6,
                     class_b_to_class_c = 1.1,
                     class_c_to_class_b = 0.9,
                     class_c_to_class_c = 1.2))

  # members sorted infected, escaped
  house_3 <- data.frame(
    id = sprintf("m%s", seq(1, 4)),
    house_id = "h1",
    status = c(rep(1, 3), 0),
    class = c("a", "b", rep("c", 2)),
    community = c(1, grisks$b$com, rep(grisks$c$com, 2)),
    susceptibility = c(1, grisks$b$sus, rep(grisks$c$sus, 2)),
    infectivity = c(1, grisks$b$inf, rep(grisks$c$inf, 2))
  )
  house_1 <- house_3[1:2, ]
  house_1$status <- c(1, 0)

  # parameters matching automatic names
  par_1 <- c(list(alpha = 0.15, beta = 0.22, gamma = 1, hhsize_ref = 4),
             crisks_1,
             list(com_class_b = grisks$b$com,
                  sus_class_b = grisks$b$sus,
                  inf_class_b = grisks$b$inf))
  par_3 <- c(par_1, crisks_3[-1],
             list(com_class_c = grisks$c$com,
                  sus_class_c = grisks$c$sus,
                  inf_class_c = grisks$c$inf))

  list(groups = groups, contacts = contacts,
       house_1 = house_1, house_3 = house_3,
       crisks_1 = crisks_1, crisks_3 = crisks_3,
       par_1 = par_1, par_3 = par_3)
})

# contact risk matrix of a house, the reference reading 1
make_ctr_mtx <- function(house, crisks) {
  n <- nrow(house)
  output <- matrix(1., nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      key <- paste0("class_", house$class[i], "_to_class_", house$class[j])
      if (i == j) output[i, j] <- 0
      else if (key %in% names(crisks)) output[i, j] <- crisks[[key]]
    }
  }
  output
}

# model over a house carrying the fixed risks and parameters
build_house_model <- function(house, par) {
  new_housepi_model(house, "final-size") |>
    add_group_risk(import = mk_house$groups) |>
    add_contact_risk(import = mk_house$contacts) |>
    # FIXME: this overwrites group risk parameters, opaque to user and risky
    add_parameter(import = par, append = FALSE)
}

# log-likelihood of a house through the chain binomial debugger
chain_binomial_log_lik <- function(house, par, crisks) {
  n_full <- nrow(house)
  n_infec <- sum(house$status > 0L)
  commu_infec_rate <- par$alpha * house$community
  p2p_infec_rate <- par$beta * (par$hhsize_ref / n_full) ^ par$gamma *
    house$infectivity %o% house$susceptibility * make_ctr_mtx(house, crisks)
  if (n_infec == 0) return(-sum(commu_infec_rate))

  infec <- seq_len(n_infec)
  escaped <- (n_infec + 1):n_full
  # escape probabilities restricted to the infected members
  B <- exp(-commu_infec_rate[infec])
  Q <- exp(-p2p_infec_rate[infec, infec, drop = FALSE])
  -sum(commu_infec_rate[escaped]) - sum(p2p_infec_rate[infec, escaped]) +
    log(eval_chain_binomial_debugger_(B, Q))
}

######################################
## FINAL SIZE HH MODEL : GUARDS     ##
######################################

test_that("Chain binomial guards against misusage", {
  # not discretizable
  dummy_household |>
    new_housepi_model("final-size", discr = TRUE) |>
    expect_error("No visible discrete time option")

  # no support for missing status
  data.frame(id = c("m1", "m2"), house_id = "h1", status = c(1, NA)) |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    run_inference() |>
    expect_error("missing status.*incompatible")

  # no infectivity profile
  flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_infectivity_profile(import = "template") |>
    expect_error("not supported.*dropped")
})

test_that("Debugger cannot be called on large numbers of infections", {
  n_infec <- 10L
  eval_chain_binomial_debugger_(
    runif(n_infec), matrix(runif(n_infec * n_infec), ncol = n_infec)
  ) |>
    expect_error("Maximum number")
})

##############################################
## FINAL SIZE HH MODEL : LOG-LIKELIHOOD     ##
##############################################

test_that("Single infected member log-likelihood is correct", {
  build_house_model(mk_house$house_1, mk_house$par_1) |>
    evaluate_log_lik() |>
    expect_equal(chain_binomial_log_lik(mk_house$house_1, mk_house$par_1,
                                        mk_house$crisks_1))
})

test_that("Three infected members likelihood is correct", {
  # no infected
  house <- dplyr::mutate(mk_house$house_3, status = 0)
  build_house_model(house, mk_house$par_3) |>
    evaluate_log_lik() |>
    expect_equal(chain_binomial_log_lik(house, mk_house$par_3,
                                        mk_house$crisks_3))

  # mixed infected and escaped
  build_house_model(mk_house$house_3, mk_house$par_3) |>
    evaluate_log_lik() |>
    expect_equal(chain_binomial_log_lik(mk_house$house_3, mk_house$par_3,
                                        mk_house$crisks_3))
})

test_that("A free baseline under an absolute risk is pinned on real data", {
  # households carrying a sibling, so every member level of the strata is seen
  hhmodel <- flu_2014_matsumoto |>
    dplyr::filter(house_id %in% utils::head(
      unique(house_id[member == "sibling"]), 500
    )) |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    add_group_risk("community", by = "member") |>
    add_contact_risk(
      strata = data.frame(
        from_member = c("student", "student", "sibling", "mother", "mother",
                        "father", "other", NA),
        to_member = c("student", "sibling", "sibling", "student", "sibling",
                      "mother", "other", NA),
        param_name = c("ctr_child_child", "ctr_child_child", "ctr_child_child",
                       "ctr_mother_child", "ctr_mother_child",
                       ".reference", "ctr_other_other",
                       "ctr_cross_generational")
      )
    )

  # the template declares a free alpha, which the absolute community risk
  # leaves no scale for: it is pinned, and the anomaly is reported
  hhmodel |>
    evaluate_log_lik() |>
    expect_warning("alpha.*turned.*reference")
})

test_that("Chain binomial returns infinite if extreme values are provided", {
  hhmodel <- data.frame(id = c(1, 2), house_id = "h1",
                        status = c(1, 1), age = c("child", "adult")) |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = 1e-20, beta = 1e-20, gamma = 1, hhsize_ref = 4)

  # homogeneous branch
  expect_false(is.finite(evaluate_log_lik(hhmodel)))

  # heterogeneous branch
  hhmodel |>
    add_group_risk("susceptibility", by = "age") |>
    evaluate_log_lik() |>
    is.finite() |>
    expect_false()
})

test_that("Chain binomial raises if the log-likelihood is not a number", {
  # a null person-to-person risk against an overflowing household-size scaling
  # leaves the escape rate undefined, and the log-likelihood with it
  data.frame(id = c(1, 2), house_id = "h1", status = c(1, 1)) |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = 0.1, beta = 0, gamma = 1e10, hhsize_ref = 4) |>
    evaluate_log_lik() |>
    expect_error("returned NaN")
})
