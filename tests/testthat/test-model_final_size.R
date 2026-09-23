library(testthat)

# infected member in class `a`, escaped member in class `b`: the homogeneous
# branch (single risk partition). Class `a` is the reference class.
mk_homo <- local({
  groups <- list(community = list("class", reference = "a"),
                 susceptibility = list("class", reference = "a"),
                 infectivity = list("class", reference = "a"))
  contacts <- list(by = "class", symmetric = FALSE,
                   reference = c(from_class = "a", to_class = "b"))
  grisks <- list(a = list(com = 1., sus = 1., inf = 1.),
                 b = list(com = 2., sus = 3., inf = 0.2))
  crisks <- list(class_b_to_class_a = 2.3)

  house <- data.frame(
    id = c("m1", "m2"), hhid = "h1", status = c(1, 0),
    class = c("a", "b"),
    community = c(1, grisks$b$com),
    susceptibility = c(1, grisks$b$sus),
    infectivity = c(1, grisks$b$inf)
  )

  # parameters matching automatic names
  par <- c(list(alpha = 0.15, beta = 0.22, gamma = 1, hhsize_ref = 4),
           crisks,
           list(com_class_b = grisks$b$com,
                sus_class_b = grisks$b$sus,
                inf_class_b = grisks$b$inf))

  list(groups = groups, contacts = contacts, house = house,
       crisks = crisks, par = par)
})

# three risk classes, members sorted infected then escaped: `a`, `b`, `c` each
# infected in a distinct class, taking the heterogeneous branch (fully
# split). Class `a` is the reference class.
mk_heter <- local({
  groups <- list(community = list("class", reference = "a"),
                 susceptibility = list("class", reference = "a"),
                 infectivity = list("class", reference = "a"))
  contacts <- list(by = "class", symmetric = FALSE,
                   reference = c(from_class = "a", to_class = "b"))
  grisks <- list(a = list(com = 1., sus = 1., inf = 1.),
                 b = list(com = 2., sus = 3., inf = 0.2),
                 c = list(com = 0.5, sus = 0.8, inf = 1.5))
  crisks <- list(class_b_to_class_a = 2.3,
                 class_a_to_class_c = 0.7,
                 class_c_to_class_a = 0.6,
                 class_b_to_class_c = 1.1,
                 class_c_to_class_b = 0.9,
                 class_c_to_class_c = 1.2)

  house <- data.frame(
    id = sprintf("m%s", seq(1, 4)), hhid = "h1",
    status = c(rep(1, 3), 0),
    class = c("a", "b", rep("c", 2)),
    community = c(1, grisks$b$com, rep(grisks$c$com, 2)),
    susceptibility = c(1, grisks$b$sus, rep(grisks$c$sus, 2)),
    infectivity = c(1, grisks$b$inf, rep(grisks$c$inf, 2))
  )

  # parameters matching automatic names
  par <- c(list(alpha = 0.15, beta = 0.22, gamma = 1, hhsize_ref = 4),
           crisks,
           list(com_class_b = grisks$b$com,
                sus_class_b = grisks$b$sus,
                inf_class_b = grisks$b$inf,
                com_class_c = grisks$c$com,
                sus_class_c = grisks$c$sus,
                inf_class_c = grisks$c$inf))

  list(groups = groups, contacts = contacts, house = house,
       crisks = crisks, par = par)
})

# same three risk classes, but two of the three infected members share class
# `c`: neither homogeneous nor fully split, taking the partitioned branch.
mk_part <- local({
  groups <- list(community = list("class", reference = "a"),
                 susceptibility = list("class", reference = "a"),
                 infectivity = list("class", reference = "a"))
  contacts <- list(by = "class", symmetric = FALSE,
                   reference = c(from_class = "a", to_class = "b"))
  grisks <- list(a = list(com = 1., sus = 1., inf = 1.),
                 b = list(com = 2., sus = 3., inf = 0.2),
                 c = list(com = 0.5, sus = 0.8, inf = 1.5))
  crisks <- list(class_b_to_class_a = 2.3,
                 class_a_to_class_c = 0.7,
                 class_c_to_class_a = 0.6,
                 class_b_to_class_c = 1.1,
                 class_c_to_class_b = 0.9,
                 class_c_to_class_c = 1.2)

  house <- data.frame(
    id = sprintf("m%s", seq(1, 4)), hhid = "h1",
    status = c(rep(1, 3), 0),
    class = c("a", rep("c", 2), "b"),
    community = c(1, rep(grisks$c$com, 2), grisks$b$com),
    susceptibility = c(1, rep(grisks$c$sus, 2), grisks$b$sus),
    infectivity = c(1, rep(grisks$c$inf, 2), grisks$b$inf)
  )

  # parameters matching automatic names
  par <- c(list(alpha = 0.15, beta = 0.22, gamma = 1, hhsize_ref = 4),
           crisks,
           list(com_class_b = grisks$b$com,
                sus_class_b = grisks$b$sus,
                inf_class_b = grisks$b$inf,
                com_class_c = grisks$c$com,
                sus_class_c = grisks$c$sus,
                inf_class_c = grisks$c$inf))

  list(groups = groups, contacts = contacts, house = house,
       crisks = crisks, par = par)
})

# contact risk matrix of a house, the reference reading 1
mk_ctr_mtx <- function(house, crisks) {
  n <- nrow(house)
  output <- matrix(1., nrow = n, ncol = n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    key <- paste0("class_", house$class[i], "_to_class_", house$class[j])
    if (i == j) output[i, j] <- 0
    else if (key %in% names(crisks)) output[i, j] <- crisks[[key]]
  }
  output
}

# model over a house carrying the fixed risks and parameters
build_house_model <- function(house, par, groups, contacts) {
  new_housepi_model(house, "final-size") |>
    add_group_risk(import = groups) |>
    add_contact_risk(import = contacts) |>
    # FIXME: this overwrites group risk parameters, opaque to user and risky
    add_parameter(import = par, append = FALSE)
}

# log-likelihood of a house through the chain binomial debugger
chain_binomial_log_lik <- function(house, par, ctr) {
  n_full <- nrow(house)
  n_infec <- sum(house$status > 0L)
  commu_infec_rate <- par$alpha * house$community
  p2p_infec_rate <- par$beta * (par$hhsize_ref / n_full) ^ par$gamma *
    house$infectivity %o% house$susceptibility * ctr
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
###  FINAL SIZE HH MODEL : GUARDS  ###
######################################

test_that("Chain binomial guards against misusage", {
  # not discretizable
  dummy_household |>
    new_housepi_model("final-size", discr = TRUE) |>
    expect_error("No visible discrete time option")

  # no support for missing status
  data.frame(id = c("m1", "m2"), hhid = "h1", status = c(1, NA)) |>
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
###  FINAL SIZE HH MODEL : LOG-LIKELIHOOD  ###
##############################################

test_that("Single infected member log-likelihood is correct", {
  build_house_model(mk_homo$house, mk_homo$par,
                    mk_homo$groups, mk_homo$contacts) |>
    evaluate_log_lik() |>
    expect_equal(chain_binomial_log_lik(
      mk_homo$house, mk_homo$par,
      mk_ctr_mtx(mk_homo$house, mk_homo$crisks)
    ))
})

test_that("Three infected members likelihood is correct", {
  # no infected
  house <- dplyr::mutate(mk_heter$house, status = 0)
  build_house_model(house, mk_heter$par,
                    mk_heter$groups, mk_heter$contacts) |>
    evaluate_log_lik() |>
    expect_equal(chain_binomial_log_lik(
      house, mk_heter$par, mk_ctr_mtx(house, mk_heter$crisks)
    ))

  # mixed infected and escaped
  build_house_model(mk_heter$house, mk_heter$par,
                    mk_heter$groups, mk_heter$contacts) |>
    evaluate_log_lik() |>
    expect_equal(chain_binomial_log_lik(
      mk_heter$house, mk_heter$par,
      mk_ctr_mtx(mk_heter$house, mk_heter$crisks)
    ))
})

test_that("Likelihood is correct when infected members fall into risk partitions", {
  # two of the three infected members are risk-equivalent, so the house
  # is neither homogeneous nor fully split, taking the partitioned branch
  build_house_model(mk_part$house, mk_part$par,
                    mk_part$groups, mk_part$contacts) |>
    evaluate_log_lik() |>
    expect_equal(chain_binomial_log_lik(
      mk_part$house, mk_part$par,
      mk_ctr_mtx(mk_part$house, mk_part$crisks)
    ))
})

test_that("A contact-risk mismatch keeps risk-equivalent members apart", {
  # m2 and m3 share class `c`, so group risk alone would merge them: only the
  # contact-risk check can still tell them apart. `member` is a fresh, unique
  # key so the override lands on one member without touching its class.
  house <- transform(mk_part$house, member = id)
  par_base <- mk_part$par[!names(mk_part$par) %in% names(mk_part$crisks)]

  build_break <- function(strata, override) {
    house |>
      new_housepi_model("final-size") |>
      add_group_risk(import = mk_part$groups) |>
      add_contact_risk(strata = strata, symmetric = FALSE) |>
      add_parameter(import = c(par_base, override), append = FALSE)
  }

  # the pair's own contact breaks the match directly: m2 -> m3 stands out
  # against the class-`c` reference, so same_contacts fails at the first check
  strata_direct <- data.frame(from_member = c("m2", NA),
                              to_member = c("m3", NA),
                              param_name = c("ctr_m2_m3", ".reference"))
  ctr_direct <- mk_ctr_mtx(house, list())
  ctr_direct[2, 3] <- 2
  build_break(strata_direct, list(ctr_m2_m3 = 2)) |>
    evaluate_log_lik() |>
    expect_equal(chain_binomial_log_lik(house, mk_part$par, ctr_direct))

  # the pair itself still matches: only m2's contact toward the escaped
  # member (m4) differs from m3's, caught while scanning the other members
  strata_third <- data.frame(from_member = c("m2", NA), to_member = c("m4", NA),
                             param_name = c("ctr_m2_m4", ".reference"))
  ctr_third <- mk_ctr_mtx(house, list())
  ctr_third[2, 4] <- 2
  build_break(strata_third, list(ctr_m2_m4 = 2)) |>
    evaluate_log_lik() |>
    expect_equal(chain_binomial_log_lik(house, mk_part$par, ctr_third))
})

test_that("A free baseline under an absolute risk is pinned on real data", {
  # households carrying a sibling, so every member level of the strata is seen
  hhdata <- flu_2014_matsumoto |>
    dplyr::filter(hhid %in% utils::head(
      unique(hhid[member == "sibling"]), 500
    ))

  # the template declares a free alpha, which the absolute community risk
  # leaves no scale for: it is pinned, and the anomaly is reported
  expect_warning(
    hhmodel <- hhdata |>
      new_housepi_model("final-size") |>
      add_parameter(import = "template") |>
      add_group_risk("community", by = "member"),
    "alpha.*turned.*reference"
  )
  params <- get_parameters(hhmodel)
  expect_equal(params$value[params$name == "alpha"], 1)

  # the contact table names a reference, so `beta` keeps its scale
  hhmodel |>
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
    ) |>
    get_parameters() |>
    subset(name == "beta", "const") |>
    unlist() |>
    expect_false()
})

test_that("Chain binomial returns infinite if extreme values are provided", {
  hhmodel <- data.frame(id = c(1, 2), hhid = "h1",
                        status = c(1, 1), age = c("child", "adult")) |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = 1e-20, beta = 1e-20, gamma = 1, hhsize_ref = 4)

  # homogeneous branch
  expect_false(is.finite(evaluate_log_lik(hhmodel)))

  # heterogeneous branch: the risk names a reference, so `beta` stays extreme
  hhmodel |>
    add_group_risk("susceptibility", by = "age", reference = ".auto") |>
    evaluate_log_lik() |>
    is.finite() |>
    expect_false()
})

test_that("Chain binomial raises if the log-likelihood is not a number", {
  # a null person-to-person risk against an overflowing household-size scaling
  # leaves the escape rate undefined, and the log-likelihood with it
  data.frame(id = c(1, 2), hhid = "h1", status = c(1, 1)) |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = 0.1, beta = 0, gamma = 1e10, hhsize_ref = 4) |>
    evaluate_log_lik() |>
    expect_error("returned NaN")
})
