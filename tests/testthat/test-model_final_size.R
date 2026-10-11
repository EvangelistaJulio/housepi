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
    id = c("m1", "m2"), hh_id = "h1", status = c(1, 0),
    class = c("a", "b"),
    community = c(1, grisks$b$com),
    susceptibility = c(1, grisks$b$sus),
    infectivity = c(1, grisks$b$inf)
  )

  # parameters matching automatic names
  par <- c(list(alpha = 0.15, beta = 0.22, gamma = 1, hh_size_ref = 4),
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
    id = sprintf("m%s", seq(1, 4)), hh_id = "h1",
    status = c(rep(1, 3), 0),
    class = c("a", "b", rep("c", 2)),
    community = c(1, grisks$b$com, rep(grisks$c$com, 2)),
    susceptibility = c(1, grisks$b$sus, rep(grisks$c$sus, 2)),
    infectivity = c(1, grisks$b$inf, rep(grisks$c$inf, 2))
  )

  # parameters matching automatic names
  par <- c(list(alpha = 0.15, beta = 0.22, gamma = 1, hh_size_ref = 4),
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
    id = sprintf("m%s", seq(1, 4)), hh_id = "h1",
    status = c(rep(1, 3), 0),
    class = c("a", rep("c", 2), "b"),
    community = c(1, rep(grisks$c$com, 2), grisks$b$com),
    susceptibility = c(1, rep(grisks$c$sus, 2), grisks$b$sus),
    infectivity = c(1, rep(grisks$c$inf, 2), grisks$b$inf)
  )

  # parameters matching automatic names
  par <- c(list(alpha = 0.15, beta = 0.22, gamma = 1, hh_size_ref = 4),
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

# probability of every transmission chain among infected members, given their
# community escapes `B` and person-to-person escapes `Q` (infector by row)
# a chain orders the infected into generations: the first is infected by the
# community, a later member escapes the community and every generation before
# the previous one, and is infected by the previous one
enumerate_chains <- function(B, Q) {
  # ordered set partitions of `remaining`, one integer vector per generation
  ordered_partitions <- function(remaining) {
    if (!length(remaining)) return(list(list()))
    unlist(lapply(seq_along(remaining), function(n) {
      firsts <- utils::combn(seq_along(remaining), n,
                             function(i) remaining[i], simplify = FALSE)
      unlist(lapply(firsts, function(first) {
        lapply(ordered_partitions(setdiff(remaining, first)),
               \(rest) c(list(first), rest))
      }), recursive = FALSE)
    }), recursive = FALSE)
  }
  chain_proba <- function(gens) {
    log_p <- sum(log1p(-B[gens[[1L]]]))
    for (g in seq_along(gens)[-1L]) {
      earlier <- unlist(gens[seq_len(g - 2L)])
      for (i in gens[[g]]) {
        log_p <- log_p + log(B[i]) + sum(log(Q[earlier, i])) +
          log1p(-prod(Q[gens[[g - 1L]], i]))
      }
    }
    exp(log_p)
  }

  chains <- ordered_partitions(seq_along(B))
  list(chains = chains, proba = vapply(chains, chain_proba, numeric(1L)))
}

# log-likelihood of a house, summing its transmission chains
chain_binomial_log_lik <- function(house, par, ctr) {
  esc <- escape_probas(house, par, ctr)
  B <- esc$community_escape
  Q <- esc$p2p_escape
  infec <- house$status > 0L
  if (!any(infec)) return(sum(log(B)))
  sum(log(B[!infec])) + sum(log(Q[infec, !infec])) +
    log(sum(enumerate_chains(B[infec], Q[infec, infec, drop = FALSE])$proba))
}

# escape probabilities of a house, named after its members
escape_probas <- function(house, par, ctr) {
  n <- nrow(house)
  p2p_rate <- par$beta * (par$hh_size_ref / n) ^ par$gamma *
    house$infectivity %o% house$susceptibility * ctr
  list(
    community_escape = setNames(exp(-par$alpha * house$community), house$id),
    p2p_escape = matrix(exp(-p2p_rate), n, n,
                        dimnames = list(house$id, house$id))
  )
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
  data.frame(id = c("m1", "m2"), hh_id = "h1", status = c(1, NA)) |>
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

test_that("Likelihood is correct when using risk partitions", {
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

test_that("Likelihood sums to one over every final outcome of a household", {
  # one house per status vector of the members of `mk_part`: the infected take
  # the homogeneous (`c`, `c`), partitioned (`a`, `c`, `c`) and heterogeneous
  # (`a`, `b`) branches
  house <- mk_part$house
  outcomes <- as.matrix(expand.grid(rep(list(0:1), nrow(house))))
  houses <- do.call(rbind, lapply(seq_len(nrow(outcomes)), function(k) {
    transform(house, id = paste0(id, "_", k), hh_id = paste0("h", k),
              status = unname(outcomes[k, ]))
  }))
  hh_model <- houses |>
    build_house_model(mk_part$par, mk_part$groups, mk_part$contacts)
  evaluate_log_lik(hh_model)

  unique(houses$hh_id) |>
    vapply(\(hh_id) inspect_house_state(hh_model, hh_id)$transmission_log_lik,
           numeric(1)) |>
    exp() |>
    sum() |>
    expect_equal(1)
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
  strata_direct <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "m2",         "m3",       "ctr_m2_m3",
    NA,           NA,         ".reference"
  )
  ctr_direct <- mk_ctr_mtx(house, list())
  ctr_direct[2, 3] <- 2
  build_break(strata_direct, list(ctr_m2_m3 = 2)) |>
    evaluate_log_lik() |>
    expect_equal(chain_binomial_log_lik(house, mk_part$par, ctr_direct))

  # the pair itself still matches: only m2's contact toward the escaped
  # member (m4) differs from m3's, caught while scanning the other members
  strata_third <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "m2",         "m4",       "ctr_m2_m4",
    NA,           NA,         ".reference"
  )
  ctr_third <- mk_ctr_mtx(house, list())
  ctr_third[2, 4] <- 2
  build_break(strata_third, list(ctr_m2_m4 = 2)) |>
    evaluate_log_lik() |>
    expect_equal(chain_binomial_log_lik(house, mk_part$par, ctr_third))
})

test_that("A free baseline under an absolute risk is pinned on real data", {
  # households carrying a sibling, so every member level of the strata is seen
  hhdata <- flu_2014_matsumoto |>
    dplyr::filter(hh_id %in% utils::head(
      unique(hh_id[member == "sibling"]), 500
    ))

  # the template declares a free alpha, which the absolute community risk
  # leaves no scale for: it is pinned, and the anomaly is reported
  expect_warning(
    hh_model <- hhdata |>
      new_housepi_model("final-size") |>
      add_parameter(import = "template") |>
      add_group_risk("community", by = "member"),
    "alpha.*turned.*reference"
  )
  params <- get_parameters(hh_model)
  params$value[params$name == "alpha"] |> expect_equal(1)

  # the contact table names a reference, so `beta` keeps its scale
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "student",    "sibling",  "ctr_child_child",
    "mother",     "student",  "ctr_mother_child",
    "mother",     "sibling",  "ctr_mother_child",
    "father",     "mother",   ".reference",
    "other",      "other",    "ctr_other_other",
    NA,           NA,         "ctr_cross_generational"
  )
  hh_model |>
    add_contact_risk(strata = strata) |>
    get_parameters() |>
    subset(name == "beta", "const") |>
    unlist() |>
    expect_false()
})

test_that("Infection rates below the machine epsilon keep their probability", {
  hh_model <- data.frame(id = c(1, 2), hh_id = "h1",
                         status = c(1, 1), age = c("child", "adult")) |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = 1e-20, beta = 1e-20, gamma = 1, hh_size_ref = 4)

  # homogeneous branch: both infected from the community,
  # or one from the community and the other from it, at rate beta * 4 / 2
  p_commu <- -expm1(-1e-20)
  p_p2p <- -expm1(-2e-20)
  evaluate_log_lik(hh_model) |>
    expect_equal(log(p_commu^2 + 2 * p_commu * (1 - p_commu) * p_p2p))

  # single infected member, the other escaping both routes
  data.frame(id = c(1, 2), hh_id = "h1", status = c(1, 0)) |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = 1e-20, beta = 1e-20, gamma = 1, hh_size_ref = 4) |>
    evaluate_log_lik() |>
    expect_equal(log(p_commu) - 1e-20 - 2e-20)

  # heterogeneous branch: the risk names a reference, so `beta` stays extreme
  # UPDATE: the recursion takes `1 - sum` over probabilities and returns -Inf
  hh_model |>
    add_group_risk("susceptibility", by = "age", reference = ".auto") |>
    evaluate_log_lik() |>
    is.finite() |>
    expect_false()
})

test_that("Chain binomial raises if the log-likelihood is not a number", {
  # a null person-to-person risk against an overflowing household-size scaling
  # leaves the escape rate undefined, and the log-likelihood with it
  data.frame(id = c(1, 2), hh_id = "h1", status = c(1, 1)) |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = 0.1, beta = 0, gamma = 1e10, hh_size_ref = 4) |>
    evaluate_log_lik() |>
    expect_error("[Nn]ot a number")
})

################################################
###  FINAL SIZE HH MODEL : MODEL QUANTITIES  ###
################################################

test_that("Escape probabilities follow the transmission formulas per member", {
  # `h2` lists the members of `h1` in another order: one equivalence class
  house <- mk_heter$house
  house_2 <- transform(house[c(4, 2, 1, 3), ],
                       id = sprintf("m%s", 5:8), hh_id = "h2")
  hh_model <- rbind(house, house_2) |>
    build_house_model(mk_heter$par, mk_heter$groups, mk_heter$contacts)
  evaluate_log_lik(hh_model)
  summary <- summarise_house_cases(hh_model)
  summary$count[summary$metric == "classes"] |> expect_equal(1)

  for (hhdata in list(house, house_2)) {
    ctr <- mk_ctr_mtx(hhdata, mk_heter$crisks)
    inspect_house_state(hh_model, hhdata$hh_id[1])$model_quantities |>
      expect_equal(escape_probas(hhdata, mk_heter$par, ctr))
  }
})

test_that("Escape probabilities follow the last accepted parameters of a run", {
  # the last move of a run may be rejected
  hh_model <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")
  for (i in 1:3) {
    run_inference(hh_model, steps = 10L, thin = 1L, adapt = 0L, chains = 1L,
                  progress = FALSE, warn = FALSE)

    # three members, default group and contact risks
    snapshot <- inspect_house_state(hh_model, "119")
    n <- nrow(snapshot$group_risks)
    snapshot$model_quantities |>
      expect_equal(escape_probas(snapshot$group_risks,
                                 as.list(snapshot$parameters), 1 - diag(n)))
  }
})
