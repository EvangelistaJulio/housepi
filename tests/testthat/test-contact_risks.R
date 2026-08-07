library(testthat)

# a dummy dataset for the plain logic-table helpers
dummy_df <- data.frame(house_id = 1,
                       age = c(rep("child", 2L), rep("adult", 2L)),
                       sex = c("M", "F", "M", "M"))

# model over the dummy features
mk_risk <- local({
  hhmodel <- data.frame(
    house_id = 1L, status = 0L,
    age = c("child", "child", "adult", "adult"),
    sex = c("M", "F", "M", "M")
  ) |>
    new_housepi_model("final_size")
  list(model = hhmodel)
})

#################################
## CONTACT RISKS : PAIR COUNTS ##
#################################

test_that("Contact risks pair counts are computed correctly", {
  df_hh <- data.frame(house_id = c(rep("hh1", 2L), rep("hh2", 2L)),
                      age = rep(c("child", "adult"), 2L),
                      stringsAsFactors = FALSE)
  tbl <- contact_risk_logic_table(df_hh, "age")

  cc <- subset(tbl, from_age == "child" & to_age == "child")$n
  ca <- subset(tbl, from_age == "child" & to_age == "adult")$n
  ac <- subset(tbl, from_age == "adult" & to_age == "child")$n
  aa <- subset(tbl, from_age == "adult" & to_age == "adult")$n

  # total: child->adult = 2, adult->child = 2
  # no child->child = 0, adult->adult = 0
  expect_equal(ca, 2L)
  expect_equal(ac, 2L)
  expect_equal(cc, integer(0))
  expect_equal(aa, integer(0))
})

######################################
## CONTACT RISKS : AUTO DEFINITIONS ##
######################################

test_that("Auto mode imposes symmetric contacts by default", {
  ltab <- mk_risk$model |>
    add_contact_risk(by = "age",
                     reference = c(from_age = "adult", to_age = "adult"),
                     prior = "lognormal", sdlog = 1.) |>
    get_contact_risks()

  # two-level feature: child, adult -> 4 pairs
  expect_equal(nrow(ltab), 4L)

  # only 1 reference
  expect_equal(sum(ltab$param_name == ".reference"), 1)
  # only 2 unique non-reference params (cross + child->child)
  expect_length(setdiff(unique(ltab$param_name), ".reference"), 2)

  # reference row is adult->adult
  ref_row <- ltab[ltab$param_name == ".reference", ]
  expect_equal(ref_row$from_age, "adult")
  expect_equal(ref_row$to_age, "adult")

  # symmetric: adult->child and child->adult share the same parameter
  expect_equal(
    subset(ltab, from_age == "adult" & to_age == "child")$param_name,
    subset(ltab, from_age == "child" & to_age == "adult")$param_name
  )
})

test_that("Auto mode without reference picks most abundant pair", {
  ltab <- mk_risk$model |>
    add_contact_risk(by = "age") |>
    get_contact_risks()

  # automatic .reference
  expect_true(".reference" %in% ltab$param_name)
})

test_that("Auto mode symmetric = FALSE gives distinct parameters", {
  ltab <- mk_risk$model |>
    add_contact_risk(by = "age", symmetric = FALSE,
                     reference = c(from_age = "adult", to_age = "adult"),
                     prior = "lognormal", sdlog = 1.) |>
    get_contact_risks()

  expect_equal(nrow(ltab), 4L)
  expect_length(setdiff(unique(ltab$param_name), ".reference"), 3L)

  # asymmetric: adult->child and child->adult are distinct params
  ac_row <- subset(ltab, from_age == "adult" & to_age == "child")
  ca_row <- subset(ltab, from_age == "child" & to_age == "adult")
  expect_true(ac_row$param_name != ca_row$param_name)
})

test_that("Non-joint produces independent per-feature tables", {
  # 2 features, non-joint -> 2 rule groups
  tidy <- mk_risk$model |>
    add_contact_risk(by = c("age", "sex"), joint = FALSE,
                     reference = c(age = "adult", sex = "M"),
                     prior = "lognormal", sdlog = 1.) |>
    get_contact_risks()

  expect_equal(length(unique(tidy$rule_group)), 2L)

  # all age pairs are observed: 2 x 2
  expect_equal(sum(tidy$rule_group == 1), 4L)
  # 1 sex pair is missing: 2 x 2 - 1
  expect_equal(sum(tidy$rule_group == 2), 3L)

  # each sub-table has its own .reference
  for (grp in unique(tidy$rule_group)) {
    sub <- subset(tidy, rule_group == grp)
    expect_equal(sum(sub$param_name == ".reference"), 1L)
  }

  # symmetric default: 2 + 2 - 1 = 3 non-reference params
  expect_length(setdiff(unique(tidy$param_name), ".reference"), 3L)

  # tidy table merges sub-tables with NA for missing columns
  expect_equal(nrow(tidy), 7L)
})

test_that("Auto mode with two features behaves correctly", {
  ltab <- mk_risk$model |>
    add_contact_risk(by = c("age", "sex"), joint = TRUE,
                     reference = c(from_age = "adult", from_sex = "M",
                                   to_age = "adult", to_sex = "M"),
                     prior = "lognormal", sdlog = 1.) |>
    get_contact_risks()

  # NOTE: 3 observed types: <child,M>, <child,F>, <adult,M>
  # 7 observed pairs: 9 pairs - <child,M><child,F> - (<->)
  expect_equal(nrow(ltab), 7L)
  expect_equal(length(unique(ltab$rule_group)), 1L)
  expect_equal(sum(ltab$param_name == ".reference"), 1)

  # NOTE: symmetric pairs by default:
  # 3 heterogeneous: <child,M><child,F> <adult,M><child,F> <adult,M><child,M>
  # 1 homogeneous: <adult,M><adult,M> (reference)
  expect_length(setdiff(unique(ltab$param_name), ".reference"), 3L)
})

########################################
## CONTACT RISKS : CUSTOM DEFINITIONS ##
########################################

test_that("Contact risk custom mode defualt behaves correctly", {
  # logic table
  param_name <- c("a2c", "c2a", "c2c")
  strata <- data.frame(from_age = c(rep("adult", 2L), rep("child", 2L)),
                       to_age = rep(c("adult", "child"), 2L),
                       param_name = c(".reference", param_name))
  # custom fine grained parameter override
  overrides <- list(
    list(value = 1, prior = "uniform", min = 1., max = 5.,
         proposal = "rw", proposal_scale = 0.5),
    2.34
  )
  fancy_param_1 <- param_name[1] # to override as variable
  fancy_param_2 <- param_name[2] # to override as constant
  default_param <- param_name[3] # default definitions
  names(overrides) <- c(fancy_param_1, fancy_param_2)

  hhmodel <- mk_risk$model |>
    add_contact_risk(strata = strata, prior = "lognormal", sdlog = 1.8,
                     params = overrides)

  # the custom logic table drives the tidy table
  tidy <- get_contact_risks(hhmodel)
  expect_true(is.data.frame(tidy))
  expect_equal(nrow(tidy), 4L)
  expect_setequal(tidy$param_name, strata$param_name)

  # parameters generated for non-reference
  params <- get_parameters(hhmodel)
  expect_setequal(params$name, param_name)

  # default parameter definitions
  default <- subset(params, name == default_param)
  expect_equal(default$prior, "lognormal")
  expect_equal(default$args[[1]], c(sdlog = 1.8))
  # overrides as variable
  fancy_1 <- subset(params, name == fancy_param_1)
  expect_equal(fancy_1$prior, "uniform")
  expect_equal(fancy_1$min, 1)
  # overrides as constant
  fancy_2 <- subset(params, name == fancy_param_2)
  expect_equal(fancy_2$value, 2.34)
  expect_true(fancy_2$const)
})

test_that("Custom mode with no reference row raises error", {
  strata <- data.frame(
    from_age = c("adult", "adult", "child", "child"),
    to_age   = c("adult", "child", "adult", "child"),
    param_name = c("a_to_a", "a_to_c", "c_to_a", "c_to_c")
  )

  mk_risk$model |>
    add_contact_risk(strata = strata) |>
    expect_error("[Cc]ustom logic table must contain a ..reference. row")
})

####################################################
## CONTACT RISKS : CUSTOM (SYMMETRIC) DEFINITIONS ##
####################################################

test_that("Custom mode allows symmetric definition automatic fixes", {
  param_name <- c("a_to_c", "c_to_c")
  strata <- data.frame(from_age = c("adult", "adult", "child"),
                       to_age = c("adult", "child", "child"),
                       param_name = c(".reference", param_name))
  overrides <- list(a_to_c = 2)

  hhmodel <- mk_risk$model |>
    add_contact_risk(strata = strata, params = overrides)

  # symmetric fill adds child->adult with same param_name as adult->child
  filled <- get_contact_risks(hhmodel)
  expect_equal(nrow(filled), 4)
  ca_row <- filled[filled$from_age == "child" & filled$to_age == "adult", ]
  expect_equal(ca_row$param_name, "a_to_c")

  # both directions share the override value
  params <- get_parameters(hhmodel)
  expect_setequal(params$name, param_name)
  a_to_c <- subset(params, name == "a_to_c")
  expect_equal(a_to_c$value, 2)
  expect_true(a_to_c$const)
})

test_that("Custom mode with symmetric = FALSE does not fill gaps", {
  param_name <- c("a_to_c", "c_to_c")
  strata <- data.frame(from_age = c("adult", "adult", "child"),
                       to_age = c("adult", "child", "child"),
                       param_name = c(".reference", param_name))

  tidy <- mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    get_contact_risks()

  expect_equal(nrow(tidy), 3L)
})

#########################################
## CONTACT RISKS : RULE GROUPING LOGIC ##
#########################################

test_that("compute_rule_groups: ambiguous subset raises error", {
  data.frame(
    from_age = c("child", "child", "child"),
    from_sex = c("M", NA, NA),
    to_age = c(NA, "child", NA),
    param_name = c("metal_gear_1", "metal_gear_2", "big_boss")
  ) |>
    compute_rule_groups() |>
    expect_error("[Aa]mbiguous")
})

test_that("Custom mode absorbs all-NA default into single group", {
  # Example 1: 4 rules in one group, last one (all-NA) is the catch-all
  strata <- data.frame(
    from_member = c("student", "mother", "father", NA),
    to_member = c("student", "student", "mother", NA),
    param_name = c("ctr_child_child", "ctr_mother_child",
                   ".reference", "ctr_cross_gen")
  )
  tidy <- data.frame(house_id = 1, status = 0,
                     member = c("student", "mother", "father", "student")) |>
    new_housepi_model("final_size") |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    get_contact_risks()
  expect_equal(length(unique(tidy$rule_group)), 1L)
})

test_that("Custom mode with two maximal patterns yields two groups", {
  # Example 2: disjoint age vs sex patterns -> 2 distinct rule_groups
  strata <- data.frame(
    from_age = c("child", NA, NA),
    from_sex = c(NA, "M", "F"),
    to_age = c("child", NA, NA),
    to_sex = c(NA, "M", "F"),
    param_name = c("age_child_to_age_child", ".reference", "sex_F_to_sex_F")
  )
  tidy <- mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    get_contact_risks()
  expect_equal(length(unique(tidy$rule_group)), 2L)
})

test_that("Auto mode joint = FALSE yields one rule_group per sub-table", {
  tidy <- mk_risk$model |>
    add_contact_risk(by = c("age", "sex"), joint = FALSE,
                     reference = c(age = "adult", sex = "M")) |>
    get_contact_risks()
  expect_equal(length(unique(tidy$rule_group)), 2L)
})

test_that("Auto mode joint = TRUE yields a single rule_group", {
  tidy <- mk_risk$model |>
    add_contact_risk(by = c("age", "sex"), joint = TRUE,
                     reference = c(from_age = "adult", from_sex = "M",
                                   to_age = "adult", to_sex = "M")) |>
    get_contact_risks()
  expect_equal(length(unique(tidy$rule_group)), 1L)
})

test_that("Non-joint sub-tables have correct columns", {
  tidy <- mk_risk$model |>
    add_contact_risk(by = c("age", "sex"), joint = FALSE,
                     prior = "lognormal", sdlog = 1) |>
    get_contact_risks()

  # first sub-table carries age, and leaves sex unset
  t1 <- subset(tidy, rule_group == 1)
  expect_true(all(!is.na(t1$from_age) & !is.na(t1$to_age)))
  expect_true(all(is.na(t1$from_sex) & is.na(t1$to_sex)))

  # second sub-table carries sex, and leaves age unset
  t2 <- subset(tidy, rule_group == 2)
  expect_true(all(!is.na(t2$from_sex) & !is.na(t2$to_sex)))
  expect_true(all(is.na(t2$from_age) & is.na(t2$to_age)))
})

#####################################################
## CONTACT RISKS : C++ RESOLUTION (best_match.h)   ##
#####################################################

# These tests drive update_best_match / resolve_best_match in src/parameter.cpp
# through the full model build, so the C++ matrix can be read back and
# inspected via export_house_snapshot.

# Cell-by-cell readout helper.
ctr_cell_ <- function(hhsnapshot, from, to, key = "member") {
  ids <- hhsnapshot$data_wide[[key]]
  hhsnapshot$contact_risks[match(from, ids), match(to, ids)]
}

test_that("Most-specific rule wins over reference (single feature)", {
  # one specific pair gets a non-reference param; everything else falls
  # through to .reference (specificity 0, hits update_best_match's empty
  # branch then never gets beaten).
  strata <- data.frame(
    from_member = c("father", NA),
    to_member = c("mother", NA),
    param_name = c("ctr_fm", ".reference"),
    stringsAsFactors = FALSE
  )

  # house with one of each member type
  hhmodel <- data.frame(
    id = c("m1", "m2", "m3", "m4"),
    house_id = "h1",
    status = c(0L, 0L, 0L, 0L),
    member = c("father", "mother", "student", "sibling"),
    stringsAsFactors = FALSE
  ) |>
    new_housepi_model("final_size") |>
    add_parameter(import = "template") |>
    add_contact_risk(strata = strata, symmetric = FALSE,
                     params = list(ctr_fm = 2))
  snap <- hhmodel$export_house_snapshot("h1")

  # the explicit pair wins
  expect_equal(ctr_cell_(snap, "father", "mother"), 2)
  # symmetric = FALSE: reverse falls back to reference
  expect_equal(ctr_cell_(snap, "mother", "father"), 1)
  # unrelated pair falls back to reference (spec=0)
  expect_equal(ctr_cell_(snap, "student", "sibling"), 1)
  # diagonal stays 0 (no self-contact)
  unname(diag(snap$contact_risks)) |>
    expect_equal(rep(0., nrow(snap$contact_risks)))
})

test_that("Distinct rule_groups multiply at the cell", {
  # two-feature household so 'compute_rule_groups' yields TWO maximal
  # patterns ({from_age,to_age}, {from_class,to_class}) and the per-group
  # winners get init_product'd into one Attribute.
  hh <- data.frame(
    id = c("m1", "m2", "m3", "m4"),
    house_id = "h1",
    status = c(0L, 0L, 0L, 0L),
    age = c("adult", "adult", "child", "adult"),
    class = c("a", "a", "b", "b"),
    stringsAsFactors = FALSE
  )

  # disjoint patterns -> two rule_groups; prime values catch any
  # accidental sum / overwrite. .reference is dropped (".none") so the
  # two maximal patterns don't conflict in compute_rule_groups; the C++
  # side still fills unmatched cells with VALUE_REFERENCE (= 1).
  strata <- data.frame(
    from_age = c("adult", NA),
    to_age = c("adult", NA),
    from_class = c(NA, "a"),
    to_class = c(NA, "a"),
    param_name = c("ctr_aa", "ctr_class_aa"),
    stringsAsFactors = FALSE
  )

  hhmodel <- new_housepi_model(hh, "final_size")
  hhmodel <- add_parameter(hhmodel, import = "template")
  hhmodel <- add_contact_risk(hhmodel,
    strata = strata, symmetric = FALSE, reference = ".none",
    params = list(ctr_aa = 3, ctr_class_aa = 5)
  )
  snap <- hhmodel$export_house_snapshot("h1")

  # m1 (adult,a) → m2 (adult,a): BOTH rules match → 3 * 5 = 15
  # (init_product over two Attribute* — parameter.cpp's multi-group loop)
  expect_equal(ctr_cell_(snap, "m1", "m2", key = "id"), 15)
  # m1 (adult,a) → m4 (adult,b): only age rule matches → 3
  expect_equal(ctr_cell_(snap, "m1", "m4", key = "id"), 3)
  # m1 (adult,a) → m3 (child,b): neither rule matches → reference (1)
  expect_equal(ctr_cell_(snap, "m1", "m3", key = "id"), 1)
  # diagonal is 0
  expect_equal(unname(diag(snap$contact_risks)), rep(0., nrow(hh)))
})

test_that("Contact risk definition work as expected", {
  hhid <- "hh2"
  hhmodel <- data.frame(
    house_id = c(rep("hh1", 3), rep(hhid, 6)),
    status = c(1, 0, 0, 1, 0, 0, 0, 0, 0),
    member = c("student", "other", "other",
               "student", "sibling", "mother", "father", "other", "student")
  ) |>
    new_housepi_model("final_size") |>
    add_parameter(import = "template")

  # no contact risk defs -> empty matrix
  hhsnapshot <- hhmodel$export_house_snapshot(hhid)
  cr_mtx <- hhsnapshot$contact_risks
  expect_true(is.matrix(cr_mtx))
  expect_equal(nrow(cr_mtx), 0)
  expect_equal(ncol(cr_mtx), 0)

  # mutually exclusive strata
  strata <- data.frame(
    from_member = c(
      "student", "student", "sibling", "mother",
      "mother", "father", "other", NA
    ),
    to_member = c(
      "student", "sibling", "sibling", "student",
      "sibling", "mother", "other", NA
    ),
    param_name = c(
      "child_child", "child_child", "child_child", "mother_child",
      "mother_child", ".reference", "other_other", "cross_gen"
    )
  )
  # prime numbers as param values to ensure non duplicate pair productes
  params <- c(2, 3, 5, 7, 11, 13, 17, 19)
  names(params) <- strata$param_name
  hhmodel <- hhmodel |>
    add_contact_risk(strata = strata, params = params)

  # matrix is n x n
  hhsnapshot <- hhmodel$export_house_snapshot(hhid)
  cr_mtx <- hhsnapshot$contact_risks
  n <- nrow(hhsnapshot$data_wide)

  expect_true(is.matrix(cr_mtx))
  expect_equal(nrow(cr_mtx), n)
  expect_equal(ncol(cr_mtx), n)
  # diagonal is 0
  expect_equal(unname(diag(cr_mtx)), rep(0., n))

  # every off-diagonal cell resolves to the reference (1)
  # or to exactly one param value: no multiplications
  expect_true(all(cr_mtx %in% c(0, 1, params)))

  # the template parameters added previously survive alongside the risk ones
  template_pars <- c("alpha", "beta", "gamma")
  expect_true(all(template_pars %in% hhmodel$parameters$name))
})

test_that("Contact risks expand on their own covariates only", {
  df <- sarscov2_2021_icofs |>
    dplyr::select(-start_date) |>
    dplyr::filter(house_id %in% utils::head(unique(house_id), 50))

  # a data column the risk does not read leaves the expansion untouched
  risks <- df |>
    new_housepi_model("case_ascertained") |>
    add_contact_risk(by = "age") |>
    get_contact_risks()
  expect_setequal(unique(risks$from_age), c("adult", "child"))
  expect_true("age_child_to_age_child" %in% risks$param_name)

  # the covariate itself must be a genuine data column
  df |>
    dplyr::select(-age) |>
    new_housepi_model("case_ascertained") |>
    add_contact_risk(by = "age") |>
    expect_error("[Uu]nknown feature.*age")

  # without data the read surface warns and keeps the table shape
  expect_warning(
    risks <- new_housepi_model(model = "final_size") |>
      add_contact_risk(by = "member") |>
      get_contact_risks(),
    regexp = "without data"
  )
  expect_equal(nrow(risks), 0L)
  expect_named(risks, c("param_name", "n", "rule_group"))
})

test_that("A contact risk keyed on an observation is rejected", {
  # the risk feeds the transmission model the observation is generated from
  df <- data.frame(house_id = c(1, 1), status = c(1, 0),
                   severity = c(2, 3), test_date = c(2, 6))

  df |>
    new_housepi_model("final_size") |>
    add_contact_risk(by = "severity") |>
    add_observation("severity", model = "uniform", min = 0, max = 5) |>
    expect_error("keyed on observation.*severity")

  # a deterministic household observation stands in for the column all the same
  df |>
    new_housepi_model("final_size") |>
    add_observation(severity = min(test_date@D), level = "house") |>
    add_contact_risk(by = "severity") |>
    expect_error("keyed on observation.*severity")
})

test_that("Contact risks raise an error for incompatible models", {
  data.frame(
    house_id = 1, id = seq(4), status = seq(4),
    cohort = c("a", "a", "b", "b")
  ) |>
    new_housepi_model("test") |>
    add_contact_risk(by = "cohort") |>
    expect_error("[Nn]ot supported")
})
