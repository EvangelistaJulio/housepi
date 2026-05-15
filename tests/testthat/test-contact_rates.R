library(testthat)

aux_add <- function(self, by = NA_character_,
                    reference = NULL, 
                    joint = TRUE, symmetric = TRUE,
                    strata = NULL, params = NULL,
                    prior = NA_character_, ...,
                    min = NA_real_, max = NA_real_,
                    proposal = NA_character_,
                    proposal_scale = NA_real_,
                    prefix = "ctr",
                    .warn_over = TRUE) {
  self$add_quo(by = by,
               reference = reference,
               joint = joint, symmetric = symmetric,
               strata = strata,
               params = params,
               prior = prior, rlang::enquos(...),
               min = rlang::enquo(min),
               max = rlang::enquo(max),
               proposal = proposal,
               proposal_scale = proposal_scale,
               prefix = prefix,
               .warn_over = .warn_over)
  invisible(self)
}

# a dummy dataset and relative information for factorization
dummy_df <- data.frame(house_id = 1,
                       age = c(rep("child", 2), rep("adult", 2)),
                       sex = c("M", "F", "M", "M"))
# corresponding factors levels
dummy_factor_info <- list(age = list(levels = c("child", "adult")),
                          sex = list(levels = c("M", "F")))

#################################
## CONTACT RATES : PAIR COUNTS ##
#################################

test_that("Contact rates pair counts are computed correctly", {
  df_hh <- data.frame(house_id = c(rep("hh1", 2), rep("hh2", 2)),
                      age = rep(c("child", "adult"), 2),
                      stringsAsFactors = FALSE)
  tbl <- contact_rate_logic_table(df_hh, "age")

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
## CONTACT RATES : AUTO DEFINITIONS ##
######################################

test_that("Auto mode imposes symmetric contacts by default", {
  cr <- HousepiContactRate$new()
  aux_add(cr, by = "age",
          reference = c(from_age = "adult", to_age = "adult"),
          prior = "lognormal", sdlog = 1)
  cr$implement(dummy_df, dummy_factor_info)

  # two-level feature: child, adult -> 4 pairs
  ltab <- cr$logic_tables[[1]]
  expect_equal(nrow(ltab), 4)

  # only 1 reference
  expect_equal(sum(ltab$param_name == ".reference"), 1)
  # only 2 unique non-reference params (cross + child->child)
  expect_length(cr$parameters, 2)

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
  cr <- HousepiContactRate$new()
  aux_add(cr, by = "age")
  cr$implement(dummy_df, dummy_factor_info)

  # automatic .reference
  expect_true(".reference" %in% cr$logic_tables[[1]]$param_name)
})

test_that("Auto mode symmetric = FALSE gives distinct parameters", {
  cr <- HousepiContactRate$new()
  aux_add(cr, by = "age", symmetric = FALSE,
          reference = c(from_age = "adult", to_age = "adult"),
          prior = "lognormal", sdlog = 1)
  cr$implement(dummy_df, dummy_factor_info)

  ltab <- cr$logic_tables[[1]]
  expect_equal(nrow(ltab), 4)
  expect_length(cr$parameters, 3)

  # asymmetric: adult->child and child->adult are distinct params
  ac_row <- subset(ltab, from_age == "adult" & to_age == "child")
  ca_row <- subset(ltab, from_age == "child" & to_age == "adult")
  expect_true(ac_row$param_name != ca_row$param_name)
})

test_that("Non-joint produces independent per-feature tables", {
  cr <- HousepiContactRate$new()

  # 2 features, non-joint -> 2 sub-tables of 4 rows each
  aux_add(cr, by = c("age", "sex"), joint = FALSE,
          reference = c(age = "adult", sex = "M"),
          prior = "lognormal", sdlog = 1)
  cr$implement(dummy_df, dummy_factor_info)

  # two sub-tables
  expect_length(cr$logic_tables, 2L)

  # all age pairs are observed: 2 x 2
  expect_equal(nrow(cr$logic_tables[[1]]), 4L)
  # 1 sex pair is missing: 2 x 2 - 1
  expect_equal(nrow(cr$logic_tables[[2]]), 3L)

  # each sub-table has its own .reference
  for (tbl in cr$logic_tables) {
    expect_equal(sum(tbl$param_name == ".reference"), 1L)
  }

  # symmetric default: 2 + 2 - 1 = 3 non-reference params
  expect_length(cr$parameters, 3L)

  # tidy table merges sub-tables with NA for missing columns
  tidy <- cr$tidy_table
  expect_equal(nrow(tidy), 7L)
})

test_that("Auto mode with two features behaves correctly", {
  cr <- HousepiContactRate$new()

  aux_add(cr, by = c("age", "sex"),
          reference = c(from_age = "adult", from_sex = "M",
                        to_age = "adult", to_sex = "M"),
          prior = "lognormal", sdlog = 1)
  cr$implement(dummy_df, dummy_factor_info)

  # NOTE: 3 observed types: <child,M>, <child,F>, <adult,M>
  # 7 observed pairs: 9 pairs - <child,M><child,F> - (<->)
  ltab <- cr$logic_tables[[1]]
  expect_equal(nrow(ltab), 7L)
  expect_equal(sum(ltab$param_name == ".reference"), 1)

  # NOTE: symmetric pairs by default:
  # 3 heterogeneous: <child,M><child,F> <adult,M><child,F> <adult,M><child,M>
  # 1 homogeneous: <adult,M><adult,M> (reference)
  expect_length(cr$parameters, 3L)
})

########################################
## CONTACT RATES : CUSTOM DEFINITIONS ##
########################################

test_that("Contact rate custom mode defualt behaves correctly", {
  # logic table
  param_name <- c("a2c", "c2a", "c2c")
  strata <- data.frame(from_age = c(rep("adult", 2), rep("child", 2)),
                       to_age = rep(c("adult", "child"), 2),
                       param_name = c(".reference", param_name))
  # custom fine grained parameter override
  overrides <- list(
    list(value = 1, prior = "uniform", min = 1, max = 5,
         proposal = "normal", proposal_scale = 0.5),
    2.34
  )
  fancy_param_1 <- param_name[1] # to override as variable
  fancy_param_2 <- param_name[2] # to override as constant
  default_param <- param_name[3] # default definitions
  names(overrides) <- c(fancy_param_1, fancy_param_2)

  cr <- HousepiContactRate$new()

  expect_error(aux_add(cr, by = "age", strata = strata),
               regexp = "Impossible")

  aux_add(cr, strata = strata, prior = "lognormal", sdlog = 1.8,
          params = overrides)

  # definitions stored with mode sentinel
  expect_equal(cr$defs$mode, "custom")

  # logic table stored
  expect_true(is.data.frame(cr$logic_tables[[1]]))
  expect_equal(nrow(cr$logic_tables[[1]]), 4L)

  # parameters generated for non-reference
  expect_length(cr$parameters, 3L)
  expect_true(all(param_name %in% names(cr$parameters)))

  # default parameter definitions
  expect_equal(cr$parameters[[default_param]]$prior, "lognormal")
  expect_equal(cr$parameters[[default_param]]$sdlog, "1.8|N")
  # overrides as variable
  expect_equal(cr$parameters[[fancy_param_1]]$prior, "uniform")
  expect_equal(cr$parameters[[fancy_param_1]]$min, 1)
  # overrides as constant
  expect_equal(cr$parameters[[fancy_param_2]], 2.34)

  # implement should not clear the logic table
  cr$implement(dummy_df, dummy_factor_info)
  expect_true(is.data.frame(cr$logic_tables[[1]]))
  expect_equal(cr$logic_tables[[1]]$param_name, strata$param_name)

  # tidy table
  tidy <- cr$tidy_table
  expect_true(is.data.frame(tidy))
  expect_equal(nrow(tidy), 4L)
})

test_that("Custom mode with no reference row raises error", {
  strata <- data.frame(
    from_age = c("adult", "adult", "child", "child"),
    to_age   = c("adult", "child", "adult", "child"),
    param_name = c("a_to_a", "a_to_c", "c_to_a", "c_to_c")
  )

  cr <- HousepiContactRate$new()
  expect_error(aux_add(cr, strata = strata),
               regexp = "reference")
})

####################################################
## CONTACT RATES : CUSTOM (SYMMETRIC) DEFINITIONS ##
####################################################

test_that("Custom mode allows symmetric definition automatic fixes", {
  param_name <- c("a_to_c", "c_to_c")
  strata <- data.frame(from_age = c("adult", "adult", "child"),
                       to_age = c("adult", "child", "child"),
                       param_name = c(".reference", param_name))
  overrides <- list(a_to_c = 2)

  cr <- HousepiContactRate$new()
  aux_add(cr, strata = strata, params = overrides)
  cr$implement(dummy_df, dummy_factor_info)

  # symmetric fill adds child->adult with same param_name as adult->child
  filled <- cr$logic_tables[[1]]
  expect_equal(nrow(filled), 4)
  ca_row <- filled[filled$from_age == "child" & filled$to_age == "adult", ]
  expect_equal(ca_row$param_name, "a_to_c")

  # both directions share the override value
  expect_true(all(param_name %in% names(cr$parameters)))
  expect_equal(cr$parameters[["a_to_c"]], 2L)

  # all-NA row (.reference) is its own reverse and must not be duplicated;
  # reverse of dad_child and mum_child are each added once
  strata <- data.frame(from_age = c("adult", "adult", NA),
                       from_sex = c("M", "F", NA),
                       to_age = c("child", "child", NA),
                       to_sex = c(NA, NA, NA),
                       param_name = c("dad_child", "mum_child", ".reference"))
  cr <- HousepiContactRate$new()
  aux_add(cr, strata = strata)

  cr$implement(dummy_df, dummy_factor_info)
  expect_equal(nrow(cr$logic_tables[[1]]), 5L)
})

test_that("Custom mode with symmetric = FALSE does not fill gaps", {
  param_name <- c("a_to_c", "c_to_c")
  strata <- data.frame(from_age = c("adult", "adult", "child"),
                       to_age = c("adult", "child", "child"),
                       param_name = c(".reference", param_name))

  cr <- HousepiContactRate$new()
  aux_add(cr, strata = strata, symmetric = FALSE)
  cr$implement(dummy_df, dummy_factor_info)

  expect_equal(nrow(cr$logic_tables[[1]]), 3L)
})

################################################
## CONTACT RATES : CUSTOM (JOINT) DEFINITIONS ##
################################################

test_that("joint = TRUE explicit matches default behavior", {
  ctr_default <- HousepiContactRate$new()
  aux_add(ctr_default, by = c("age", "sex"),
          reference = c(age = "adult", sex = "M"),
          prior = "lognormal", sdlog = 1)
  ctr_default$implement(dummy_df, dummy_factor_info)

  ctr_joint <- HousepiContactRate$new()
  aux_add(ctr_joint, by = c("age", "sex"),
          joint = TRUE,
          reference = c(age = "adult", sex = "M"),
          prior = "lognormal", sdlog = 1)
  ctr_joint$implement(dummy_df, dummy_factor_info)

  # same single sub-table
  expect_length(ctr_default$logic_tables, 1)
  expect_length(ctr_joint$logic_tables, 1)
  expect_equal(ctr_default$logic_tables[[1]],
               ctr_joint$logic_tables[[1]])
})

#########################################
## CONTACT RATES : RULE GROUPING LOGIC ##
#########################################

test_that("compute_rule_groups: ambiguous subset raises error", {
  tbl <- data.frame(
    from_age = c("child", "child", "child"),
    from_sex = c("M", NA, NA),
    to_age = c(NA, "child", NA),
    param_name = c("metal_gear_1", "metal_gear_2", "big_boss")
  )
  expect_error(compute_rule_groups(tbl), regexp = "[Aa]mbiguous")
})

test_that("Custom mode absorbs all-NA default into single group", {
  # Example 1: 4 rules in one group, last one (all-NA) is the catch-all
  strata <- data.frame(
    from_member = c("student", "mother", "father", NA),
    to_member = c("student", "student", "mother", NA),
    param_name = c("ctr_child_child", "ctr_mother_child",
                   ".reference", "ctr_cross_gen")
  )
  cr <- HousepiContactRate$new()
  aux_add(cr, strata = strata, symmetric = FALSE)
  tidy <- cr$tidy_table
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
  cr <- HousepiContactRate$new()
  aux_add(cr, strata = strata, symmetric = FALSE)
  expect_equal(length(unique(cr$tidy_table$rule_group)), 2L)
})

test_that("Auto mode joint = FALSE yields one rule_group per sub-table", {
  cr <- HousepiContactRate$new()
  aux_add(cr, by = c("age", "sex"), joint = FALSE,
          reference = c(age = "adult", sex = "M"))
  cr$implement(dummy_df, dummy_factor_info)
  tidy <- cr$tidy_table
  expect_equal(length(unique(tidy$rule_group)), 2L)
})

test_that("Auto mode joint = TRUE yields a single rule_group", {
  cr <- HousepiContactRate$new()
  aux_add(cr, by = c("age", "sex"),
          reference = c(from_age = "adult", from_sex = "M",
                        to_age = "adult", to_sex = "M"))
  cr$implement(dummy_df, dummy_factor_info)
  tidy <- cr$tidy_table
  expect_equal(length(unique(tidy$rule_group)), 1L)
})

test_that("Non-joint sub-tables have correct columns", {
  cr <- HousepiContactRate$new()
  aux_add(cr, by = c("age", "sex"),
          joint = FALSE,
          prior = "lognormal", sdlog = 1)
  cr$implement(dummy_df, dummy_factor_info)

  # first table has from_age/to_age
  t1 <- cr$logic_tables[[1]]
  expect_true(all(c("from_age", "to_age") %in% colnames(t1)))
  expect_false("from_sex" %in% colnames(t1))

  # second table has from_sex/to_sex
  t2 <- cr$logic_tables[[2]]
  expect_true(all(c("from_sex", "to_sex") %in% colnames(t2)))
  expect_false("from_age" %in% colnames(t2))
})

#####################################################
## CONTACT RATES : C++ RESOLUTION (best_match.h)   ##
#####################################################

# These tests drive update_best_match / resolve_best_match in src/rates.cpp
# through the full model build, so the C++ matrix can be read back and
# inspected via export_house_snapshot.

# Cell-by-cell readout helper.
ctr_cell_ <- function(hhsnapshot, from, to, key = "member") {
  ids <- hhsnapshot$data_wide[[key]]
  hhsnapshot$contact_rate[match(from, ids), match(to, ids)]
}

test_that("Most-specific rule wins over reference (single feature)", {
  # house with one of each member type
  hh <- data.frame(
    id = c("m1", "m2", "m3", "m4"),
    house_id = "h1",
    status = c(0L, 0L, 0L, 0L),
    member = c("father", "mother", "student", "sibling"),
    stringsAsFactors = FALSE
  )

  # one specific pair gets a non-reference param; everything else falls
  # through to .reference (specificity 0, hits update_best_match's empty
  # branch then never gets beaten).
  strata <- data.frame(
    from_member = c("father", NA),
    to_member = c("mother", NA),
    param_name = c("ctr_fm", ".reference"),
    stringsAsFactors = FALSE
  )

  hhmodel <- new_housepi_model(hh, "chain_binomial")
  hhmodel$add_parameter_template()
  hhmodel$add_contact_rate(strata = strata, symmetric = FALSE,
                           params = list(ctr_fm = 2))
  snap <- hhmodel$export_house_snapshot("h1")

  # the explicit pair wins
  expect_equal(ctr_cell_(snap, "father", "mother"), 2)
  # symmetric = FALSE: reverse falls back to reference
  expect_equal(ctr_cell_(snap, "mother", "father"), 1)
  # unrelated pair falls back to reference (spec=0)
  expect_equal(ctr_cell_(snap, "student", "sibling"), 1)
  # diagonal stays 0 (no self-contact)
  expect_equal(unname(diag(snap$contact_rate)), rep(0., nrow(hh)))
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

  hhmodel <- new_housepi_model(hh, "chain_binomial")
  hhmodel$add_parameter_template()
  hhmodel$add_contact_rate(
    strata = strata, symmetric = FALSE, reference = ".none",
    params = list(ctr_aa = 3, ctr_class_aa = 5)
  )
  snap <- hhmodel$export_house_snapshot("h1")

  # m1 (adult,a) → m2 (adult,a): BOTH rules match → 3 * 5 = 15
  # (init_product over two Attribute* — rates.cpp's multi-group loop)
  expect_equal(ctr_cell_(snap, "m1", "m2", key = "id"), 15)
  # m1 (adult,a) → m4 (adult,b): only age rule matches → 3
  expect_equal(ctr_cell_(snap, "m1", "m4", key = "id"), 3)
  # m1 (adult,a) → m3 (child,b): neither rule matches → reference (1)
  expect_equal(ctr_cell_(snap, "m1", "m3", key = "id"), 1)
  # diagonal is 0
  expect_equal(unname(diag(snap$contact_rate)), rep(0., nrow(hh)))
})

