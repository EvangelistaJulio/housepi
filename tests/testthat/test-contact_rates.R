library(testthat)

aux_add <- function(self, by = NA_character_,
                    reference = NULL, prefix = "ctr",
                    joint = TRUE, symmetric = TRUE,
                    logic_table = NULL, params = NULL,
                    prior = NA_character_, ...,
                    min = NA_real_, max = NA_real_,
                    proposal = NA_character_,
                    proposal_scale = NA_real_,
                    .warn_over = TRUE) {
  self$add_quo(by = by, prefix = prefix,
               reference = reference,
               joint = joint, symmetric = symmetric,
               logic_table = logic_table,
               params = params,
               prior = prior, rlang::enquos(...),
               min = rlang::enquo(min),
               max = rlang::enquo(max),
               proposal = proposal,
               proposal_scale = proposal_scale,
               .warn_over = .warn_over)
  invisible(self)
}

#################################
## CONTACT RATES : PAIR COUNTS ##
#################################

test_that("Contact rates pair counts are computed correctly", {
  df_hh <- data.frame(house_id = c(rep("hh1", 2), rep("hh2", 2)),
                      age = rep(c("child", "adult"), 2),
                      stringsAsFactors = FALSE)
  tbl <- contact_rate_logic_table(df_hh, "age")

  ca <- tbl[tbl$from_age == "child" & tbl$to_age == "adult", "n"]
  ac <- tbl[tbl$from_age == "adult" & tbl$to_age == "child", "n"]
  cc <- tbl[tbl$from_age == "child" & tbl$to_age == "child", "n"]
  aa <- tbl[tbl$from_age == "adult" & tbl$to_age == "adult", "n"]

  # total: child->adult = 2, adult->child = 2
  # no child->child = 0, adult->adult = 0
  expect_equal(ca, 2L)
  expect_equal(ac, 2L)
  expect_equal(cc, 0L)
  expect_equal(aa, 0L)
})

######################################
## CONTACT RATES : AUTO DEFINITIONS ##
######################################

# a dummy dataset and relative information for factorization
dummy_df <- data.frame(house_id = 1,
                       age = c(rep("child", 2), rep("adult", 2)),
                       sex = c("M", "F", "M", "M"),
                       stringsAsFactors = FALSE)
dummy_factor_info <- list(age = list(levels = c("child", "adult")),
                          sex = list(levels = c("M", "F")))

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

test_that("Auto mode with two features behaves correctly", {
  cr <- HousepiContactRate$new()

  aux_add(cr, by = c("age", "sex"),
          reference = c(from_age = "adult", from_sex = "M",
                        to_age = "adult", to_sex = "M"),
          prior = "lognormal", sdlog = 1)
  cr$implement(dummy_df, dummy_factor_info)

  # NOTE: dummy_df observation:
  # <child,M>, <child,F>, <adult,M>
  # 3 observed types -> 9 pairs
  ltab <- cr$logic_tables[[1]]
  expect_equal(nrow(ltab), 9)
  expect_equal(sum(ltab$param_name == ".reference"), 1)

  # NOTE: symmetric by default):
  # unique parameters = 3!/2! pairs + 3 self-pairs - 1 reference = 5
  expect_length(cr$parameters, 5)
})

test_that("Non-joint produces independent per-feature tables", {
  cr <- HousepiContactRate$new()

  # 2 features, non-joint -> 2 sub-tables of 4 rows each
  aux_add(cr, by = c("age", "sex"), joint = FALSE,
          reference = c(age = "adult", sex = "M"),
          prior = "lognormal", sdlog = 1)
  cr$implement(dummy_df, dummy_factor_info)

  # two sub-tables
  expect_length(cr$logic_tables, 2)

  # each sub-table has 4 rows (2 from x 2 to)
  expect_equal(nrow(cr$logic_tables[[1]]), 4)
  expect_equal(nrow(cr$logic_tables[[2]]), 4)

  # each sub-table has its own .reference
  for (tbl in cr$logic_tables) {
    expect_equal(sum(tbl$param_name == ".reference"), 1)
  }

  # symmetric default: 2 + 2 = 4 non-reference params
  expect_length(cr$parameters, 4)

  # tidy table merges sub-tables with NA for missing columns
  tidy <- cr$tidy_table
  expect_equal(nrow(tidy), 8)
})

########################################
## CONTACT RATES : CUSTOM DEFINITIONS ##
########################################

test_that("Contact rate custom mode defualt behaves correctly", {
  # logic table
  param_name <- c("ctr_a2c", "ctr_c2a", "ctr_c2c")
  lt <- data.frame(from_age = c(rep("adult", 2), rep("child", 2)),
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

  expect_error(aux_add(cr, by = "age", logic_table = lt),
               regexp = "Impossible")

  aux_add(cr, logic_table = lt, prior = "lognormal", sdlog = 1.8,
          params = overrides)

  # definitions stored with mode sentinel
  expect_equal(cr$defs$mode, "custom")

  # logic table stored
  expect_true(is.data.frame(cr$logic_tables[[1]]))
  expect_equal(nrow(cr$logic_tables[[1]]), 4)

  # parameters generated for non-reference
  expect_length(cr$parameters, 3)
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
  expect_equal(cr$logic_tables[[1]]$param_name, lt$param_name)

  # tidy table
  tidy <- cr$tidy_table
  expect_true(is.data.frame(tidy))
  expect_equal(nrow(tidy), 4)
})

test_that("Custom mode allows symmetric definition automatic fixes", {
  param_name <- c("ctr_a2c", "ctr_c2c")
  lt <- data.frame(from_age = c("adult", "adult", "child"),
                   to_age = c("adult", "child", "child"),
                   param_name = c(".reference", param_name))
  overrides <- list(ctr_a2c = 2)

  cr <- HousepiContactRate$new()
  aux_add(cr, logic_table = lt, params = overrides)
  cr$implement(dummy_df, dummy_factor_info)

  # symmetric fill adds child->adult with same param_name as adult->child
  filled <- cr$logic_tables[[1]]
  expect_equal(nrow(filled), 4)
  ca_row <- filled[filled$from_age == "child" & filled$to_age == "adult", ]
  expect_equal(ca_row$param_name, "ctr_a2c")

  # both directions share the override value
  expect_true(all(param_name %in% names(cr$parameters)))
  expect_equal(cr$parameters[["ctr_a2c"]], 2)
})

test_that("Custom mode with no .reference row raises error", {
  lt <- data.frame(
    from_age = c("adult", "adult", "child", "child"),
    to_age   = c("adult", "child", "adult", "child"),
    param_name = c("ctr_a2a", "ctr_a2c", "ctr_c2a", "ctr_c2c"),
    stringsAsFactors = FALSE
  )

  cr <- HousepiContactRate$new()
  expect_error(
    aux_add(cr, logic_table = lt, prior = "lognormal", sdlog = 1),
    regexp = "reference"
  )
})

###########################################################################
# 4. joint = FALSE: independent per-feature tables
###########################################################################

test_that("joint = TRUE explicit matches default behavior", {
  ctr_default <- HousepiContactRate$new()
  aux_add(ctr_default, by = c("age", "sex"), prefix = "cr",
          reference = c(age = "adult", sex = "M"),
          prior = "lognormal", sdlog = 1)
  ctr_default$implement(dummy_df, dummy_factor_info)

  ctr_joint <- HousepiContactRate$new()
  aux_add(ctr_joint, by = c("age", "sex"), prefix = "cr",
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

test_that("Non-joint sub-tables have correct columns", {
  cr <- HousepiContactRate$new()
  aux_add(cr, by = c("age", "sex"), prefix = "cr",
          joint = FALSE,
          prior = "lognormal", sdlog = 1)
  cr$implement(dummy_df, dummy_factor_info)

  # first table has from_age/to_age
  t1 <- cr$logic_tables[[1]]
  expect_true(all(c("from_age", "to_age") %in%
                    colnames(t1)))
  expect_false("from_sex" %in% colnames(t1))

  # second table has from_sex/to_sex
  t2 <- cr$logic_tables[[2]]
  expect_true(all(c("from_sex", "to_sex") %in% colnames(t2)))
  expect_false("from_age" %in% colnames(t2))
})

###########################################################################
# 5. Custom mode symmetric gap-filling
###########################################################################

test_that("Custom mode with symmetric fills missing reverse pairs", {
  # only define adult->child; child->adult should be auto-filled
  lt <- data.frame(
    from_age = c("adult", "adult", "child"),
    to_age   = c("adult", "child", "child"),
    param_name = c(".reference", "ctr_a2c", "ctr_c2c"),
    stringsAsFactors = FALSE
  )

  cr <- HousepiContactRate$new()
  aux_add(cr, logic_table = lt, symmetric = TRUE,
          prior = "lognormal", sdlog = 1)
  cr$implement(dummy_df, dummy_factor_info)

  filled <- cr$logic_tables[[1]]
  # child->adult should now exist with same param as adult->child
  expect_equal(nrow(filled), 4)
  ca_row <- filled[filled$from_age == "child" & filled$to_age == "adult", ]
  expect_equal(nrow(ca_row), 1)
  expect_equal(ca_row$param_name, "ctr_a2c")
})

test_that("Custom mode with symmetric = FALSE does not fill gaps", {
  lt <- data.frame(
    from_age = c("adult", "adult", "child"),
    to_age   = c("adult", "child", "child"),
    param_name = c(".reference", "ctr_a2c", "ctr_c2c"),
    stringsAsFactors = FALSE
  )

  cr <- HousepiContactRate$new()
  aux_add(cr, logic_table = lt, symmetric = FALSE,
          prior = "lognormal", sdlog = 1)
  cr$implement(dummy_df, dummy_factor_info)

  expect_equal(nrow(cr$logic_tables[[1]]), 3)
})
