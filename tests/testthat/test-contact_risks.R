library(testthat)

# model over two covariates
mk_risk <- local({
  hhmodel <- data.frame(house_id = 1L, status = 0L,
                        age = c(rep("child", 2), rep("adult", 2)),
                        sex = c("M", "F", "M", "M")) |>
    new_housepi_model("final-size")
  list(model = hhmodel)
})

# house holding one of each member type, with template parameters
mk_house <- local({
  model <- data.frame(id = c("m1", "m2", "m3", "m4"),
                      house_id = "h1", status = 0L,
                      member = c("father", "mother", "student", "sibling")) |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")
  list(model = model)
})

# cell-by-cell readout of the contact risk matrix of a household snapshot
ctr_cell <- function(snapshot, from, to, key = "member") {
  ids <- snapshot$data_wide[[key]]
  snapshot$contact_risks[match(from, ids), match(to, ids)]
}

# one pinned father-to-mother pair over `mk_house`
pin_ctr_strata <- function(value = c(0.4, NA), const = c(TRUE, NA), ...) {
  strata <- data.frame(from_member = c("father", NA),
                       to_member = c("mother", NA),
                       param_name = c("ctr_fm", ".reference"),
                       value = value, const = const)
  mk_house$model  |>
    add_contact_risk(strata = strata, ...)
}

#################################
## CONTACT RISKS : PAIR COUNTS ##
#################################

test_that("Contact risks pair counts are computed correctly", {
  tbl <- data.frame(house_id = c(rep("h1", 2L), rep("h2", 2L)),
                    age = rep(c("child", "adult"), 2L)) |>
    contact_risk_logic_table("age")
  count <- function(from, to) tbl$n[tbl$from_age == from & tbl$to_age == to]

  # total: child->adult = 2, adult->child = 2, no same-age pair
  expect_equal(count("child", "adult"), 2L)
  expect_equal(count("adult", "child"), 2L)
  expect_equal(count("child", "child"), integer(0))
  expect_equal(count("adult", "adult"), integer(0))

  # no household structure, no pair to count
  contact_risk_logic_table(data.frame(age = c("child", "adult")), "age") |>
    expect_null()
})

test_that("Pair counts drop the member meeting itself", {
  # 2 children in one household imply 2 ordered child->child pairs, not 4
  tbl <- data.frame(house_id = c(rep("h1", 3), rep("h2", 2)),
                    age = c(rep("child", 2), rep("adult", 3))) |>
    contact_risk_logic_table("age")
  count <- function(from, to) tbl$n[tbl$from_age == from & tbl$to_age == to]

  expect_equal(count("child", "child"), 2L)
  expect_equal(count("adult", "adult"), 2L)
  expect_equal(count("child", "adult"), 2L)
  expect_equal(count("adult", "child"), 2L)
})

######################################
## CONTACT RISKS : AUTO DEFINITIONS ##
######################################

test_that("Auto mode imposes symmetric contacts by default", {
  risks_of <- function(...) {
    mk_risk$model |>
      add_contact_risk(by = "age", prior = "lognormal", sdlog = 1., ...,
                       reference = c(from_age = "adult", to_age = "adult")) |>
      get_contact_risks()
  }
  tidy <- risks_of()

  # two-level covariate: child, adult -> 4 pairs
  expect_equal(nrow(tidy), 4L)

  # only 1 reference
  expect_equal(sum(tidy$param_name == ".reference"), 1)
  # only 2 unique non-reference params (cross + child->child)
  expect_length(setdiff(unique(tidy$param_name), ".reference"), 2)

  # reference row is adult->adult
  ref_row <- tidy[tidy$param_name == ".reference", ]
  ref_row$from_age |> expect_equal("adult")
  ref_row$to_age |> expect_equal("adult")

  # symmetric: adult->child and child->adult share the same parameter
  expect_equal(
    subset(tidy, from_age == "adult" & to_age == "child")$param_name,
    subset(tidy, from_age == "child" & to_age == "adult")$param_name
  )

  # symmetric = FALSE gives distinct parameters to the two directions
  tidy <- risks_of(symmetric = FALSE)
  expect_equal(nrow(tidy), 4L)
  expect_length(setdiff(unique(tidy$param_name), ".reference"), 3L)
  expect_false(identical(
    subset(tidy, from_age == "adult" & to_age == "child")$param_name,
    subset(tidy, from_age == "child" & to_age == "adult")$param_name
  ))
})

test_that("Auto mode reference is opt-in and picks most abundant pair", {
  tidy <- mk_risk$model |>
    add_contact_risk(by = "age", reference = ".auto") |>
    get_contact_risks()

  # automatic .reference on the most abundant pair,
  # tagged in both directions since the contact is symmetric
  ref <- subset(tidy, param_name == ".reference")
  expect_gte(nrow(ref), 1L)
  expect_true(all(ref$n == max(tidy$n)))

  # absolute by default: every pair keeps its own parameter
  ctr <- mk_risk$model |>
    add_contact_risk(by = "age") |>
    get_contact_risks()
  nrow(subset(ctr, param_name == ".reference")) |>
    expect_equal(0L)
})

test_that("Non-joint produces independent per-covariate tables", {
  # 2 covariates, non-joint -> 2 rule groups
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
    sum(tidy$param_name[tidy$rule_group == grp] == ".reference") |>
      expect_equal(1L)
  }

  # symmetric default: 2 + 2 - 1 = 3 non-reference params
  setdiff(unique(tidy$param_name), ".reference") |>
    expect_length(3L)

  # first sub-table carries age, and leaves sex unset
  first <- subset(tidy, rule_group == 1)
  expect_true(all(!is.na(first$from_age) & !is.na(first$to_age)))
  expect_true(all(is.na(first$from_sex) & is.na(first$to_sex)))

  # second sub-table carries sex, and leaves age unset
  second <- subset(tidy, rule_group == 2)
  expect_true(all(!is.na(second$from_sex) & !is.na(second$to_sex)))
  expect_true(all(is.na(second$from_age) & is.na(second$to_age)))
})

test_that("Auto mode with two joint covariates behaves correctly", {
  tidy <- mk_risk$model |>
    add_contact_risk(by = c("age", "sex"), joint = TRUE,
                     reference = c(from_age = "adult", from_sex = "M",
                                   to_age = "adult", to_sex = "M"),
                     prior = "lognormal", sdlog = 1.) |>
    get_contact_risks()

  # NOTE: 3 observed types: <child,M>, <child,F>, <adult,M>
  # 7 observed pairs: 9 pairs - <child,M><child,F> - (<->)
  expect_equal(nrow(tidy), 7L)
  expect_equal(length(unique(tidy$rule_group)), 1L)
  expect_equal(sum(tidy$param_name == ".reference"), 1)

  # NOTE: symmetric pairs by default:
  # 3 heterogeneous: <child,M><child,F> <adult,M><child,F> <adult,M><child,M>
  # 1 homogeneous: <adult,M><adult,M> (reference)
  expect_length(setdiff(unique(tidy$param_name), ".reference"), 3L)
})

test_that("Absoluteness follows the reference across every definition mode", {
  # the default prior tracks whether a risk carries absolute values
  risk_prior <- function(...) {
    params <- get_parameters(add_contact_risk(mk_risk$model, ...))
    unique(params$prior[!params$const])
  }
  tagged <- data.frame(from_age = c("adult", "adult"),
                       to_age = c("adult", "child"),
                       param_name = c(".reference", "a_to_c"))
  untagged <- tagged
  untagged$param_name[1] <- "a_to_a"

  # auto mode: absolute unless a reference is resolved
  risk_prior(by = "age") |>
    expect_equal("uniform")
  risk_prior(by = "age", reference = ".auto") |>
    expect_equal("lognormal")
  risk_prior(by = "age", reference = c(age = "adult")) |>
    expect_equal("lognormal")
  # custom mode: the tag alone decides
  risk_prior(strata = untagged) |>
    expect_equal("uniform")
  risk_prior(strata = tagged) |>
    expect_equal("lognormal")
})

test_that("Every rule group keeps its anchoring in auto", {
  risk_refer <- function(...) {
    sarscov2_2021_icofs |>
      dplyr::filter(house_id %in% utils::head(unique(house_id), 50)) |>
      new_housepi_model("final-size") |>
      add_contact_risk(by = c("age", "sex"), ...)
  }

  # reference anchors and sets lognormal prior
  hhmodel <- risk_refer(reference = c(from_age = "adult", to_age = "adult"))
  params <- get_parameters(hhmodel)
  contact <- params[grepl("^(age|sex)_", params$name), ]
  expect_true(all(contact["prior"] == "lognormal"))

  # joint = FALSE by default so each covariate owns a rule group
  get_contact_risks(hhmodel)$rule_group |>
    expect_setequal(c(1L, 2L))

  # no reference makes it scale free and and sets uniform prior
  params <- risk_refer() |>
    get_parameters()
  contact <- params[grepl("^(age|sex)_", params$name), ]
  expect_true(all(contact["prior"] == "uniform"))
})

test_that("Spare prior arguments on a contact risk are reported", {
  mk_risk$model |>
    add_contact_risk(by = "age", sdlog = 1.8) |>
    expect_error("[Ss]pare prior arguments.*contact risk")
})

test_that("An invalid contact reference is rejected", {
  mk_risk$model |>
    add_contact_risk(by = "age",
                     reference = c(from_age = "elderly", to_age = "adult")) |>
    expect_error("age = .elderly. <-> age = .adult. is never observed")
  mk_risk$model |>
    add_contact_risk(by = "age", reference = c(age_from = "adult")) |>
    expect_error("[Uu]nrecognized choice for reference: .age_from.")
})

test_that("Contact risks expand on their own covariates only", {
  hhdata <- sarscov2_2021_icofs |>
    dplyr::select(-start_date) |>
    dplyr::filter(house_id %in% utils::head(unique(house_id), 50))

  # the covariate itself must be a genuine data column
  hhdata |>
    dplyr::select(-age) |>
    new_housepi_model("final-size") |>
    add_contact_risk(by = "age") |>
    expect_error("[Uu]nknown covariate.*age")

  # without data the read surface warns and keeps the table shape
  expect_warning(
    risks <- new_housepi_model(model = "final-size") |>
      add_contact_risk(by = "member") |>
      get_contact_risks(),
    regexp = "without data"
  )
  nrow(risks) |> expect_equal(0L)
  risks |> expect_named(c("param_name", "n", "rule_group"))
})

########################################
## CONTACT RISKS : CUSTOM DEFINITIONS ##
########################################

test_that("Contact risk custom mode default behaves correctly", {
  # logic table
  strata <- data.frame(from_age = c(rep("adult", 2L), rep("child", 2L)),
                       to_age = rep(c("adult", "child"), 2L),
                       param_name = c(".reference", "a2c", "c2a", "c2c"))
  # a2c overridden as variable, c2a as constant, c2c left to the defaults
  overrides <- list(
    a2c = list(value = 1, prior = "uniform", min = 1., max = 5.,
               proposal = "rw", proposal_scale = 0.5),
    c2a = 2.34
  )

  hhmodel <- mk_risk$model |>
    add_contact_risk(strata = strata, prior = "lognormal", sdlog = 1.8,
                     params = overrides)

  # the custom logic table drives the tidy table
  tidy <- get_contact_risks(hhmodel)
  nrow(tidy) |> expect_equal(4L)
  tidy$param_name |> expect_setequal(strata$param_name)

  # parameters generated for non-reference
  get_parameters(hhmodel) |>
    getElement("name") |>
    expect_setequal(c("a2c", "c2a", "c2c"))

  # default parameter definitions
  default <- get_parameters(hhmodel, "c2c")
  default$prior |> expect_equal("lognormal")
  default$args[[1]] |> expect_equal(c(sdlog = 1.8, meanlog = 0))
  # overrides as variable
  variable <- get_parameters(hhmodel, "a2c")
  variable$prior |> expect_equal("uniform")
  variable$min |> expect_equal(1)
  # overrides as constant
  constant <- get_parameters(hhmodel, "c2a")
  constant$value |> expect_equal(2.34)
  expect_true(constant$const)
})

test_that("A reference has no meaning against a custom contact table", {
  strata <- data.frame(from_age = c("adult", "adult", "child", "child"),
                       to_age   = c("adult", "child", "adult", "child"),
                       param_name = c("a_to_a", "a_to_c", "c_to_a", "c_to_c"))
  mk_risk$model |>
    add_contact_risk(strata = strata, reference = ".auto") |>
    expect_error("[Nn]ot used with a custom")
})

test_that("Custom contact strata carry at most one reference per rule group", {
  # a symmetric pair reaches one cell through either direction
  strata <- data.frame(from_age = c("adult", "child", "child"),
                       to_age   = c("child", "adult", "child"),
                       param_name = c(".reference", ".reference", "c_to_c"))
  mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = TRUE) |>
    expect_no_error()
  # without symmetry the two directions are two competing scales
  mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    expect_error("[Rr]ows.*one rule group")
})

test_that("Custom mode fills the symmetric counterpart of a pair", {
  strata <- data.frame(from_age = c("adult", "adult", "child"),
                       to_age = c("adult", "child", "child"),
                       param_name = c(".reference", "a_to_c", "c_to_c"))

  hhmodel <- mk_risk$model |>
    add_contact_risk(strata = strata, params = list(a_to_c = 2))

  # symmetric fill adds child->adult with same param_name as adult->child
  filled <- get_contact_risks(hhmodel)
  nrow(filled) |> expect_equal(4L)
  filled$param_name[filled$from_age == "child" & filled$to_age == "adult"] |>
    expect_equal("a_to_c")

  # both directions share the override value
  get_parameters(hhmodel) |>
    getElement("name") |>
    expect_setequal(c("a_to_c", "c_to_c"))
  a_to_c <- get_parameters(hhmodel, "a_to_c")
  expect_equal(a_to_c$value, 2)
  expect_true(a_to_c$const)

  # symmetric = FALSE does not fill gaps
  mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    get_contact_risks() |>
    nrow() |>
    expect_equal(3L)
})

test_that("An unobserved contact level is dropped", {
  strata <- data.frame(from_member = c("mother", NA),
                       to_member = c("grandparent", NA),
                       param_name = c("rho", ".reference"))
  expect_warning(
    hhmodel <- mk_house$model |>
      add_contact_risk(strata = strata, symmetric = FALSE,
                       params = list(rho = 2)),
    "[Uu]nknown definition.*grandparent"
  )

  # since mother->grandparent is not present, `rho`` is not present too
  risks <- get_contact_risks(hhmodel)
  expect_setequal(risks$param_name, ".reference")
  expect_false("rho" %in% get_parameters(hhmodel)$name)

  # every pair reads the reference, the mother's contacts included
  snapshot <- inspect_house_state(hhmodel, "h1")
  ctr_cell(snapshot, "mother", "student") |> expect_equal(1.)
  ctr_cell(snapshot, "mother", "father") |> expect_equal(1.)
  ctr_cell(snapshot, "student", "sibling") |> expect_equal(1.)
})

test_that("A wildcard contact rule counts every pair it reaches", {
  # 4 members: 12 ordered pairs, one of them father -> mother
  strata <- data.frame(from_member = c("father", NA),
                       to_member = c("mother", NA),
                       param_name = c("rho", ".reference"))
  hhmodel <- mk_house$model |>
    add_contact_risk(strata = strata, symmetric = FALSE,
                     params = list(rho = 2))
  risks <- hhmodel |>
    get_contact_risks()

  ctr <- inspect_house_state(hhmodel, "h1")$contact_risks
  risks$n[risks$param_name == "rho"] |> expect_equal(sum(ctr == 2))
  risks$n[risks$param_name == ".reference"] |> expect_equal(sum(ctr == 1))
})

test_that("A dropped contact rule leaves no mirrored counterpart", {
  # the symmetric fill mirrors what survives: were the grandparent rule kept,
  # its two halves would tie on the mother meeting herself
  strata <- data.frame(from_member = c("mother", "mother", NA),
                       to_member = c("student", "grandparent", NA),
                       param_name = c("rho", "rho", ".reference"))
  expect_warning(
    hhmodel <- mk_house$model |>
      add_contact_risk(strata = strata, params = list(rho = 2)),
    "[Uu]nknown definition.*grandparent"
  )

  # the observed pair stands, in both directions
  snapshot <- inspect_house_state(hhmodel, "h1")
  expect_equal(ctr_cell(snapshot, "mother", "student"), 2)
  expect_equal(ctr_cell(snapshot, "student", "mother"), 2)
  # every other pair falls through to the reference
  expect_equal(ctr_cell(snapshot, "father", "sibling"), 1)
})

test_that("A contact table emptied by unknown levels defines no risk", {
  # the only rule names a level no household carries: nothing is left to
  # stratify, and the model ends up carrying no contact risk at all
  expect_warning(
    hhmodel <- mk_house$model |>
      add_contact_risk(strata = data.frame(from_member = "grandparent",
                                           to_member = "grandparent",
                                           param_name = "rho")),
    "[Uu]nknown definition.*grandparent"
  )

  expect_equal(nrow(get_contact_risks(hhmodel)), 0L)
  expect_false("rho" %in% get_parameters(hhmodel)$name)
})

#########################################
## CONTACT RISKS : RULE GROUPING LOGIC ##
#########################################

test_that("An ambiguous subset of rules raises an error", {
  data.frame(from_age = c("child", "child", "child"),
             from_sex = c("M", NA, NA),
             to_age = c(NA, "child", NA),
             param_name = c("rho", "rho_2", "rho_3")) |>
    compute_rule_groups() |>
    expect_error("[Aa]mbiguous")
})

test_that("Custom mode absorbs all-NA default into single group", {
  # 4 rules in one group, last one (all-NA) is the catch-all
  strata <- data.frame(
    from_member = c("student", "mother", "father", NA),
    to_member = c("student", "student", "mother", NA),
    param_name = c("ctr_child_child", "ctr_mother_child",
                   ".reference", "ctr_cross_gen")
  )
  data.frame(house_id = 1, status = 0,
             member = c("student", "mother", "father", "student")) |>
    new_housepi_model("final-size") |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    get_contact_risks() |>
    getElement("rule_group") |>
    unique() |>
    expect_length(1L)
})

test_that("Rule groups follow the maximal patterns of the strata", {
  # disjoint age vs sex patterns -> 2 distinct rule_groups
  strata <- data.frame(
    from_age = c("child", NA, NA),
    from_sex = c(NA, "M", "F"),
    to_age = c("child", NA, NA),
    to_sex = c(NA, "M", "F"),
    param_name = c("age_child_to_age_child", ".reference", "sex_F_to_sex_F")
  )
  mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    get_contact_risks() |>
    getElement("rule_group") |>
    unique() |>
    expect_length(2L)
})

##################################
## CONTACT RISKS : END-TO-END   ##
##################################

test_that("Most-specific rule wins over reference (single covariate)", {
  # one specific pair gets a non-reference param; everything else falls
  # through to .reference (specificity 0, never beaten)
  strata <- data.frame(from_member = c("father", NA),
                       to_member = c("mother", NA),
                       param_name = c("ctr_fm", ".reference"))
  snapshot <- mk_house$model |>
    add_contact_risk(strata = strata, symmetric = FALSE,
                     params = list(ctr_fm = 2)) |>
    inspect_house_state("h1")

  # the explicit pair wins
  expect_equal(ctr_cell(snapshot, "father", "mother"), 2)
  # symmetric = FALSE: reverse falls back to reference
  expect_equal(ctr_cell(snapshot, "mother", "father"), 1)
  # unrelated pair falls back to reference
  expect_equal(ctr_cell(snapshot, "student", "sibling"), 1)
  # diagonal stays 0 (no self-contact)
  expect_equal(unname(diag(snapshot$contact_risks)), rep(0., 4))
})

test_that("Rules tied at equal specificity are rejected as ambiguous", {
  # student -> mother is reached by two half-keyed rules, one specificity each,
  # and the fully-keyed pair makes {from, to} the only maximal pattern: the
  # three rules share a group, so the engine has to compare the tied two
  strata <- data.frame(from_member = c("student", NA, "father"),
                       to_member = c(NA, "mother", "mother"),
                       param_name = c("rho", "rho_2", "rho_3"))
  mk_house$model |>
    add_contact_risk(strata = strata, symmetric = FALSE,
                     params = list(rho = 2, rho_2 = 3, rho_3 = 4)) |>
    inspect_house_state("h1") |>
    expect_error("[Aa]mbiguous contact risk.*tied rules")
})

test_that("A wildcard rule survives a missing cell in its covariate column", {
  # the catch-all row keys on nothing, so it must reach every pair, including
  # the ones whose covariate cell carries a value
  snapshot <- mk_house$model |>
    add_data_wide(transform(mk_house$model$data_wide,
                            member = c("father", "mother", "student", NA))) |>
    add_contact_risk(strata = data.frame(from_member = c("father", NA),
                                         to_member = c("mother", NA),
                                         param_name = c("ctr_fm", "ctr_any")),
                     symmetric = FALSE,
                     params = list(ctr_fm = 2, ctr_any = 3)) |>
    inspect_house_state("h1")

  # the specific pair still wins on specificity
  expect_equal(ctr_cell(snapshot, "m1", "m2", key = "id"), 2)
  # the catch-all reaches a pair of members that both carry a role
  expect_equal(ctr_cell(snapshot, "m3", "m1", key = "id"), 3)
  # and the member with no role reads it too
  expect_equal(ctr_cell(snapshot, "m4", "m1", key = "id"), 3)
  expect_equal(ctr_cell(snapshot, "m1", "m4", key = "id"), 3)
})

test_that("Distinct rule_groups multiply at the cell", {
  # two disjoint maximal patterns ({from_age,to_age}, {from_class,to_class})
  # make two rule groups whose winners multiply into one cell. The table tags
  # no .reference, and unmatched cells read the reference value 1
  snapshot <- data.frame(
    id = c("m1", "m2", "m3", "m4"), house_id = "h1", status = 0L,
    age = c("adult", "adult", "child", "adult"),
    class = c("a", "a", "b", "b")
  ) |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    add_contact_risk(
      strata = data.frame(from_age = c("adult", NA), to_age = c("adult", NA),
                          from_class = c(NA, "a"), to_class = c(NA, "a"),
                          param_name = c("ctr_aa", "ctr_class_aa")),
      symmetric = FALSE,
      # prime values catch any accidental sum or overwrite
      params = list(ctr_aa = 3, ctr_class_aa = 5)
    ) |>
    inspect_house_state("h1")

  # m1 (adult,a) -> m2 (adult,a): both rules match -> 3 * 5 = 15
  expect_equal(ctr_cell(snapshot, "m1", "m2", key = "id"), 15)
  # m1 (adult,a) -> m4 (adult,b): only age rule matches -> 3
  expect_equal(ctr_cell(snapshot, "m1", "m4", key = "id"), 3)
  # m1 (adult,a) -> m3 (child,b): neither rule matches -> reference (1)
  expect_equal(ctr_cell(snapshot, "m1", "m3", key = "id"), 1)
  # diagonal is 0
  expect_equal(unname(diag(snapshot$contact_risks)), rep(0., 4))
})

test_that("Contact risk definition work as expected", {
  hhmodel <- data.frame(
    house_id = c(rep("h1", 3), rep("h2", 6)),
    status = c(1, 0, 0, 1, 0, 0, 0, 0, 0),
    member = c("student", "other", "other",
               "student", "sibling", "mother", "father", "other", "student")
  ) |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")

  # no contact risk defs -> empty matrix
  contact_risks <- inspect_house_state(hhmodel, "h2")$contact_risks
  expect_true(is.matrix(contact_risks))
  expect_equal(dim(contact_risks), c(0L, 0L))

  # mutually exclusive strata
  strata <- data.frame(
    from_member = c("student", "student", "sibling", "mother",
                    "mother", "father", "other", NA),
    to_member = c("student", "sibling", "sibling", "student",
                  "sibling", "mother", "other", NA),
    param_name = c("child_child", "child_child", "child_child", "mother_child",
                   "mother_child", ".reference", "other_other", "cross_gen")
  )
  # prime numbers as param values to ensure non duplicate pair products
  params <- stats::setNames(c(2, 3, 5, 7, 11, 13, 17, 19), strata$param_name)
  hhmodel <- add_contact_risk(hhmodel, strata = strata, params = params)

  # matrix is n x n with a null diagonal
  snapshot <- inspect_house_state(hhmodel, "h2")
  n <- nrow(snapshot$data_wide)
  expect_equal(dim(snapshot$contact_risks), c(n, n))
  expect_equal(unname(diag(snapshot$contact_risks)), rep(0., n))

  # every off-diagonal cell resolves to the reference (1)
  # or to exactly one param value: no multiplications
  expect_true(all(snapshot$contact_risks %in% c(0, 1, params)))

  # the template parameters added previously survive alongside the risk ones
  expect_true(all(c("alpha", "beta", "gamma") %in% hhmodel$parameters$name))
})

test_that("A contact risk keyed on an observation is rejected", {
  # the risk feeds the transmission model the observation is generated from
  hhmodel <- data.frame(house_id = c(1, 1), status = c(1, 0),
                        datum = c(2, 3), datum_2 = c(2, 6)) |>
    new_housepi_model("final-size")

  hhmodel |>
    add_contact_risk(by = "datum") |>
    add_observation("datum", model = "uniform", min = 0, max = 5) |>
    expect_error("keyed on observation.*datum")

  # a deterministic household observation stands in for the column all the same
  hhmodel |>
    add_observation(datum = min(datum_2@D), level = "house") |>
    add_contact_risk(by = "datum") |>
    expect_error("keyed on observation.*datum")
})

test_that("Contact risks raise an error for incompatible models", {
  data.frame(
    house_id = 1, id = seq(4), status = seq(4),
    datum = c("a", "a", "b", "b")
  ) |>
    new_housepi_model("test") |>
    add_contact_risk(by = "datum") |>
    expect_error("[Nn]ot supported")
})

#################################
## CONTACT RISKS : STRATA PINS ##
#################################

test_that("A contact strata value pins the stratum and keeps its prior", {
  pinned <- pin_ctr_strata(const = c(NA, NA)) |>
    get_parameters("ctr_fm")
  expect_equal(pinned$value, 0.4)
  expect_false(pinned$const)
  expect_equal(pinned$prior, "lognormal")

  # a fixed stratum is a constant parameter: it carries no prior
  hhmodel <- pin_ctr_strata()
  fixed <- get_parameters(hhmodel, "ctr_fm")
  expect_equal(fixed$value, 0.4)
  expect_true(fixed$const)
  expect_true(is.na(fixed$prior))

  # the read table keeps the columns typed
  tidy <- get_contact_risks(hhmodel)
  expect_type(tidy$value, "double")
  expect_type(tidy$const, "logical")
})

test_that("Invalid contact strata pins are rejected", {
  # a fixed stratum needs the value it is fixed at
  pin_ctr_strata(value = rep(NA_real_, 2)) |>
    expect_error("ctr_fm")
  # an all-NA column pins nothing, whichever type it was written with
  unpinned <- pin_ctr_strata(value = c(NA, NA), const = c(NA, NA)) |>
    get_parameters("ctr_fm")
  expect_false(unpinned$const)
  expect_equal(unpinned$prior, "lognormal")
  # the reference is 1 by definition
  pin_ctr_strata(value = c(0.4, 0.9), const = c(NA, NA)) |>
    expect_error("[Rr]eference")
  # a pinned stratum takes no entry in params
  pin_ctr_strata(params = list(ctr_fm = 2)) |>
    expect_error("ctr_fm")
  # a stratum spread over two pairs carries a single pin, even when the data
  # carries no grandparent to keep the second pair
  strata <- data.frame(from_member = c("father", "father", NA),
                       to_member = c("mother", "grandparent", NA),
                       param_name = c("ctr_fm", "ctr_fm", ".reference"),
                       value = c(0.4, 0.9, NA))
  mk_house$model  |>
    add_contact_risk(strata = strata) |>
    expect_error("[Cc]onflicting.*value.*ctr_fm.*rows 1 and 2")
})

test_that("A symmetric pin reaches the mirrored pair", {
  snapshot <- pin_ctr_strata(symmetric = TRUE) |>
    inspect_house_state("h1")

  # the mirrored row inherits both the parameter and its pin
  expect_equal(ctr_cell(snapshot, "father", "mother"), 0.4)
  expect_equal(ctr_cell(snapshot, "mother", "father"), 0.4)
  # every other pair falls through to the reference
  expect_equal(ctr_cell(snapshot, "student", "sibling"), 1)

  # symmetric = FALSE leaves the reverse direction at the reference
  snapshot <- pin_ctr_strata(symmetric = FALSE) |>
    inspect_house_state("h1")
  expect_equal(ctr_cell(snapshot, "father", "mother"), 0.4)
  expect_equal(ctr_cell(snapshot, "mother", "father"), 1)
})
