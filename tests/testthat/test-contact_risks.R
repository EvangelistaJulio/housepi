library(testthat)

# model over two covariates
mk_risk <- local({
  hhmodel <- data.frame(hhid = 1L, status = 0L,
                        age = c(rep("child", 2), rep("adult", 2)),
                        sex = c("M", "F", "M", "M")) |>
    new_housepi_model("final-size")
  list(model = hhmodel)
})

# house holding one of each member type, with template parameters
mk_house <- local({
  model <- data.frame(id = c("m1", "m2", "m3", "m4"),
                      hhid = "h1", status = 0L,
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
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,  ~value,   ~const,
    "father",     "mother",   "ctr_fm",     value[1], const[1],
    NA,           NA,         ".reference", value[2], const[2]
  )
  mk_house$model  |>
    add_contact_risk(strata = strata, ...)
}

#################################
## CONTACT RISKS : PAIR COUNTS ##
#################################

test_that("Contact risks pair counts are computed correctly", {
  tbl <- data.frame(hhid = c(rep("h1", 2L), rep("h2", 2L)),
                    age = rep(c("child", "adult"), 2L)) |>
    contact_risk_logic_table("age")
  count <- \(from, to) tbl$n[tbl$from_age == from & tbl$to_age == to]

  # total: child->adult = 2, adult->child = 2, no same-age pair
  count("child", "adult") |> expect_identical(2L)
  count("adult", "child") |> expect_identical(2L)
  count("child", "child") |> expect_identical(integer(0))
  count("adult", "adult") |> expect_identical(integer(0))

  # no household structure, no pair to count
  contact_risk_logic_table(data.frame(age = c("child", "adult")), "age") |>
    expect_null()
})

test_that("Pair counts drop the member meeting itself", {
  # 2 children in one household imply 2 ordered child->child pairs, not 4
  tbl <- data.frame(hhid = c(rep("h1", 3), rep("h2", 2)),
                    age = c(rep("child", 2), rep("adult", 3))) |>
    contact_risk_logic_table("age")
  count <- \(from, to) tbl$n[tbl$from_age == from & tbl$to_age == to]

  count("child", "child") |> expect_identical(2L)
  count("adult", "adult") |> expect_identical(2L)
  count("child", "adult") |> expect_identical(2L)
  count("adult", "child") |> expect_identical(2L)
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
  tidy |> expect_nrow(4L)

  # only 1 reference
  sum(tidy$param_name == ".reference") |> expect_identical(1L)
  # only 2 unique non-reference params (cross + child->child)
  setdiff(unique(tidy$param_name), ".reference") |> expect_length(2)

  # reference row is adult->adult
  ref_row <- tidy[tidy$param_name == ".reference", ]
  ref_row$from_age |> expect_identical("adult")
  ref_row$to_age |> expect_identical("adult")

  # symmetric: adult->child and child->adult share the same parameter
  subset(tidy, from_age == "adult" & to_age == "child")$param_name |>
    expect_identical(
      subset(tidy, from_age == "child" & to_age == "adult")$param_name
    )

  # symmetric = FALSE gives distinct parameters to the two directions
  tidy <- risks_of(symmetric = FALSE)
  tidy |> expect_nrow(4L)
  setdiff(unique(tidy$param_name), ".reference") |> expect_length(3L)
  identical(
    subset(tidy, from_age == "adult" & to_age == "child")$param_name,
    subset(tidy, from_age == "child" & to_age == "adult")$param_name
  ) |>
    expect_false()
})

test_that("Auto mode reference is opt-in and picks most abundant pair", {
  tidy <- mk_risk$model |>
    add_contact_risk(by = "age", reference = ".auto") |>
    get_contact_risks()

  # automatic .reference on the most abundant pair,
  # tagged in both directions since the contact is symmetric
  ref <- subset(tidy, param_name == ".reference")
  nrow(ref) |> expect_gte(1L)
  ref$n |> expect_all_equal(max(tidy$n))

  # absolute by default: every pair keeps its own parameter
  ctr <- mk_risk$model |>
    add_contact_risk(by = "age") |>
    get_contact_risks()
  subset(ctr, param_name == ".reference") |>
    expect_nrow(0L)
})

test_that("Non-joint produces independent per-covariate tables", {
  # 2 covariates, non-joint -> 2 rule groups
  tidy <- mk_risk$model |>
    add_contact_risk(by = c("age", "sex"), joint = FALSE,
                     reference = c(age = "adult", sex = "M"),
                     prior = "lognormal", sdlog = 1.) |>
    get_contact_risks()

  length(unique(tidy$rule_group)) |> expect_identical(2L)

  # all age pairs are observed: 2 x 2
  sum(tidy$rule_group == 1) |> expect_identical(4L)
  # 1 sex pair is missing: 2 x 2 - 1
  sum(tidy$rule_group == 2) |> expect_identical(3L)

  # each sub-table has its own .reference
  for (grp in unique(tidy$rule_group)) {
    sum(tidy$param_name[tidy$rule_group == grp] == ".reference") |>
      expect_identical(1L)
  }

  # symmetric default: 2 + 2 - 1 = 3 non-reference params
  setdiff(unique(tidy$param_name), ".reference") |>
    expect_length(3L)

  # first sub-table carries age, and leaves sex unset
  first <- subset(tidy, rule_group == 1)
  first$from_age |> expect_no_na()
  first$to_age |> expect_no_na()
  first$from_sex |> expect_all_na()
  first$to_sex |> expect_all_na()

  # second sub-table carries sex, and leaves age unset
  second <- subset(tidy, rule_group == 2)
  second$from_sex |> expect_no_na()
  second$to_sex |> expect_no_na()
  second$from_age |> expect_all_na()
  second$to_age |> expect_all_na()
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
  tidy |> expect_nrow(7L)
  unique(tidy$rule_group) |> expect_length(1L)
  sum(tidy$param_name == ".reference") |> expect_identical(1L)

  # NOTE: symmetric pairs by default:
  # 3 heterogeneous: <child,M><child,F> <adult,M><child,F> <adult,M><child,M>
  # 1 homogeneous: <adult,M><adult,M> (reference)
  setdiff(unique(tidy$param_name), ".reference") |> expect_length(3L)
})

test_that("Absoluteness follows the reference across every definition mode", {
  # the default prior tracks whether a risk carries absolute values
  risk_prior <- function(...) {
    params <- get_parameters(add_contact_risk(mk_risk$model, ...))
    unique(params$prior[!params$const])
  }

  # auto mode: absolute unless a reference is provided
  risk_prior(by = "age") |> expect_identical("uniform")
  risk_prior(by = "age", reference = ".auto") |> expect_identical("lognormal")
  risk_prior(by = "age", reference = "adult") |> expect_identical("lognormal")

  # custom mode: the reference parameter decides
  anchored <- tibble::tribble(
    ~from_age, ~to_age, ~param_name,
    "adult",   "adult", ".reference",
    "adult",   "child", "a_to_c",
    "child",   "child", "c_to_c"
  )
  unanchored <- anchored
  unanchored$param_name[1] <- "a_to_a"

  risk_prior(strata = unanchored) |> expect_identical("uniform")
  risk_prior(strata = anchored) |> expect_identical("lognormal")
})

test_that("Two absolute rule groups are not identifiable", {
  mk_risk$model |>
    add_contact_risk(by = c("age", "sex")) |>
    add_parameter(alpha = 0.1, gamma = 1, hhsize_ref = 4) |>
    evaluate_log_lik() |>
    expect_error("2 absolute rule groups.*non-identifiability")
})

test_that("Every rule group keeps its anchoring in auto", {
  by <- c("age", "sex")
  hhmodel <-  sarscov2_2021_icofs |>
    dplyr::filter(hhid %in% utils::head(unique(hhid), 50)) |>
    new_housepi_model("final-size")

  contact_prior <- function(hhmodel) {
    params <- get_parameters(hhmodel)
    params$prior[grepl("^(age|sex)_", params$name)]
  }

  # no reference makes it scale free and sets uniform prior
  hhmodel |>
    add_contact_risk(by = by) |>
    contact_prior() |>
    expect_all_equal("uniform")

  # the sex rule group left without reference makes the risk absolute
  hhmodel |>
    add_contact_risk(by = by,
                     reference = c(from_age = "adult", to_age = "adult")) |>
    contact_prior() |>
    expect_all_equal("uniform")

  # a reference in every rule group anchors the risk: lognormal prior
  hhmodel <- hhmodel |>
    add_contact_risk(by = by,
                     reference = c(from_age = "adult", to_age = "adult",
                                   from_sex = "F", to_sex = "F"))
  contact_prior(hhmodel) |> expect_all_equal("lognormal")

  # joint = FALSE by default so each covariate owns a rule group
  get_contact_risks(hhmodel)$rule_group |>
    expect_setequal(c(1L, 2L))
})

test_that("Spare prior arguments on a contact risk are reported", {
  mk_risk$model |>
    add_contact_risk(by = "age", sdlog = 1.8) |>
    expect_error("[Ss]pare prior arguments.*contact risk")
})

test_that("A reference naming part of its pair is refused", {
  # one end of the pair
  mk_risk$model |>
    add_contact_risk(by = "age", reference = c(from_age = "adult")) |>
    expect_error("[Rr]eference names from_age but leaves to_age unnamed")

  # part of a joint cell, prefixed or bare
  mk_risk$model |>
    add_contact_risk(by = c("age", "sex"), joint = TRUE,
                     reference = c(from_age = "adult", to_age = "adult")) |>
    expect_error("leaves from_sex and to_sex unnamed")
  mk_risk$model |>
    add_contact_risk(by = c("age", "sex"), joint = TRUE,
                     reference = c(age = "adult")) |>
    expect_error("leaves from_sex and to_sex unnamed")
})

test_that("An automatic contact risk refuses pairs missing a covariate", {
  # the member with no role meets the two others both ways: 4 ordered pairs
  data.frame(hhid = 1, status = 0, member = c("father", "mother", NA)) |>
    new_housepi_model("final-size") |>
    add_contact_risk(by = "member") |>
    expect_error("leaves 4 ordered pairs unmatched")

  # no member carries the covariate: the table holds no rule
  data.frame(hhid = 1, status = 0, member = rep(NA_character_, 2)) |>
    new_housepi_model("final-size") |>
    add_contact_risk(by = "member") |>
    expect_error("leaves 2 ordered pairs unmatched")
})

test_that("Households of singles give an automatic risk no rule", {
  data.frame(hhid = c("h1", "h2"), status = 0, member = c("a", "b")) |>
    new_housepi_model("final-size") |>
    add_contact_risk(by = "member") |>
    get_contact_risks() |>
    expect_nrow(0L)
})

test_that("An invalid contact reference is rejected", {
  mk_risk$model |>
    add_contact_risk(by = "age",
                     reference = c(from_age = "elderly", to_age = "adult")) |>
    expect_error("age = .elderly. <-> age = .adult. is never observed")
  mk_risk$model |>
    add_contact_risk(by = "age", reference = c(age_from = "adult")) |>
    expect_error("[Uu]nrecognized choice for reference: .age_from.")
  # two covariates leave an unnamed reference ambiguous
  mk_risk$model |>
    add_contact_risk(by = c("age", "sex"), reference = c("adult", "M")) |>
    expect_error("must be named using covariates")
  mk_risk$model |>
    add_contact_risk(by = c("age", "sex"), reference = "adult") |>
    expect_error("must be named using covariates")
})

test_that("An unnamed reference names the single covariate", {
  tidy <- mk_risk$model |>
    add_contact_risk(by = "age", reference = "adult") |>
    get_contact_risks()

  ref <- subset(tidy, param_name == ".reference")
  c(ref$from_age, ref$to_age) |> expect_identical(c("adult", "adult"))
})

test_that("Contact risks expand on their own covariates only", {
  hhdata <- sarscov2_2021_icofs |>
    dplyr::select(-start_date) |>
    dplyr::filter(hhid %in% utils::head(unique(hhid), 50))

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
  risks |> expect_nrow(0L)
  risks |> expect_named(c("param_name", "n", "rule_group"))
})

########################################
## CONTACT RISKS : CUSTOM DEFINITIONS ##
########################################

test_that("Contact risk custom mode default behaves correctly", {
  # logic table
  strata <- tibble::tribble(
    ~from_age, ~to_age, ~param_name,
    "adult",   "adult", ".reference",
    "adult",   "child", "a2c",
    "child",   "adult", "c2a",
    "child",   "child", "c2c"
  )
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
  tidy |> expect_nrow(4L)
  tidy$param_name |> expect_setequal(strata$param_name)

  # parameters generated for non-reference
  get_parameters(hhmodel) |>
    getElement("name") |>
    expect_setequal(c("a2c", "c2a", "c2c"))

  # default parameter definitions
  default <- get_parameters(hhmodel, "c2c")
  default$prior |> expect_identical("lognormal")
  default$args[[1]] |> expect_equal(c(sdlog = 1.8, meanlog = 0))
  # overrides as variable
  variable <- get_parameters(hhmodel, "a2c")
  variable$prior |> expect_identical("uniform")
  variable$min |> expect_equal(1)
  # overrides as constant
  constant <- get_parameters(hhmodel, "c2a")
  constant$value |> expect_equal(2.34)
  constant$const |> expect_true()
})

test_that("A contact strata row without a param_name is refused", {
  strata <- tibble::tribble(~from_age, ~to_age, ~param_name,
                            NA,        NA,      ".reference",
                            "child",   "child", "")

  mk_risk$model |>
    add_contact_risk(strata = strata) |>
    expect_error("[Ss]trata row without param_name: 2")
})

test_that("Contact strata leaving pairs unmatched are refused", {
  # the two child -> child pairs match no rule
  strata <- tibble::tribble(~from_age, ~to_age, ~param_name,
                            "adult",   "adult", "a_to_a",
                            "adult",   "child", "a_to_c")
  mk_risk$model |>
    add_contact_risk(strata = strata) |>
    expect_error("[Cc]ontact risk leaves 2 ordered pairs unmatched")
})

test_that("A reference has no meaning against a custom contact table", {
  strata <- tibble::tribble(~from_age, ~to_age, ~param_name,
                            "adult",   "adult", "a_to_a",
                            "adult",   "child", "a_to_c",
                            "child",   "adult", "c_to_a",
                            "child",   "child", "c_to_c")

  mk_risk$model |>
    add_contact_risk(strata = strata, reference = ".auto") |>
    expect_error("[Nn]ot used with a custom")
})

test_that("Custom contact strata carry at most one reference per rule group", {
  # a symmetric pair reaches one cell through either direction
  strata <- tibble::tribble(~from_age, ~to_age, ~param_name,
                            "adult",   "child", ".reference",
                            "child",   "adult", ".reference",
                            "child",   "child", "c_to_c",
                            "adult",   "adult", "a_to_a")
  mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = TRUE) |>
    get_contact_risks() |>
    getElement("param_name") |>
    expect_setequal(c(".reference", "c_to_c", "a_to_a"))
  # without symmetry the two directions are two competing scales
  mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    expect_error("[Rr]ows.*one rule group")
})

test_that("Custom mode fills the symmetric counterpart of a pair", {
  strata <- tibble::tribble(~from_age, ~to_age, ~param_name,
                            "adult",   "adult", ".reference",
                            "adult",   "child", "a_to_c",
                            "child",   "child", "c_to_c")

  hhmodel <- mk_risk$model |>
    add_contact_risk(strata = strata, params = list(a_to_c = 2))

  # symmetric fill adds child->adult with same param_name as adult->child
  filled <- get_contact_risks(hhmodel)
  filled |> expect_nrow(4L)
  filled$param_name[filled$from_age == "child" & filled$to_age == "adult"] |>
    expect_identical("a_to_c")

  # both directions share the override value
  get_parameters(hhmodel) |>
    getElement("name") |>
    expect_setequal(c("a_to_c", "c_to_c"))
  a_to_c <- get_parameters(hhmodel, "a_to_c")
  a_to_c$value |> expect_equal(2)
  a_to_c$const |> expect_true()

  # symmetric = FALSE does not fill gaps: child->adult pairs match no rule
  mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    expect_error("leaves 4 ordered pairs unmatched")
})

test_that("An unobserved contact level is refused", {
  strata <- tibble::tribble(~from_member, ~to_member, ~param_name,
                            "mother",     "uncle",    "rho",
                            NA,           NA,         ".reference")
  mk_house$model |>
    add_contact_risk(strata = strata, symmetric = FALSE,
                     params = list(rho = 2)) |>
    expect_error("no ordered pair resolves.*uncle.*rho")
})

test_that("A catch-all contact rule counts every pair it reaches", {
  # 4 members: 12 ordered pairs, one of them father -> mother
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "father",     "mother",   "rho",
    NA,           NA,         ".reference"
  )
  hhmodel <- mk_house$model |>
    add_contact_risk(strata = strata, symmetric = FALSE,
                     params = list(rho = 2))
  risks <- hhmodel |>
    get_contact_risks()

  ctr <- inspect_house_state(hhmodel, "h1")$contact_risks
  risks$n[risks$param_name == "rho"] |> expect_identical(sum(ctr == 2))
  risks$n[risks$param_name == ".reference"] |> expect_identical(sum(ctr == 1))
})

test_that("A contact table setting no covariate value is refused", {
  # a single risk over every pair duplicates the baseline
  hhmodel <- data.frame(hhid = c(1, 1, 2, 2, 2), status = c(1, 0, 1, 1, 0),
                        age = c("child", "adult", "child", "adult", "adult")) |>
    new_housepi_model("final-size")
  hhmodel |>
    add_contact_risk(strata = data.frame(param_name = "ctr")) |>
    expect_error("[Ss]trata table sets no covariate value")
  hhmodel |>
    add_contact_risk(strata = data.frame(from_age = NA, to_age = NA,
                                         param_name = "ctr")) |>
    expect_error("[Ss]trata table sets no covariate value")
})

test_that("Households that never form a pair leave every contact rule dead", {
  # two singleton households: no member ever meets another
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "father",     "mother",   "rho",
    NA,           NA,         ".reference"
  )
  data.frame(hhid = c("h1", "h2"), status = 0,
             member = c("father", "mother")) |>
    new_housepi_model("final-size") |>
    add_contact_risk(strata = strata, symmetric = FALSE,
                     params = list(rho = 2)) |>
    expect_error("2 rules no ordered pair resolves to")
})

test_that("An unknown contact level is refused as written, not mirrored", {
  # the rule is checked before the symmetric fill: the message names the
  # direction the table carries, not its mirror
  strata <- tibble::tribble(
    ~from_member, ~to_member,    ~param_name,
    "mother",     "student",     "rho",
    "mother",     "grandparent", "rho",
    NA,           NA,            ".reference"
  )
  mk_house$model |>
    add_contact_risk(strata = strata, params = list(rho = 2)) |>
    expect_error(paste("1 rule no ordered pair resolves to:",
                       "from_member = .mother. and to_member = .grandparent."))
})

test_that("A params entry names a parameter of the table", {
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "father",     "mother",   "ctr_fm"
  )
  mk_house$model |>
    add_contact_risk(strata = strata, params = list(unknown = 2.4)) |>
    expect_error("[Uu]nknown parameter")
})

test_that("An automatic definition takes no params list", {
  # `by` derives its own parameters, so the list has nothing to override,
  # and with no reference pair the risk also takes the scale of `beta`
  mk_house$model |>
    add_contact_risk("member", params = list(ctr_fm = 2.4)) |>
    expect_warning("[Pp]arams.*not used with") |>
    expect_warning("[Bb]aseline.*constant reference")
})

#############################################
###  CONTACT RISKS : RULE GROUPING LOGIC  ###
#############################################

test_that("A rule fitting two rule groups is refused as ambiguous", {
  # rho_3 keys on from_age alone: it fits {from_age, from_sex} and
  # {from_age, to_age} equally
  strata <- tibble::tribble(
    ~from_age, ~from_sex, ~to_age, ~to_sex, ~param_name,
    "child",   "M",       NA,      NA,      "rho",
    "child",   NA,        "child", NA,      "rho_2",
    "child",   NA,        NA,      NA,      "rho_3"
  )
  mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    expect_error("[Aa]mbiguous rule 3 .*rho_3.*fitting 2")
})

test_that("A contact column without its partner is refused", {
  strata <- tibble::tribble(~from_age, ~param_name,
                            "child",   "rho",
                            NA,        ".reference")
  mk_risk$model |>
    add_contact_risk(strata = strata) |>
    expect_error("[Cc]ontact strata column from_age has no partner")

  # "from children to anyone" names its catch-all partner
  strata$to_age <- NA
  mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    get_contact_risks() |>
    getElement("param_name") |>
    expect_setequal(c("rho", ".reference"))
})

test_that("Custom mode absorbs all-NA default into single group", {
  # 4 rules in one group, last one (all-NA) is the catch-all
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "student",    "student",  "ctr_child_child",
    "mother",     "student",  "ctr_mother_child",
    "father",     "mother",   ".reference",
    NA,           NA,         "ctr_cross_gen"
  )
  data.frame(hhid = 1, status = 0,
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
  # a half-keyed rule stays in the group of the pattern containing it
  strata <- tibble::tribble(
    ~from_age, ~from_sex, ~to_age, ~to_sex, ~param_name,
    "child",   NA,        "child", NA,      "age_child_to_age_child",
    "child",   NA,        "adult", NA,      "age_child_to_age_adult",
    "adult",   NA,        NA,      NA,      "age_adult_to_any",
    NA,        "M",       NA,      "M",     ".reference",
    NA,        "M",       NA,      "F",     "sex_M_to_sex_F",
    NA,        "F",       NA,      "M",     "sex_F_to_sex_M"
  )
  mk_risk$model |>
    add_contact_risk(strata = strata, symmetric = FALSE) |>
    get_contact_risks() |>
    getElement("rule_group") |>
    unique() |>
    expect_length(2L)
})

####################################
###  CONTACT RISKS : END-TO-END  ###
####################################

test_that("Most-specific rule wins over reference (single covariate)", {
  # one specific pair gets a non-reference param; everything else falls
  # through to .reference (specificity 0, never beaten)
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "father",     "mother",   "ctr_fm",
    NA,           NA,         ".reference"
  )
  snapshot <- mk_house$model |>
    add_contact_risk(strata = strata, symmetric = FALSE,
                     params = list(ctr_fm = 2)) |>
    inspect_house_state("h1")

  # the explicit pair wins
  ctr_cell(snapshot, "father", "mother") |> expect_equal(2)
  # symmetric = FALSE: reverse falls back to reference
  ctr_cell(snapshot, "mother", "father") |> expect_equal(1)
  # unrelated pair falls back to reference
  ctr_cell(snapshot, "student", "sibling") |> expect_equal(1)
  # diagonal stays 0 (no self-contact)
  unname(diag(snapshot$contact_risks)) |> expect_equal(rep(0., 4))
})

test_that("Rules tied at equal specificity are rejected as ambiguous", {
  # student -> mother is reached by two half-keyed rules, one specificity each,
  # and the fully-keyed pair makes {from, to} the only maximal pattern: the
  # rules share a group, so the engine has to compare the tied two
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "student",    NA,         "rho",
    NA,           "mother",   "rho_2",
    "father",     "mother",   "rho_3",
    NA,           NA,         ".reference"
  )
  mk_house$model |>
    add_contact_risk(strata = strata, symmetric = FALSE,
                     params = list(rho = 2, rho_2 = 3, rho_3 = 4)) |>
    inspect_house_state("h1") |>
    expect_error("[Aa]mbiguous contact risk.*tied rules")
})

test_that("A catch-all rule survives a missing cell in its covariate column", {
  # the catch-all row keys on nothing, so it must reach every pair, including
  # the ones whose covariate cell carries a value
  # the table tags no .reference, so `beta` is collapsed onto it
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "father",     "mother",   "ctr_fm",
    NA,           NA,         "ctr_any"
  )
  snapshot <- mk_house$model |>
    add_data_wide(transform(mk_house$model$data_wide,
                            member = c("father", "mother", "student", NA))) |>
    add_contact_risk(strata = strata, symmetric = FALSE,
                     params = list(ctr_fm = 2, ctr_any = 3)) |>
    inspect_house_state("h1") |>
    suppressWarnings() # suppress baseline collapse warning

  # the specific pair still wins on specificity
  ctr_cell(snapshot, "m1", "m2", key = "id") |> expect_equal(2)
  # the catch-all reaches a pair of members that both carry a role
  ctr_cell(snapshot, "m3", "m1", key = "id") |> expect_equal(3)
  # and the member with no role reads it too
  ctr_cell(snapshot, "m4", "m1", key = "id") |> expect_equal(3)
  ctr_cell(snapshot, "m1", "m4", key = "id") |> expect_equal(3)
})

test_that("Distinct rule_groups multiply at the cell", {
  # two disjoint maximal patterns ({from_age,to_age}, {from_class,to_class})
  # make two rule groups whose winners multiply into one cell
  # a half-keyed rule stays in the group of the pattern containing it
  strata <- tibble::tribble(
    ~from_age, ~to_age, ~from_class, ~to_class, ~param_name,
    "adult",   "adult", NA,          NA,        "ctr_aa",
    "adult",   "child", NA,          NA,        ".reference",
    "child",   NA,      NA,          NA,        "ctr_child",
    NA,        NA,      "a",         "a",       "ctr_class_aa",
    NA,        NA,      "a",         "b",       ".reference",
    NA,        NA,      "b",         NA,        "ctr_class_b"
  )
  snapshot <- data.frame(
    id = c("m1", "m2", "m3", "m4"), hhid = "h1", status = 0L,
    age = c("adult", "adult", "child", "adult"),
    class = c("a", "a", "b", "b")
  ) |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    add_contact_risk(
      strata = strata, symmetric = FALSE,
      # prime values catch any accidental sum or overwrite
      params = list(ctr_aa = 3, ctr_class_aa = 5, ctr_child = 7,
                    ctr_class_b = 11)
    ) |>
    inspect_house_state("h1")

  # m1 (adult,a) -> m2 (adult,a): both rules match -> 3 * 5 = 15
  ctr_cell(snapshot, "m1", "m2", key = "id") |> expect_equal(15)
  # m1 (adult,a) -> m4 (adult,b): age rule and class reference -> 3
  ctr_cell(snapshot, "m1", "m4", key = "id") |> expect_equal(3)
  # m1 (adult,a) -> m3 (child,b): both references -> 1
  ctr_cell(snapshot, "m1", "m3", key = "id") |> expect_equal(1)
  # m3 (child,b) -> m1 (adult,a): both half-keyed rules -> 7 * 11 = 77
  ctr_cell(snapshot, "m3", "m1", key = "id") |> expect_equal(77)
  # diagonal is 0
  unname(diag(snapshot$contact_risks)) |> expect_equal(rep(0., 4))
})

test_that("Contact risk definition work as expected", {
  hhmodel <- data.frame(
    hhid = c(rep("h1", 3), rep("h2", 6)),
    status = c(1, 0, 0, 1, 0, 0, 0, 0, 0),
    member = c("student", "other", "other",
               "student", "sibling", "mother", "father", "other", "student")
  ) |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")

  # no contact risk defs -> empty matrix
  contact_risks <- inspect_house_state(hhmodel, "h2")$contact_risks
  is.matrix(contact_risks) |> expect_true()
  dim(contact_risks) |> expect_identical(c(0L, 0L))

  # mutually exclusive strata
  strata <- tibble::tribble(
    ~from_member, ~to_member, ~param_name,
    "student",    "student",  "child_child",
    "student",    "sibling",  "child_child",
    "mother",     "student",  "mother_child",
    "mother",     "sibling",  "mother_child",
    "father",     "mother",   ".reference",
    "other",      "other",    "other_other",
    NA,           NA,         "cross_gen"
  )
  # prime numbers as param values to ensure non duplicate pair products
  params <- stats::setNames(c(2, 3, 5, 7, 11, 13, 17), strata$param_name)
  params <- params[names(params) != ".reference"]
  hhmodel <- add_contact_risk(hhmodel, strata = strata, params = params)

  # matrix is n x n with a null diagonal
  snapshot <- inspect_house_state(hhmodel, "h2")
  n <- nrow(snapshot$data_wide)
  dim(snapshot$contact_risks) |> expect_equal(c(n, n))
  unname(diag(snapshot$contact_risks)) |> expect_equal(rep(0., n))

  # every off-diagonal cell resolves to the reference (1)
  # or to exactly one param value: no multiplications
  snapshot$contact_risks |> expect_in(c(0, 1, params))

  # the template parameters added previously survive alongside the risk ones
  hhmodel$parameters$name |> expect_contains(c("alpha", "beta", "gamma"))
})

test_that("A contact risk keyed on an observation is rejected", {
  # the risk feeds the transmission model the observation is generated from
  hhmodel <- data.frame(hhid = c(1, 1), status = c(1, 0),
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
  data.frame(hhid = 1, id = seq(4), status = seq(4),
             datum = c("a", "a", "b", "b")) |>
    new_housepi_model("test") |>
    add_contact_risk(by = "datum") |>
    expect_error("[Nn]ot supported")
})

#####################################
###  CONTACT RISKS : STRATA PINS  ###
#####################################

test_that("A contact strata value pins the stratum and keeps its prior", {
  pinned <- pin_ctr_strata(const = c(NA, NA)) |>
    get_parameters("ctr_fm")
  pinned$value |> expect_equal(0.4)
  pinned$const |> expect_false()
  pinned$prior |> expect_identical("lognormal")

  # a fixed stratum is a constant parameter: it carries no prior
  hhmodel <- pin_ctr_strata()
  fixed <- get_parameters(hhmodel, "ctr_fm")
  fixed$value |> expect_equal(0.4)
  fixed$const |> expect_true()
  fixed$prior |> expect_all_na()

  # the read table keeps the columns typed
  tidy <- get_contact_risks(hhmodel)
  tidy$value |> expect_type("double")
  tidy$const |> expect_type("logical")
})

test_that("Invalid contact strata pins are rejected", {
  # a fixed stratum needs the value it is fixed at
  pin_ctr_strata(value = rep(NA_real_, 2)) |>
    expect_error("ctr_fm")
  # an all-NA column pins nothing, whichever type it was written with
  unpinned <- pin_ctr_strata(value = c(NA, NA), const = c(NA, NA)) |>
    get_parameters("ctr_fm")
  unpinned$const |> expect_false()
  unpinned$prior |> expect_identical("lognormal")
  # the reference is 1 by definition
  pin_ctr_strata(value = c(0.4, 0.9), const = c(NA, NA)) |>
    expect_error("[Rr]eference")
  # a pinned stratum takes no entry in params
  pin_ctr_strata(params = list(ctr_fm = 2)) |>
    expect_error("ctr_fm")
  # a stratum spread over two pairs carries a single pin, even when the data
  # carries no grandparent to keep the second pair
  strata <- tibble::tribble(
    ~from_member, ~to_member,    ~param_name,  ~value,
    "father",     "mother",      "ctr_fm",     0.4,
    "father",     "grandparent", "ctr_fm",     0.9,
    NA,           NA,            ".reference", NA
  )
  mk_house$model  |>
    add_contact_risk(strata = strata) |>
    expect_error("[Cc]onflicting.*value.*ctr_fm.*rows 1 and 2")
})

test_that("A symmetric pin reaches the mirrored pair", {
  snapshot <- pin_ctr_strata(symmetric = TRUE) |>
    inspect_house_state("h1")

  # the mirrored row inherits both the parameter and its pin
  ctr_cell(snapshot, "father", "mother") |> expect_equal(0.4)
  ctr_cell(snapshot, "mother", "father") |> expect_equal(0.4)
  # every other pair falls through to the reference
  ctr_cell(snapshot, "student", "sibling") |> expect_equal(1)

  # symmetric = FALSE leaves the reverse direction at the reference
  snapshot <- pin_ctr_strata(symmetric = FALSE) |>
    inspect_house_state("h1")
  ctr_cell(snapshot, "father", "mother") |> expect_equal(0.4)
  ctr_cell(snapshot, "mother", "father") |> expect_equal(1)
})
