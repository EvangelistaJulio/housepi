library(testthat)

####  Helpers  ####

# Model over two covariates. Covariate columns stay character on purpose:
# levels follow order of first appearance, which fixes the reference stratum.
mk_rate <- local({
  hh_data <- data.frame(hh_id = 1, status = 0,
                        age = c("child", "child", "adult", "adult"),
                        sex = c("M", "F", "M", "M"),
                        isolation = c(FALSE, FALSE, TRUE, TRUE))
  strata <- tibble::tribble(~isolation, ~param_name,
                            FALSE,      ".reference",
                            TRUE,       "inf_isol")
  list(data = hh_data, strata = strata,
       model = new_housepi_model(hh_data, "final-size"))
})

# five households over a "member" covariate, "student" being the only role
# present in every household
mk_flu <- local({
  hh_data <- data.frame(
    hh_id = c(rep(1, 2), rep(2, 2), rep(3, 3), rep(4, 2), rep(5, 3)),
    status = c(1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    member = c("student", "other",
               "student", "sibling",
               "student", "mother", "father",
               "student", "other",
               "student", "mother", "father")
  )
  hh_model <- hh_data |>
    new_housepi_model("final-size")
  pinned <- hh_model |>
    add_parameter(alpha = 0.01, beta = 0.1, gamma = 1, hh_size_ref = 4)
  list(data = hh_data, model = hh_model, pinned = pinned)
})

# real dataset for validation: add_group_rate() copies before mutating
mk_icofs <- local({
  list(model = new_housepi_model(sarscov2_2021_icofs, "final-size"))
})

# two observed infected statuses, plus a data column duplicating them on the
# infected members only: the covariate a status stratification replaces
mk_status <- local({
  status <- c("escaped", "symptomatic", "asymptomatic", "symptomatic",
              "escaped", "asymptomatic", "symptomatic", "escaped",
              "asymptomatic", "symptomatic", "asymptomatic", "escaped")
  hh_data <- data.frame(hh_id = rep(1:4, each = 3), status = status,
                        datum = ifelse(status == "escaped", NA, status))
  list(data = hh_data, model = new_housepi_model(hh_data, "final-size"))
})

# numeric covariates a rate can stratify on, and an observation can read
mk_observed <- local({
  hh_model <- data.frame(hh_id = c(1, 1), status = c(1, 0),
                         datum = c(2, 3), datum_2 = c(2, 6)) |>
    new_housepi_model("final-size")
  list(model = hh_model)
})

# a simulation table giving every free parameter of the model a value
pin_all_params <- function(hh_model, value = 0.5) {
  defs <- get_parameters(hh_model)
  names <- defs$name[is.na(defs$const)]
  as.data.frame(stats::setNames(as.list(rep(value, length(names))), names))
}

# custom strata over `mk_flu`, pinning two strata
# a catch-all stratum covers the remaining members
pin_strata <- function(const = c(NA, 1.2, 0.5)) {
  tibble::tribble(
    ~member,    ~param_name,  ~const,
    "student",  ".reference", const[1],
    "mother",   "sus_mum",    const[2],
    "father",   "sus_dad",    const[3],
    NA,         "sus_rest",   const[4]
  )
}

####  Logic table  ####

test_that("Group rates logic table computes combinations and counts", {
  tbl <- group_rate_logic_table(mk_rate$data, c("age", "sex"))

  # counts sum to the number of rows, one count per observed combination
  sum(tbl$n) |> expect_identical(nrow(mk_rate$data))
  sort(tbl$n) |> expect_identical(c(1L, 1L, 2L))
  colnames(tbl) |> expect_contains(c("age", "sex", "n"))
})

test_that("Associate param names decorates table and marks reference", {
  tbl <- group_rate_logic_table(mk_rate$data, c("age", "sex"))
  res <- associate_group_param_names(tbl, prefix = "sus", reference = ".auto")

  # exactly one .reference row (largest group)
  sum(res$param_name == ".reference") |> expect_identical(1L)
  # other param_names start with the provided prefix
  others <- res$param_name[res$param_name != ".reference"]
  others |> expect_match("^sus_")

  # absolute by default: every stratum keeps its own parameter
  abs_res <- associate_group_param_names(tbl, prefix = "sus", reference = NULL)
  abs_res$param_name |> expect_disjoint(".reference")
  abs_res$param_name |> expect_match("^sus_")
})

test_that("An empty logic table passes through grouping and ordering", {
  tbl <- rate_rule_groups_(
    data.frame(age = character(), param_name = character(), n = integer())
  )
  tbl |> expect_nrow(0L)
  colnames(tbl) |> expect_contains("rule_group")
  reorder_rates_logic_table(tbl) |> expect_identical(tbl)
})

test_that("Custom strata logic table is validated", {
  # must be a data.frame
  mk_rate$model |>
    add_group_rate("susceptibility", strata = list(a = 1)) |>
    expect_error("data\\.frame")
  # must have param_name column
  mk_rate$model |>
    add_group_rate("susceptibility", strata = data.frame(x = 1)) |>
    expect_error("param_name")
  # must have at least one row
  mk_rate$model |>
    add_group_rate("susceptibility",
                   strata = data.frame(param_name = character(0))) |>
    expect_error("at least one row")
  # reserved covariate column names rejected
  mk_rate$model |>
    add_group_rate("susceptibility",
                   strata = data.frame(param_name = "a", id = 1)) |>
    expect_error("[Rr]eserved")
})

####  Standard calls  ####

test_that("Group rates integrate with data", {
  hh_model <- mk_rate$model |>
    add_group_rate("community", by = c("age", "sex"), joint = TRUE)

  tidy <- get_group_rates(hh_model)
  tidy$group_rate |> expect_all_equal("community")
  # joint expansion: one rule group over both covariates
  unique(tidy$rule_group) |> expect_length(1L)
  colnames(tidy) |> expect_contains(c("age", "sex"))
  # absolute by default: no reference stratum
  sum(tidy$param_name == ".reference") |> expect_identical(0L)

  # mcmc configuration is the absolute default
  params <- get_parameters(hh_model)
  params <- params[is.na(params$const), ]
  params$min |> expect_all_equal(0)
  params$max |> expect_all_equal(1)
  params$proposal |> expect_all_na()
})

test_that("Reference turns rates relative with its own mcmc configuration", {
  hh_model <- mk_rate$model |>
    add_group_rate("community", by = c("age", "sex"), joint = TRUE,
                   reference = ".auto")

  tidy <- get_group_rates(hh_model)
  sum(tidy$param_name == ".reference") |> expect_identical(1L)
  # the reference is the largest stratum
  tidy$n[tidy$param_name == ".reference"] |> expect_identical(max(tidy$n))

  # mcmc configuration is default (Cauchemez et al., 2009)
  params <- get_parameters(hh_model)
  params$args[[1]] |> expect_equal(c(meanlog = 0, sdlog = 3))
  params$min |> expect_all_equal(0)
  params$max |> expect_all_equal(NA_real_)
  params$proposal |> expect_all_na()
})

test_that("Group rates hold without data and re-derive on add_dataset", {
  # rates ingest data: a data-less model holds them, no rate params yet
  hh_model <- new_housepi_model(model = "final-size") |>
    add_parameter(import = "template") |>
    # relative: the template's baseline stays free alongside the rate
    add_group_rate("susceptibility", by = "member", reference = ".auto")
  get_parameters(hh_model)$name |> expect_no_match("^susceptibility_member")

  # attaching data materialises the rate parameters and runs
  hh_model_2 <- add_dataset(hh_model, mk_flu$data)
  full_params <- get_parameters(hh_model_2)$name
  full_params |> expect_match("^susceptibility_member", all = FALSE)
  evaluate_log_lik(hh_model_2) |> expect_all_finite()

  # replacing data re-derives: a single-level subset drops the stale rate params
  hh_model_3 <- add_dataset(hh_model_2,
                            mk_flu$data[mk_flu$data$member == "student", ])
  get_parameters(hh_model_3)$name |> expect_no_match("^susceptibility_member")
  evaluate_log_lik(hh_model_3) |> expect_all_finite()

  # re-attaching the full data restores the exact original parameter set
  get_parameters(add_dataset(hh_model_3, mk_flu$data))$name |>
    expect_setequal(full_params)
})

test_that("Group rates expand on their own covariates only", {
  hh_data <- sarscov2_2021_icofs |>
    dplyr::filter(hh_id %in% utils::head(unique(hh_id), 50))

  # a data column the rate does not read leaves the expansion untouched
  rates <- hh_data |>
    new_housepi_model("case-ascertained") |>
    add_group_rate("susceptibility", by = "age") |>
    get_group_rates()
  rates$group_rate |> expect_setequal("susceptibility")
  rates$param_name |> expect_contains("susceptibility_age_child")

  # without data the read surface warns and keeps the table shape
  expect_warning(
    rates <- new_housepi_model(model = "final-size") |>
      add_group_rate("susceptibility", by = "member") |>
      get_group_rates(),
    regexp = "without data"
  )
  rates |> expect_nrow(0L)
  rates |> expect_named(c("group_rate", "param_name", "n", "rule_group"))
})

test_that("Group rates sentinels work correctly", {
  # empty rate definition
  mk_icofs$model |>
    add_group_rate("susceptibility", by = NA) |>
    expect_error("[Ii]nvalid rate definition")

  # wrong type for covariate column choice
  mk_icofs$model |>
    add_group_rate("susceptibility", by = 4) |>
    expect_error("[Aa]rgument .by. must be of type character")

  # not implemented group
  mk_icofs$model |>
    add_group_rate("unknown", by = "sex") |>
    expect_error("[Uu]nused group rate.*unknown")

  # incomplete custom prior definition
  mk_icofs$model |>
    add_group_rate("susceptibility", by = "sex", max = 3) |>
    expect_error("[Ss]pare prior arguments")

  mk_icofs$model |>
    add_group_rate("susceptibility", by = "test_date") |>
    expect_error("[Mm]ust be a unidimensional type")
})

test_that("A rate axis of many levels warns once per level count", {
  hh_data <- data.frame(hh_id = "h1", status = 0L,
                        datum = c(letters[1:11], "a"))
  hh_model <- new_housepi_model(hh_data, "final-size")
  expect_warning(
    hh_model <- add_group_rate(hh_model, "community", by = "datum"),
    "[Ss]tratifies on 11 levels.*datum", class = "housepi_warning_levels"
  ) |>
    expect_s3_class("housepi_warning")

  # same level count, a copy included -> silent
  hh_model |>
    add_parameter(beta = 0.1) |>
    expect_no_warning(class = "housepi_warning_levels")
  hh_model$copy() |>
    add_parameter(beta = 0.1) |>
    expect_no_warning(class = "housepi_warning_levels")

  # new level count -> warned again
  hh_data$datum[12] <- "l"
  hh_model |>
    add_dataset(hh_data) |>
    expect_warning("[Ss]tratifies on 12 levels",
                   class = "housepi_warning_levels")
})

test_that("An automatic reference names an observed stratum of its axes", {
  mk_rate$model |>
    add_group_rate("susceptibility", by = "age",
                   reference = list(age = "uncle")) |>
    expect_error("[Ii]ndividual.*age = .uncle.*observed")
  mk_rate$model |>
    add_group_rate("susceptibility", by = c("age", "sex"),
                   reference = list(age = "adult")) |>
    expect_error("by. and .reference. should have equal dimensions")
  mk_rate$model |>
    add_group_rate("susceptibility", by = "age",
                   reference = list(sex = "M")) |>
    expect_error("[Uu]nknown covariate.*reference.*sex.*[Nn]ame each entry")
})

test_that("Group rates reference support behaves as expected", {
  # an untagged custom table is absolute
  strata <- tibble::tribble(~sex, ~param_name,
                            "F",  "ladies",
                            "M",  "gentlemen")
  mk_icofs$model |>
    add_group_rate("susceptibility", strata = strata) |>
    get_group_rates() |>
    getElement("param_name") |>
    expect_identical(c("ladies", "gentlemen"))
  # reference has no meaning against a custom table
  mk_icofs$model |>
    add_group_rate("susceptibility", reference = ".auto", strata = strata) |>
    expect_error("[Nn]ot used with a custom")

  # a tagged row makes the same table relative
  strata$param_name[1] <- ".reference"
  mk_icofs$model |>
    add_group_rate("susceptibility", strata = strata) |>
    get_group_rates() |>
    subset(param_name == ".reference") |>
    expect_nrow(1L)
})

test_that("Custom strata carry at most one reference per rule group", {
  # two references within the same rule group are two competing scales
  strata <- tibble::tribble(~sex, ~param_name,
                            "F",  ".reference",
                            "M",  ".reference")
  mk_icofs$model |>
    add_group_rate("susceptibility", strata = strata) |>
    expect_error("[Rr]ows.*one rule group")

  # distinct rule groups each carry their own reference
  strata <- tibble::tribble(
    ~sex, ~age,     ~param_name,
    "F",  NA,       ".reference",
    "M",  NA,       "gentlemen",
    NA,   "child",  ".reference",
    NA,   "adult",  "grown_ups"
  )
  mk_icofs$model |>
    add_group_rate("susceptibility", strata = strata) |>
    get_group_rates() |>
    subset(param_name == ".reference") |>
    getElement("rule_group") |>
    expect_setequal(1:2)
})

test_that("A group rate cannot combine 'by' with a custom 'strata' table", {
  strata <- tibble::tribble(~sex, ~param_name,
                            "F",  "ladies")
  mk_icofs$model |>
    add_group_rate("susceptibility", by = "sex", strata = strata) |>
    expect_error("[Ii]mpossible.*at the same time")
})

test_that("A group rate the model does not read is refused", {
  mk_rate$model |>
    add_group_rate(infectivty = list(by = "age")) |>
    expect_error("[Uu]nused group rate.*[Dd]id you mean .infectivity")
  # no close name points to the model summary
  mk_rate$model |>
    add_group_rate(unrelated = list(by = "age")) |>
    expect_error("[Uu]nused group rate.*[Pp]rint the model")
})

test_that("A rate prior outside the naming convention is refused", {
  mk_rate$model |>
    add_group_rate("susceptibility", by = "age", prior = "uniform|x") |>
    expect_error("[Pp]rior .uniform.x. contained invalid")
  mk_rate$model |>
    add_group_rate("susceptibility", by = "age", prior = "uniform",
                   proposal = "rw|x") |>
    expect_error("[Pp]roposal .rw.x. contained invalid")
})

test_that("Parameter options work correctly", {
  # the collapsed baseline is not one of the rate parameters
  sampled_params <- function(hh_model) {
    params <- get_parameters(hh_model)
    params[is.na(params$const), ]
  }

  params <- mk_rate$model |>
    add_group_rate("susceptibility", by = c("age", "sex"),
                   prior = "uniform", min = 0, max = 10, joint = TRUE) |>
    sampled_params()
  params$prior |> expect_all_equal("uniform")
  params$min |> expect_all_equal(0)
  params$max |> expect_all_equal(10)
  # a sampled rate parameter declares no proposal
  params$proposal |> expect_all_na()

  params <- mk_rate$model |>
    add_group_rate("infectivity", by = c("age", "sex"),
                   prior = "gamma", shape = 2, scale = 1, joint = TRUE) |>
    sampled_params()
  params$prior |> expect_all_equal("gamma")
  params$args[[1]] |> expect_equal(c(shape = 2, scale = 1))
  params$min |> expect_all_equal(0)
  params$max |> expect_all_na()
  params$proposal |> expect_all_na()
})

test_that("Group rates load from a list of definitions", {
  tidy <- mk_rate$model |>
    add_group_rate(
      community = list(by = c("age", "sex"), prefix = "ct",
                       prior = "uniform", min = 0, max = 10),
      susceptibility = list(by = c("age", "sex"),
                            prior = "lognormal", sdlog = 3),
      # an empty prefix leaves the parameter names bare
      infectivity = list(by = "sex", prefix = "")
    ) |>
    get_group_rates()

  unique(tidy$group_rate) |>
    expect_setequal(c("community", "susceptibility", "infectivity"))
  tidy$param_name |> expect_match("^ct_", all = FALSE)
  tidy$param_name |> expect_match("^susceptibility_", all = FALSE)
  subset(tidy, group_rate == "infectivity")$param_name |>
    expect_setequal(c("sex_M", "sex_F"))
})

test_that("An unknown covariate in a group rate points to the close one", {
  mk_rate$model |>
    add_group_rate("susceptibility", by = "agee") |>
    expect_error("[Uu]nknown covariate.*agee.*[Dd]id you mean .age")
  mk_rate$model |>
    add_group_rate("susceptibility", by = "age",
                   reference = c(agee = "adult")) |>
    expect_error("[Uu]nknown covariate.*reference.*[Dd]id you mean .age")
})

test_that("A rate parameter keeps its definition across a dataset swap", {
  expect_warning(
    hh_model <- mk_flu$model |>
      add_group_rate("susceptibility", by = "member", reference = "student") |>
      add_parameter("susceptibility_member_mother",
                    prior = "lognormal", meanlog = 2, sdlog = 0.1),
    "[Oo]verwriting"
  )

  # the same data leaves every definition as declared, the override included
  hh_model |>
    add_dataset(mk_flu$data) |>
    get_parameters() |>
    expect_identical(get_parameters(hh_model))

  # data without mothers and fathers generates neither parameter:
  # the override is dropped with a warning, the default silently
  hh_data <- subset(mk_flu$data, !member %in% c("mother", "father"))
  expect_warning(hh_model_2 <- add_dataset(hh_model, hh_data),
                 "susceptibility_member_mother") |>
    conditionMessage() |>
    expect_no_match("father")
  get_parameters(hh_model_2)$name |> expect_no_match("mother|father")
})

####  Custom strata definitions  ####

test_that("An unknown name in params points to the strata table", {
  mk_rate$model |>
    add_group_rate("infectivity", strata = mk_rate$strata,
                   params = list(inf_isl = list(prior = "uniform"))) |>
    expect_error("[Uu]nknown parameter.*inf_isl.*[Dd]id you mean .inf_isol")
  mk_rate$model |>
    add_group_rate("infectivity", strata = mk_rate$strata,
                   params = list(unrelated = list(prior = "uniform"))) |>
    expect_error("[Uu]nknown parameter.*unrelated.*param_name")
})

test_that("Group rates deal with custom logic table correctly", {
  hh_model <- mk_rate$model |>
    add_group_rate("infectivity", strata = mk_rate$strata,
                   prior = "lognormal", meanlog = 0, sdlog = 2)

  tidy <- get_group_rates(hh_model)
  tidy$group_rate |> expect_all_equal("infectivity")
  # the custom strata drives the covariate column and the rows
  colnames(tidy) |> expect_contains("isolation")
  tidy$param_name |> expect_setequal(c(".reference", "inf_isol"))

  # parameters generated for non-reference only
  params <- get_parameters(hh_model)
  params$name |> expect_identical("inf_isol")
  params$prior |> expect_identical("lognormal")
  params$args[[1]] |> expect_equal(c(meanlog = 0, sdlog = 2))
})

test_that("Group rate parameter overrides work", {
  strata <- tibble::tribble(~age,    ~param_name,
                            "child", ".reference",
                            "adult", "sus_old",
                            "uncle", "sus_uncle")
  overrides <- list(
    sus_old = 1.73,
    sus_uncle = list(prior = "lognormal", sdlog = 1.2,
                     min = NA_real_, max = 10,
                     proposal = "rw_log", proposal_scale = 1)
  )
  hh_model <- mk_rate$data |>
    transform(age = c("child", "adult", "uncle", "adult")) |>
    new_housepi_model("final-size") |>
    add_group_rate("susceptibility", strata = strata,
                   params = overrides, prior = "gamma", shape = 2, scale = 1)

  params <- get_parameters(hh_model)
  params$name |> expect_setequal(c("sus_old", "sus_uncle"))
  # overrides used instead of raw params
  old <- get_parameters(hh_model, "sus_old")
  old$const |> expect_equal(1.73)
  uncle <- get_parameters(hh_model, "sus_uncle")
  uncle$prior |> expect_identical("lognormal")
  uncle$args[[1]] |> expect_equal(c(sdlog = 1.2, meanlog = 0))

  strata <- tibble::tribble(~sex, ~param_name,
                            "M",  "inf_male",
                            "F",  ".reference")
  hh_model |>
    add_group_rate("infectivity", strata = strata, params = list(rho = 2.4)) |>
    expect_error("[Uu]nknown.*rho")
})

test_that("An automatic definition takes no params list", {
  # `by` derives its own parameters, so the list has nothing to override
  mk_rate$model |>
    add_group_rate("infectivity", by = "sex", params = list(inf_male = 2.4)) |>
    expect_warning("[Pp]arams.*not used with")
})

test_that("Group rate mixed auto and custom modes preserve strata entries", {
  hh_model <- mk_rate$model |>
    add_group_rate(community = list(by = "age", prefix = "com"),
                   infectivity = list(strata = mk_rate$strata))

  tidy <- get_group_rates(hh_model)
  unique(tidy$group_rate) |>
    expect_setequal(c("community", "infectivity"))
  # strata entry is unchanged
  subset(tidy, group_rate == "infectivity")$param_name |>
    expect_setequal(c(".reference", "inf_isol"))
  # both modes materialise their parameters
  param_names <- get_parameters(hh_model)$name
  param_names |> expect_contains("inf_isol")
  param_names |> expect_match("^com_", all = FALSE)

  add_group_rate(hh_model, "infectivity", by = "age") |>
    expect_warning("[Oo]verwriting.*infectivity")
  add_group_rate(hh_model, "community", strata = mk_rate$strata) |>
    expect_warning("[Oo]verwriting.*community")

  # a covariate the data does not carry fails as soon as data is attached
  data.frame(hh_id = 1, status = 0, age = c("child", "adult")) |>
    new_housepi_model("final-size") |>
    add_group_rate("infectivity", strata = mk_rate$strata) |>
    expect_error("[Uu]nknown covariate.*isolation")
})

test_that("Absoluteness follows the reference across every definition mode", {
  # the default prior tracks whether a rate carries absolute values
  rate_prior <- function(...) {
    params <- get_parameters(add_group_rate(mk_rate$model, ...))
    unique(params$prior[is.na(params$const)])
  }
  tagged <- data.frame(sex = c("F", "M"),
                       param_name = c(".reference", "gentlemen"))
  untagged <- data.frame(sex = c("F", "M"),
                         param_name = c("ladies", "gentlemen"))

  # auto mode: absolute unless a reference is resolved
  rate_prior("susceptibility", by = "sex") |> expect_identical("uniform")
  rate_prior("susceptibility", by = "sex", reference = ".auto") |>
    expect_identical("lognormal")
  rate_prior("susceptibility", by = "sex", reference = c(sex = "M")) |>
    expect_identical("lognormal")
  # custom mode: the tag alone decides
  rate_prior("susceptibility", strata = untagged) |> expect_identical("uniform")
  rate_prior("susceptibility", strata = tagged) |> expect_identical("lognormal")
})

test_that("A strata row without a param_name is refused", {
  strata <- tibble::tribble(~age,    ~param_name,
                            "child", ".reference",
                            "adult", NA)
  mk_rate$model |>
    add_group_rate("susceptibility", strata = strata) |>
    expect_error("[Ss]trata row without param_name.*2")
})

test_that("Strata leaving members unmatched are refused", {
  # the single rule reaches the five students only
  strata <- data.frame(member = "student", param_name = "com_student")
  mk_flu$model |>
    add_group_rate("community", strata = strata) |>
    expect_error("[Gg]roup rate.*community.*7 members unmatched")

  # a reference does not stand in for the members no rule reaches
  strata <- tibble::tribble(~member,   ~param_name,
                            "student", ".reference",
                            "mother",  "com_mother")
  mk_flu$model |>
    add_group_rate("community", strata = strata) |>
    expect_error("5 members unmatched")
})

test_that("A catch-all rule makes an untagged table absolute", {
  strata <- tibble::tribble(~member,   ~param_name,
                            "student", "com_student",
                            NA,        "com_other")
  params <- mk_flu$model |>
    add_group_rate("community", strata = strata) |>
    get_parameters()

  subset(params, name == "alpha")$const |> expect_equal(1)
  subset(params, name != "alpha")$prior |> expect_all_equal("uniform")
})

test_that("A table setting no covariate value is refused", {
  # a single rate over every member duplicates the baseline
  mk_flu$model |>
    add_group_rate("community", strata = data.frame(param_name = "com")) |>
    expect_error("[Ss]trata table sets no covariate value")
  mk_flu$model |>
    add_group_rate("community",
                   strata = data.frame(member = NA, param_name = "com")) |>
    expect_error("[Ss]trata table sets no covariate value")
})

test_that("A custom table is overwritten like any other", {
  strata <- data.frame(member = c("student", NA),
                       param_name = c("com_student", "com_other"))
  mk_flu$model |>
    add_group_rate("community", strata = strata) |>
    add_group_rate("community", by = "member") |>
    expect_warning("[Oo]verwriting 1 group rate definition")
})

test_that("An automatic rate over a covariate missing everywhere is refused", {
  transform(mk_flu$data, datum = NA_character_) |>
    new_housepi_model("final-size") |>
    add_group_rate("community", by = "datum") |>
    expect_error("[Gg]roup rate.*community.*12 members unmatched")
})

test_that("A reference no member resolves to is refused", {
  # the tagged rule names a level the data lacks
  strata <- data.frame(
    member = c("uncle", "student", "other", "sibling", "mother", "father"),
    param_name = c(".reference", "com_student", "com_other", "com_sibling",
                   "com_mother", "com_father")
  )
  mk_flu$model |>
    add_group_rate("community", strata = strata) |>
    expect_error("1 rule no member resolves to.*uncle.*reference")

  # the catch-all reference is beaten by a specific rule for every member
  strata$member[1] <- NA
  mk_flu$model |>
    add_group_rate("community", strata = strata) |>
    expect_error("1 rule no member resolves to: catch-all.*reference")
})

test_that("Rules tied at equal specificity are refused", {
  strata <- data.frame(member = c("student", "student", NA),
                       param_name = c("com_student", "com_pupil", ".reference"))
  mk_flu$model |>
    add_group_rate("community", strata = strata) |>
    expect_error("[Gg]roup rate .community. has tied rules")
})

####  Strata pins  ####

test_that("A strata const pins the stratum, with no prior", {
  hh_model <- mk_flu$model |>
    add_group_rate("susceptibility", strata = pin_strata(),
                   prior = "lognormal", sdlog = 3)
  # a pinned stratum is a constant parameter: it carries no prior
  mum <- get_parameters(hh_model, "sus_mum")
  mum$const |> expect_equal(1.2)
  mum$prior |> expect_all_na()
  # an unpinned stratum keeps the shared prior
  get_parameters(hh_model, "sus_rest")$prior |> expect_identical("lognormal")

  # the pin column surfaces on the read table only when a pin exists
  tidy <- get_group_rates(hh_model)
  colnames(tidy) |> expect_contains("const")
  # padding keeps the column typed across sub-tables
  tidy$const |> expect_type("double")
  tidy$const[tidy$param_name == "sus_mum"] |> expect_equal(1.2)
  tidy$const[tidy$param_name == ".reference"] |> expect_all_na()

  # an unpinned rate keeps the table it had
  mk_flu$model |>
    add_group_rate("susceptibility", by = "member") |>
    get_group_rates() |>
    colnames() |>
    intersect("const") |>
    expect_length(0)
})

test_that("Invalid strata pins are rejected", {
  mk <- function(strata) {
    add_group_rate(mk_flu$model, "susceptibility", strata = strata,
                   prior = "lognormal", sdlog = 3)
  }

  # the reference is 1 by definition
  mk(pin_strata(const = c(0.9, 1.2, 0.5))) |>
    expect_error("[Rr]eference")
  # the column carries the value a stratum is pinned to
  mk(pin_strata(const = c(NA, TRUE, TRUE))) |>
    expect_error("const.*of the strata table must be numeric")
  mk(transform(pin_strata(), const = as.character(const))) |>
    expect_error("numeric")

  # a stratum spread over two rules carries a single pin, even when the data
  # carries no uncle to keep the second rule
  data.frame(member = c("student", "mother", "uncle"),
             param_name = c(".reference", "sus_mum", "sus_mum"),
             const = c(NA, 1.2, 0.4)) |>
    mk() |>
    expect_error("[Cc]onflicting.*const.*sus_mum.*rows 2 and 3")
})

test_that("A pinned stratum takes no entry in params", {
  mk_flu$model |>
    add_group_rate("susceptibility", strata = pin_strata(),
                   params = list(sus_mum = 2), prior = "lognormal",
                   sdlog = 3) |>
    expect_error("sus_mum")
})

####  Status strata  ####

test_that("Status is accepted as a stratification axis", {
  # both gates let 'status' through for a group rate
  mk_status$model |>
    add_group_rate("infectivity", by = "status") |>
    get_group_rates() |>
    getElement("status") |>
    as.character() |>
    expect_setequal(c("symptomatic", "asymptomatic"))
  strata <- tibble::tribble(~status,        ~param_name,
                            "symptomatic",  "inf_symptomatic",
                            "asymptomatic", ".reference")
  mk_status$model |>
    add_group_rate("infectivity", strata = strata) |>
    get_parameters() |>
    subset(is.na(const)) |>
    getElement("name") |>
    expect_identical("inf_symptomatic")

  # one parameter per non-reference stratum, named after the label
  mk_status$model |>
    add_group_rate(infectivity = list(
      by = "status", reference = list(status = "symptomatic")
    )) |>
    get_parameters() |>
    getElement("name") |>
    expect_identical("infectivity_status_asymptomatic")
})

test_that("A status stratum leaving a level uncovered is refused", {
  strata <- tibble::tribble(~status,       ~param_name,
                            "symptomatic", "inf_symptomatic")

  # any member may be drawn asymptomatic, escapees included: all 12 unmatched
  mk_status$model |>
    add_group_rate("infectivity", strata = strata) |>
    expect_error("[Gg]roup rate .infectivity. leaves 12 members unmatched")

  # whatever status the data records
  hh_data <- mk_status$data
  hh_data$status[2] <- NA
  hh_data |>
    new_housepi_model("final-size") |>
    add_group_rate("infectivity", strata = strata) |>
    expect_error("leaves 12 members unmatched")
})

test_that("A status stratum covers the codes a status process draws", {
  # the data records code 1 only, the status process also draws code 2
  hh_data <- data.frame(hh_id = rep(1:4, each = 3), status = rep(c(1, 0, 0), 4))
  hh_model <- hh_data |>
    new_housepi_model("final-size") |>
    add_status_process(0.5, 0.5)

  # any member may be drawn at code 2: all 12 unmatched
  hh_model |>
    add_group_rate("infectivity",
                   strata = data.frame(status = 1, param_name = "inf_one")) |>
    expect_error("leaves 12 members unmatched")
  hh_model |>
    add_group_rate("infectivity", by = "status") |>
    get_group_rates() |>
    getElement("status") |>
    as.character() |>
    expect_identical(c("1", "2"))

  # with no status recorded, the status process alone declares the codes
  hh_data[["status"]] <- NULL
  hh_data |>
    new_housepi_model("final-size") |>
    add_status_process(0.5, 0.5) |>
    add_group_rate("infectivity", by = "status") |>
    get_group_rates() |>
    getElement("status") |>
    as.character() |>
    expect_identical(c("1", "2"))
})

test_that("A status stratum covers combinations the data never realised", {
  # children are never asymptomatic, and the one mother not infected escaped
  hh_data <- data.frame(hh_id = c(1, 1, 2, 2),
                        status = c("symptomatic", "asymptomatic",
                                   "symptomatic", "escaped"),
                        member = c("child", "mother", "child", "mother"))
  hh_model <- new_housepi_model(hh_data, "final-size")

  # a table covering the observed combinations only leaves 4 members exposed
  hh_model |>
    add_group_rate("infectivity", strata = tibble::tribble(
      ~member,  ~status,        ~param_name,
      "child",  "symptomatic",  ".reference",
      "mother", "asymptomatic", "inf_ma"
    )) |>
    expect_error("leaves 4 members unmatched")

  # a rule for an unrealised combination stands, with no observed member
  rates <- hh_model |>
    add_group_rate("infectivity", strata = tibble::tribble(
      ~member,  ~status,        ~param_name,
      "child",  "symptomatic",  ".reference",
      "child",  "asymptomatic", "inf_ca",
      "mother", "symptomatic",  "inf_ms",
      "mother", "asymptomatic", "inf_ma"
    )) |>
    get_group_rates()
  rates$n[rates$param_name %in% c("inf_ca", "inf_ms")] |>
    expect_identical(c(0L, 0L))

  # the automatic grid carries every declared level on each member role
  rates <- hh_model |>
    add_group_rate(infectivity = list(by = c("member", "status"),
                                      joint = TRUE, reference = ".auto")) |>
    get_group_rates()
  rates |> expect_nrow(4L)
  sum(rates$n == 0L) |> expect_identical(2L)
})

test_that("Escaped members leave a status stratum scale free", {
  # escaped members never read infectivity, so no status rule anchors it
  strata <- tibble::tribble(~datum,        ~param_name,
                            "symptomatic", "sus_symptomatic",
                            NA,            "sus_rest")
  mk_status$model |>
    add_group_rate(infectivity = list(by = "status"),
                   susceptibility = list(strata = strata)) |>
    add_parameter(alpha = 0.1, gamma = 1, hh_size_ref = 4) |>
    evaluate_log_lik() |>
    expect_error("2 absolute rule groups.*susceptibility.*infectivity")
})

test_that("Only the rates the model declares may be status stratified", {
  # infectivity is read after the member's status is settled, the others are
  # read on members that carry no infected status
  mk_status$model |>
    add_group_rate("susceptibility", by = "status") |>
    expect_error("cannot be stratified by the infected status")
  mk_status$model |>
    add_group_rate("community", by = "status") |>
    expect_error("cannot be stratified by the infected status")

  # contact rates keep no status axis at all
  mk_status$model |>
    add_contact_rate(by = "status") |>
    expect_error("reserved")

  # a model that declares no status level cannot stratify on it. The rate
  # waits for the declaration, so decorator order is free, and the wall reports
  # that it never came
  hh_model <- add_group_rate(mk_rate$model, "infectivity", by = "status")
  get_group_rates(hh_model)$group_rate |> expect_identical("infectivity")
  hh_model |>
    add_parameter(alpha = 0.01, beta = 0.2, gamma = 1, hh_size_ref = 3) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("model declares no level")
})

test_that("Status strata are declared, not read off the data", {
  # the observed labels are the model dictionary: a stratum outside it
  # is refused
  hh_data <- mk_status$data
  hh_data$status[hh_data$status == "asymptomatic"] <- "escaped"
  strata <- tibble::tribble(~status,        ~param_name,
                            "symptomatic",  ".reference",
                            "asymptomatic", "inf_asymptomatic")
  hh_data |>
    new_housepi_model("final-size") |>
    add_group_rate("infectivity", strata = strata) |>
    expect_error("[Ss]tatus label .asymptomatic. is not in the data")

  # with no status observed at all, no sample size can pick a reference
  hh_model <- mk_status$data |>
    dplyr::mutate(status = NA) |>
    new_housepi_model("final-size") |>
    add_latent(infec_time = list(model = "uniform", min = 0, max = 5)) |>
    add_status_process(symptomatic = 0.3, asymptomatic = 0.7)
  hh_model |>
    add_group_rate(infectivity = list(by = "status", reference = ".auto")) |>
    expect_error("reference = list\\(status")

  # a custom table stands all the same, its sample sizes unknown
  strata
  hh_model |>
    add_group_rate("infectivity", strata = strata) |>
    get_group_rates() |>
    getElement("n") |>
    expect_identical(c(NA_integer_, NA_integer_))
})

test_that("A status stratum reads as a covariate duplicating it", {
  # the workaround a status axis replaces: a data column carrying the status
  # of the infected members, NA elsewhere. Escapees read a catch-all rule, which
  # the likelihood never evaluates for them
  pin <- function(hh_model, name) {
    hh_model <- add_parameter(hh_model, alpha = 0.01, beta = 0.2,
                              gamma = 1, hh_size_ref = 3)
    suppressWarnings(add_parameter(hh_model, import = stats::setNames(
      as.list(rep(0.4, length(name))), name
    )))
  }
  by_status <- mk_status$model |>
    add_group_rate(infectivity = list(
      by = "status", reference = list(status = "symptomatic")
    )) |>
    pin("infectivity_status_asymptomatic")
  by_datum <- mk_status$model |>
    add_group_rate("infectivity", strata = tibble::tribble(
      ~datum,         ~param_name,
      "symptomatic",  ".reference",
      "asymptomatic", "inf_datum_asymptomatic",
      NA,             "inf_datum_escaped"
    )) |>
    pin(c("inf_datum_asymptomatic", "inf_datum_escaped"))

  # every parameter is pinned, so the likelihood is a deterministic probe
  get_parameters(by_status)$const |> expect_no_na()
  evaluate_log_lik(by_status) |> expect_identical(evaluate_log_lik(by_datum))
})

test_that("A status stratum drives simulated transmission", {
  # the drawn status, not a frozen covariate, selects the infectivity: with the
  # asymptomatic stratum near zero, a run that draws only asymptomatic cases
  # must transmit far less than one that draws only symptomatic ones
  hh_model <- sarscov2_2021_icofs |>
    dplyr::filter(hh_id %in% sample(unique(hh_id), 150)) |>
    dplyr::mutate(status = NA, end_date = 21) |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_infectivity_profile(import = "template") |>
    add_status_process(symptomatic = 0.5, asymptomatic = 0.5) |>
    add_group_rate(infectivity = list(
      by = "status", reference = list(status = "symptomatic")
    )) |>
    add_parameter(alpha = 1e-12, beta = 20, gamma = 1, hh_size_ref = 4,
                  infectivity_status_asymptomatic = 1e-3) |>
    suppressWarnings()

  # the status weights must reach the categorical as literals
  attack_rate <- function(symptomatic) {
    simulated <- hh_model |>
      add_status_process(symptomatic = !!symptomatic,
                         asymptomatic = !!(1 - symptomatic)) |>
      run_simulation() |>
      suppressWarnings()
    mean(as.character(simulated$status) != "escaped")
  }

  # alpha is negligible, so every case beyond the introductions is transmitted
  attack_rate(0.9) |> expect_gt(1.5 * attack_rate(0.1))
})

test_that("Status strata count the infected members only", {
  # the status is read off the member, not off a data column: escaped members
  # carry no infected label, so they land in no stratum
  tbl <- mk_icofs$model |>
    add_parameter(import = "template") |>
    add_group_rate(infectivity = list(by = "status", reference = ".auto")) |>
    get_group_rates()

  sum(tbl$n) |> expect_identical(sum(sarscov2_2021_icofs$status > 0))
  as.character(tbl$status) |> expect_setequal(c("1", "2"))
})

####  End-to-end  ####

test_that("An absolute rate collapses the baseline it scales", {
  # alpha is never declared: the community rate carries the whole rate
  hh_model <- mk_flu$model |>
    add_group_rate("community", by = "member") |>
    add_parameter(beta = 0.1, gamma = 1, hh_size_ref = 4)

  evaluate_log_lik(hh_model) |> expect_all_finite()
  params <- get_parameters(hh_model)
  params$const[params$name == "alpha"] |> expect_equal(1.)

  # the collapsed baseline is carried into the posterior
  draws <- hh_model |>
    run_inference(steps = 5L, thin = 1L, adapt = 0L, chains = 1L) |>
    suppressWarnings() |>
    summarise_draws()
  draws$variable |> expect_disjoint("alpha")
  draws$variable |> expect_match("^community_member", all = FALSE)
})

test_that("A baseline pinned by an absolute rate follows its redeclaration", {
  hh_model <- mk_flu$model |>
    add_group_rate("community", by = "member", reference = "student") |>
    add_parameter(beta = 0.1, gamma = 1, hh_size_ref = 4)
  # the same rate declared absolute first, then relative
  hh_model_2 <- mk_flu$model |>
    add_group_rate("community", by = "member") |>
    add_parameter(beta = 0.1, gamma = 1, hh_size_ref = 4)
  expect_warning(
    hh_model_2 <- hh_model_2 |>
      add_group_rate("community", by = "member", reference = "student"),
    "[Oo]verwriting.*community"
  )

  by_name <- \(x) x[order(x$name), , drop = FALSE]
  by_name(get_parameters(hh_model_2)) |>
    expect_equal(by_name(get_parameters(hh_model)), ignore_attr = TRUE)
  get_parameters(hh_model_2)$name |> expect_disjoint("alpha")
  hh_model_2 |>
    evaluate_log_lik() |>
    expect_error("[Mm]issing parameter definition.*alpha")
})

test_that("Two absolute rates on one baseline are identifiable once pinned", {
  hh_model <- mk_flu$model |>
    add_group_rate("susceptibility", by = "member") |>
    add_group_rate("infectivity", by = "member") |>
    add_parameter(alpha = 0.01, gamma = 1, hh_size_ref = 4)
  hh_model |>
    evaluate_log_lik() |>
    expect_error("[Nn]on-identifiability.*susceptibility.*infectivity")

  # a simulation table gives every scale a value: no product of free scales
  # is left to resolve
  run_simulation(hh_model, draws = pin_all_params(hh_model)) |>
    expect_s3_class("data.frame")
})

test_that("An absolute group rate holds its baseline at the reference", {
  # `community` names no reference, so it absorbs the scale `alpha` carried
  # even if baseline gets outside of its proposed prior domain
  expect_warning(
    hh_model <- mk_flu$model |>
      add_parameter(alpha = list(prior = "uniform", max = 0.7),
                    beta = 0.1, gamma = 1, hh_size_ref = 4) |>
      add_group_rate("community", by = "member"),
    "[A]bsolute rate.*community.*no reference"
  )

  # one stratum read from the table fixes the level of the whole group
  hh_model |>
    run_simulation(draws = data.frame(community_member_mother = 0.5)) |>
    expect_s3_class("data.frame")
  # the prior it was given leaves with the scale
  params <- get_parameters(hh_model)
  params$const[params$name == "alpha"] |> expect_equal(1.)

  # a stratum held constant anchors the same way
  hh_model |>
    add_parameter(community_member_student = 1.7) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_all_finite()
})

test_that("A declared baseline is collapsed at the decorator", {
  # the collapse is settled where the rate is declared, not at the run
  expect_warning(
    hh_model <- mk_flu$model |>
      add_parameter(alpha = list(prior = "uniform", min = 0, max = 1),
                    beta = 0.1, gamma = 1, hh_size_ref = 4) |>
      add_group_rate("community", by = "member"),
    "[Bb]aseline.*constant reference"
  )

  run_simulation(hh_model, draws = pin_all_params(hh_model)) |>
    expect_s3_class("data.frame")

  # a later definition on an untouched baseline says nothing again
  hh_model |>
    add_parameter(gamma = 2) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_all_finite()

  # the run cannot give the scale back
  hh_model |>
    run_simulation(draws = data.frame(alpha = 0.5)) |>
    expect_error("draws.*alpha.*other than 1")
  hh_model |>
    run_simulation(draws = data.frame(alpha = 1)) |>
    expect_s3_class("data.frame")
})

test_that("A collapsed baseline is reported again only on a new declaration", {
  hh_model <- mk_flu$model |>
    add_parameter(beta = 0.1, gamma = 1, hh_size_ref = 4) |>
    add_parameter(alpha = list(prior = "uniform", min = 0, max = 1)) |>
    suppressWarnings()
  expect_warning(
    hh_model <- add_group_rate(hh_model, "community", by = "member"),
    "[Bb]aseline.*constant reference"
  )

  # an untouched baseline and a copy of the model say nothing more
  hh_model |>
    add_parameter(gamma = 2) |>
    expect_no_warning(message = "[Bb]aseline") |>
    expect_warning("[Oo]verwriting.*gamma")
  hh_model$copy() |>
    add_parameter(gamma = 2) |>
    expect_no_warning(message = "[Bb]aseline") |>
    expect_warning("[Oo]verwriting.*gamma")

  # a new range for the same prior is a new declaration
  hh_model |>
    add_parameter(alpha = list(prior = "uniform", min = 0, max = 2)) |>
    expect_warning("[Bb]aseline.*constant reference") |>
    expect_warning("[Oo]verwriting.*alpha")

  # a baseline its rate released is reported again once re-collapsed
  released <- hh_model |>
    add_group_rate("community", by = "member", reference = "student") |>
    suppressWarnings()
  get_parameters(released, "alpha")$const |> expect_all_na()
  add_group_rate(released, "community", by = "member") |>
    expect_warning("[Bb]aseline.*constant reference") |>
    expect_warning("[Oo]verwriting.*community")
})

test_that("A measure table cannot give a collapsed baseline its scale back", {
  hh_model <- mk_flu$model |>
    add_parameter(beta = 0.1, gamma = 1, hh_size_ref = 4) |>
    add_group_rate("community", by = "member")
  hh_model |>
    check_predictive(
      data.frame(alpha = 0.5, beta = 0.1, gamma = 1, hh_size_ref = 4L,
                 community_member_mother = 0.5),
      sar = measure_sar_by_hh_size, thin_ess = FALSE, nrep = 2L
    ) |>
    expect_error("draws.*alpha.*other than 1")
})

test_that("A strata pin reaches simulation and stratifies on data only", {
  # every scale reads a value, and the pin columns are no stratification axis
  mk_flu$pinned |>
    add_group_rate("susceptibility", strata = pin_strata(),
                   prior = "lognormal", sdlog = 3) |>
    run_simulation() |>
    expect_s3_class("data.frame")
})

test_that("A pinned stratum holds across draws, an unpinned one is resampled", {
  mk <- function(strata) {
    mk_flu$pinned |>
      add_group_rate("susceptibility", strata = strata,
                     prior = "lognormal", sdlog = 3)
  }
  twice <- \(hh_model) c(evaluate_log_lik(hh_model), evaluate_log_lik(hh_model))

  # every stratum pinned: inference reads the values, nothing is drawn
  log_lik <- twice(mk(pin_strata(const = c(NA, 1.2, 0.5, 0.8))))
  log_lik[1] |> expect_equal(log_lik[2])

  # an unpinned stratum is drawn from its prior
  log_lik <- twice(mk(pin_strata()))
  isTRUE(all.equal(log_lik[1], log_lik[2])) |> expect_false()
})

test_that("A fixed stratum anchors the scale of its rule group", {
  # the group is anchored by the constant, so no baseline is collapsed onto it
  mk_flu$model |>
    add_parameter(alpha = list(prior = "uniform", min = 0.2, max = 0.9),
                  beta = 0.1, gamma = 1, hh_size_ref = 4) |>
    add_group_rate("community",
      strata = pin_strata(),
      prior = "lognormal", sdlog = 3
    ) |>
    evaluate_log_lik() |>
    expect_all_finite()
})

test_that("Strata counts come from the engine cross-tabulation", {
  # a cell the data does not carry forms no stratum of its own: the member
  # holding it is left unmatched
  hh_model <- data.frame(hh_id = c("h1", "h1", "h2", "h2", "h3"),
                         status = c(1, 0, 1, 1, 0),
                         datum = c("a", "b", "a", NA, "b"),
                         datum_2 = c(TRUE, FALSE, TRUE, FALSE, NA)) |>
    new_housepi_model("final-size", quiet = TRUE) |>
    add_parameter(import = "template")

  # susceptibility reads the member missing its datum: no stratum resolves it
  hh_model |>
    add_group_rate(susceptibility = list(by = "datum", reference = ".auto")) |>
    expect_error("leaves 1 member unmatched.*datum.*values are missing")

  # infectivity reads the escaped member missing its datum_2 too: a
  # simulation may infect it
  hh_model |>
    add_group_rate(infectivity = list(by = "datum_2", reference = ".auto")) |>
    expect_error("leaves 1 member unmatched")

  # a logical column keeps both of its levels
  datum_2 <- data.frame(hh_id = c("h1", "h1", "h2"), status = c(1, 0, 1),
                        datum_2 = c(TRUE, FALSE, TRUE)) |>
    new_housepi_model("final-size", quiet = TRUE) |>
    add_group_rate(infectivity = list(by = "datum_2", reference = ".auto")) |>
    get_group_rates("infectivity")
  as.character(datum_2$datum_2) |> expect_setequal(c("TRUE", "FALSE"))
  sum(datum_2$n) |> expect_identical(3L)
})

test_that("A catch-all stratum survives a missing covariate cell", {
  # NA reads as missing in a data cell and as a catch-all in a rule cell: a
  # column that happens to carry an empty cell must not collapse the two
  strata <- tibble::tribble(~datum, ~param_name,
                            "a",    ".reference",
                            NA,     "sus_rest")
  hh_model <- data.frame(id = c("m1", "m2", "m3"), hh_id = "h1",
                         status = c(1, 0, 0), datum = c("a", "b", NA)) |>
    new_housepi_model("final-size", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_group_rate("susceptibility", strata = strata,
                   params = list(sus_rest = 3))
  snapshot <- hh_model |> inspect_state("h1")

  rate <- stats::setNames(snapshot$group_rates$susceptibility,
                          snapshot$group_rates$id)

  # the specific stratum sets the unit scale
  unname(rate["m1"]) |> expect_equal(1)
  # the catch-all covers every other member, the ones carrying a value...
  unname(rate["m2"]) |> expect_equal(3)
  # ...and the ones that do not
  unname(rate["m3"]) |> expect_equal(3)

  # the catch-all stratum counts the members it resolves to
  rates <- get_group_rates(hh_model)
  rates$n[rates$param_name == ".reference"] |> expect_identical(sum(rate == 1))
  rates$n[rates$param_name == "sus_rest"] |> expect_identical(sum(rate == 3))
})

####  Observation collision  ####

test_that("A group rate keyed on an observation is rejected", {
  # the rate feeds the transmission model the observation is generated from
  mk_observed$model |>
    add_group_rate("susceptibility", by = "datum") |>
    add_observation("datum", model = "uniform", min = 0, max = 5) |>
    expect_error("keyed on observation.*datum")

  # the reverse declaration order fails just the same
  mk_observed$model |>
    add_observation("datum", model = "uniform", min = 0, max = 5) |>
    add_group_rate("susceptibility", by = "datum") |>
    expect_error("keyed on observation.*datum")

  # a custom strata table takes no loophole
  strata <- tibble::tribble(~datum, ~param_name,
                            2,      ".reference",
                            3,      "sus_high")
  mk_observed$model |>
    add_observation("datum", model = "uniform", min = 0, max = 5) |>
    add_group_rate("susceptibility", strata = strata) |>
    expect_error("keyed on observation.*datum")

  # a deterministic household observation standing in for its data column is
  # derived, an observation all the same
  mk_observed$model |>
    add_observation(datum = min(datum_2), level = "house") |>
    add_group_rate("susceptibility", by = "datum") |>
    expect_error("keyed on observation.*datum")
})
