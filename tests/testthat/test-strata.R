library(testthat)

####  Helpers  ####

# three members over `datum`, two of them sharing a value
mk_strata <- local({
  hh_data <- data.frame(datum = c("a", "a", "b"))
  list(data = hh_data, cells = member_cells_(hh_data, "datum"))
})

# fit check of one table keyed on `datum`, its rules in one rule group
fit_datum <- function(tbl, cells = mk_strata$cells) {
  tbl$rule_group <- 1L
  assert_strata_fit_(list(tbl), cells, c("param_name", "rule_group"),
                     "param_name", "The table", "unit")
}

####  Definitions  ####

test_that("A strata table needs rows and its required columns", {
  assert_strata_shape_(list(datum = "a"), "datum", "The table", "") |>
    expect_error("[Tt]he table must be a data.frame with at least one row")
  assert_strata_shape_(data.frame(datum = character()), "datum",
                       "The table", "") |>
    expect_error("[Tt]he table must be a data.frame with at least one row")
  assert_strata_shape_(data.frame(datum = "a"), c("datum", "param_name"),
                       "The table", "") |>
    expect_error("[Tt]he table has no column.*param_name")
  assert_strata_shape_(data.frame(datum = "a"), "datum", "The table", "") |>
    expect_null()
})

test_that("A covariate takes no reserved name and is a data column", {
  assert_valid_covariates(c("datum", "hh_id"), c("id", "hh_id")) |>
    expect_error("[Ii]mpossible to stratify on covariate.*hh_id")
  assert_valid_covariates(c("datum", NA), c("id", "hh_id")) |>
    expect_null()
  assert_data_covariates(c("datum", "datum_2"), "datum") |>
    expect_error("[Uu]nknown covariate.*datum_2")
})

test_that("A rule naming a value the data never carries is flagged", {
  # a catch-all cell names no value
  data.frame(datum = c("a", "z", NA)) |>
    unknown_strata_rows(list(datum = c("a", "b")), c(datum = "datum")) |>
    expect_identical(c(FALSE, TRUE, FALSE))
})

test_that("Rule groups follow the maximal patterns of the table", {
  # half-keyed rules and the catch-all fall in the one maximal pattern
  tbl <- data.frame(datum = c("a", "a", NA, NA), datum_2 = c("x", NA, "y", NA),
                    param_name = c("rho", "rho_2", "rho_3", "rho_4"))
  compute_rule_groups(tbl, c("datum", "datum_2"), "param_name")$rule_group |>
    unique() |>
    expect_length(1L)

  # disjoint patterns make one group each
  tbl <- data.frame(datum = c("a", "b", NA, NA), datum_2 = c(NA, NA, "x", "y"),
                    param_name = c("rho", "rho_2", "rho_3", "rho_4"))
  compute_rule_groups(tbl, c("datum", "datum_2"), "param_name")$rule_group |>
    unique() |>
    expect_length(2L)

  # rho_3 keys on datum alone: it fits {datum, datum_2} and {datum, datum_3}
  tbl <- data.frame(datum = "a", datum_2 = c("x", NA, NA),
                    datum_3 = c(NA, "y", NA),
                    param_name = c("rho", "rho_2", "rho_3"))
  compute_rule_groups(tbl, c("datum", "datum_2", "datum_3"), "param_name") |>
    expect_error("[Aa]mbiguous rule 3 .*rho_3.*fitting 2")
})

####  Resolution  ####

test_that("Member cells hold each covariate combination once", {
  hh_data <- data.frame(datum = c("a", "a", "b", NA))
  cells <- member_cells_(hh_data, "datum")
  cells$datum |> expect_identical(c("a", "b", NA))
  cells$.weight |> expect_identical(c(2L, 1L, 1L))
  cells$.combo |> expect_identical(1:3)

  # every infected status repeats every combination
  cells <- member_cells_(hh_data, "datum", c("symptomatic", "asymptomatic"))
  cells$status |> expect_identical(rep(c("symptomatic", "asymptomatic"),
                                       each = 3))
  cells$.combo |> expect_identical(rep(1:3, 2))

  # no covariate: one cell holding every member
  member_cells_(hh_data, character())$.weight |> expect_identical(4L)
})

test_that("A cell goes to the most specific rule it matches", {
  tbl <- data.frame(datum = c("a", NA), datum_2 = c("x", NA))
  cells <- data.frame(datum = c("a", "a", "b", NA),
                      datum_2 = c("x", "y", "x", "x"))
  rule_winners_(tbl, cells, c(1L, 1L)) |>
    expect_identical(cbind(c(TRUE, FALSE, FALSE, FALSE),
                           c(FALSE, TRUE, TRUE, TRUE)))
  # a missing data cell matches the catch-all only
  rule_hits_(tbl, cells)[4, ] |> expect_identical(c(FALSE, TRUE))
})

test_that("Rules tied at the top specificity all win the cell", {
  tbl <- data.frame(datum = c("a", NA), datum_2 = c(NA, "x"))
  cells <- data.frame(datum = "a", datum_2 = "x")
  rule_winners_(tbl, cells, c(1L, 1L)) |>
    as.vector() |>
    expect_identical(c(TRUE, TRUE))
})

test_that("A cell is unmatched when some rule group misses it", {
  tbl <- data.frame(datum = c("a", NA), datum_2 = c(NA, "x"))
  cells <- data.frame(datum = c("a", "a", "b"), datum_2 = c("x", "y", "x"))
  unmatched_cells_(tbl, cells, c(1L, 2L)) |>
    expect_identical(c(FALSE, TRUE, TRUE))
  # a table without rules matches no cell
  unmatched_cells_(tbl[0, ], cells, integer()) |>
    expect_identical(rep(TRUE, 3))
})

test_that("A member counts once per rule across the infected statuses", {
  cells <- member_cells_(mk_strata$data, "datum", c("1", "2"))
  tbl <- data.frame(datum = c(NA, "a"), status = c(NA, "1"), rule_group = 1L)
  # the catch-all wins "a" at status 2 and "b" at both
  rule_counts_(tbl, cells, "rule_group") |> expect_identical(c(3L, 2L))
})

####  Fit checks  ####

test_that("The fit check names the table and the unit it fails", {
  data.frame(datum = "a", param_name = "rho") |>
    fit_datum() |>
    expect_error("[Tt]he table leaves 1 unit unmatched")
  data.frame(datum = c("a", "b", "z"),
             param_name = c("rho", "rho_2", "rho_3")) |>
    fit_datum() |>
    expect_error("[Tt]he table has 1 rule no unit resolves to.*z.*rho_3")
  data.frame(datum = c("a", "a", "b"),
             param_name = c("rho", "rho_2", "rho_3")) |>
    fit_datum() |>
    expect_error("[Tt]he table has tied rules.*rho.*rho_2")
})

test_that("A table resolving every unit returns its winners", {
  win <- data.frame(datum = c("a", NA), param_name = c("rho", "rho_2")) |>
    fit_datum()
  win[[1]] |> expect_identical(cbind(c(TRUE, FALSE), c(FALSE, TRUE)))
})

test_that("The unmatched hint defaults to adding a covering rule", {
  stop_unmatched_(0, "The table", "unit") |> expect_null()
  stop_unmatched_(2, "The table", "unit") |>
    expect_error("[Tt]he table leaves 2 units.*matches any value")
  # braces in a rendered hint stay literal
  stop_unmatched_(2, "The table", "unit", "Rename {datum}.") |>
    expect_error("[Tt]he table leaves 2 units.*Rename \\{datum\\}")
})

####  Engine coding  ####

test_that("Covariates code as data levels, labelled statuses as engine codes", {
  tbl <- data.frame(datum = c("b", NA), status = c("asymptomatic", NA))
  encode_strata_levels_(tbl, list(datum = list(levels = c("a", "b"))),
                        "datum")$datum |>
    expect_identical(c(2L, NA))
  code_status_col_(tbl, c(symptomatic = 1L, asymptomatic = 2L))$status |>
    expect_identical(c(2L, NA))
  # integer statuses stand for themselves
  code_status_col_(data.frame(status = c(3L, NA)), c(symptomatic = 1L)) |>
    getElement("status") |>
    expect_identical(c(3L, NA))
})
