library(testthat)

# Rd source of a template help page: man/ under load_all, help db when installed
template_rd <- function(name) {
  topic <- paste0("template_", name, ".Rd")
  src <- system.file("man", topic, package = "housepi")
  if (nzchar(src)) return(paste(readLines(src), collapse = "\n"))
  paste(as.character(tools::Rd_db("housepi")[[topic]]), collapse = "")
}

################################
###  TEMPLATE : DEFINITIONS  ###
################################

test_that("set_source_name_ assigns names to vector", {
  result <- set_source_name_(c("a", "b"), "param")
  names(result) |> expect_identical(c("param", "param"))
  unname(result) |> expect_identical(c("a", "b"))

  # NULL input: returns NULL
  set_source_name_(NULL, "param") |>
    expect_null()
  # NA name: unchanged
  names(set_source_name_(c("a", "b"), NA)) |>
    expect_null()
  # empty string name: unchanged
  names(set_source_name_(c("a"), "")) |>
    expect_null()
})

test_that("Missing template definitions are found and reported", {
  required <- c(A = "a", B = "b", C = "c")
  defined <- c("a", "b")

  # partial overlap
  template_missing_find(required, defined) |>
    expect_identical("c")
  # all defined
  template_missing_find(defined, defined) |>
    expect_identical(character(0))
  # no required
  template_missing_find(character(0), c("a")) |>
    expect_null()
  # no defined: returns all required sorted
  template_missing_find(rev(defined), character(0)) |>
    expect_identical(defined)

  results <- template_missing_stop(required, defined, "parameter")
  results |> expect_match("missing", all = FALSE)

  # no missing: empty output
  template_missing_stop(c("a"), c("a"), "parameter") |>
    expect_null()
})

test_that("Unused template definitions are found and warned about", {
  required <- c(A = "a")
  defined <- c("a", "b", "c")

  # all used: returns empty (character(0))
  template_not_used_find(defined, defined) |>
    expect_length(0L)
  # all used: no warning
  template_not_used_warn(defined, defined, "parameter") |>
    expect_no_warning()

  # partial overlap
  template_not_used_find(required, defined) |>
    expect_identical(c("b", "c"))
  # no required: returns all defined sorted
  template_not_used_find(character(0), rev(defined)) |>
    expect_identical(defined)
  # no defined
  template_not_used_find(c("a"), character(0)) |>
    expect_null()

  # warning for unused
  template_not_used_warn(c("a"), defined, "parameter") |>
    expect_warning("unused")
})

test_that("Overwritten template definitions are warned about", {
  # empty: no warning
  template_overwrite_warn(character(0), "parameter") |>
    expect_no_warning()

  template_overwrite_warn(c("x"), "parameter") |>
    expect_warning("Overwriting")

  # warn handles many elements
  template_overwrite_warn(letters[1:10], "parameter") |>
    expect_warning("Overwriting 10 parameter definitions.*j")
})

test_that("fmt_assert_cli formats messages correctly", {
  # with messages (non-minor): rewrites "were" to "are"
  result <- fmt_assert_cli(
    "Parameters", c("There were 2 missing parameter definitions: 'a', 'b'.")
  )
  result |> expect_match("Parameters", all = FALSE)
  result |> expect_match("are", all = FALSE)

  # minor mode: rewrites "Required by" to "Can be specified for"
  result <- fmt_assert_cli(
    "Parameters", c("There was 1 definition.", "Required by ..."),
    minor = TRUE
  )
  result |> expect_match("Can be specified", all = FALSE)
})

test_that("Template name normalization", {
  for (f in c("Final-Size", "final_size", " final siZe ", "FINAL SIZE")) {
    spaced_choice(f) |> expect_identical("final size")
  }

  result <- find_closest_template_name("final siz")
  result |> expect_match("Did you mean", all = FALSE)
  result <- find_closest_template_name("")
  result |> expect_match("browse_template", all = FALSE)
})

#############################
###  TEMPLATE : BROWSING  ###
#############################

test_that("browse_templates branches inform correctly", {
  # print existing templates
  browse_templates() |> expect_snapshot()
  # non-character name
  browse_templates(42) |>
    expect_error("not a character")
  # unknown name: error mentions "not a known template" and a hint
  browse_templates("unknown") |>
    expect_error("not a known template")
  # unknown name close to a real one: hint should suggest "Did you mean"
  browse_templates("final siz") |>
    expect_error("Did you mean")
})

test_that("browse_templates triggers help preview when path is NA", {
  # mock the internal helper so the test never spawns a real previewer
  called_with <- NULL
  with_mocked_bindings(
    open_model_helper_ = function(path) {
      called_with <<- path
      invisible()
    },
    code = browse_templates("final-size")
  )
  called_with |> expect_identical("final_size")
})

test_that("Every template has a help page", {
  for (name in list.files(template_models_dir())) {
    template_rd(name) |> expect_match("\\\\name\\{template_")
  }
})

test_that("browse_templates copies template to a new path", {
  dest <- tempfile("template_")
  on.exit(unlink(dest, recursive = TRUE), add = TRUE)
  model_cpp <- file.path(dest, "model.cpp")

  # first copy: creates the model sources
  browse_templates("final-size", path = dest)
  file.exists(model_cpp) |> expect_true()
  file.exists(file.path(dest, "param.R")) |> expect_true()
  # final_size has no latent file: nothing should be copied
  file.exists(file.path(dest, "latent.R")) |> expect_false()
  shipped <- readLines(model_cpp)

  # second copy without replace: warns and does not overwrite
  writeLines("edited", model_cpp)
  browse_templates("final-size", path = dest) |>
    expect_warning("already exists")
  readLines(model_cpp) |> expect_identical("edited")

  # with replace = TRUE: overwrites silently
  browse_templates("final-size", path = dest, replace = TRUE)
  readLines(model_cpp) |> expect_identical(shipped)
})

test_that("browse_templates rejects '.' as destination path", {
  sandbox <- tempfile("sandbox_")
  dir.create(sandbox)
  on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)
  withr::with_dir(sandbox, {
    browse_templates("final-size", path = ".") |>
      expect_warning("current directory")
    # helper redirects "." to a folder named after the template
    dir.exists("final-size") |> expect_true()
  })
})

test_that("Model help pages name the keys the engine registers", {
  modules <- c(Data = "Data", Param = "Parameters",
               Latent = "Latent", GroupRisk = "Group Risks")
  for (name in c("final_size", "case_ascertained")) {
    cpp <- system.file("models", name, "model.cpp", package = "housepi") |>
      readLines() |>
      paste(collapse = "\n")
    rd <- template_rd(name)
    # loop over all modules
    for (k in names(modules)) {
      # EMPTY_KEYS declarations match nothing and register no key
      fmt <- "(?<=REGISTER_KEYS\\(%s,)[^)]*"
      decl <- regmatches(cpp, regexpr(sprintf(fmt, k), cpp, perl = TRUE))
      keys <- trimws(unlist(strsplit(decl, ",")))
      # section body up to the next section: its \code{} items are the keys
      fmt <- "^.*\\\\section\\{%s\\}\\{(.*?)\\\\section.*$"
      section <- sub(sprintf(fmt, modules[[k]]), "\\1", rd)
      if (identical(section, rd)) section <- ""
      pattern <- "(?<=\\\\item \\\\code\\{)[^}]+"
      listed <- regmatches(section, gregexpr(pattern, section, perl = TRUE))
      listed[[1]] |> expect_setequal(keys)
    }
  }
})
