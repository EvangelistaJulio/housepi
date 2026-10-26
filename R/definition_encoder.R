## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Reserved constant names forbidden as definition names.
.reserved_math_consts <- c("pi")

# Coerce an R Date value to its numeric day count.
coerce_date_ <- function(x) {
  if (inherits(x, "Date")) as.numeric(unclass(x))
  else x
}

# Evaluate a self-contained constant (sub)expression to a scalar number.
# * base-R constant, as.Date() or as.integer() -> folded;
# * Date result -> its day count;
# * housepi node reference -> fails to evaluate, stays symbolic.
try_fold_const_ <- function(expr, env) {
  val <- tryCatch(suppressWarnings(coerce_date_(eval(expr, envir = env))),
                  error = \(e) NULL)
  if (is.numeric(val) && length(val) == 1L && !is.na(val)) val else NULL
}

# Plain name of every definition source, keyed by its source letter.
source_labels_ <- c(A = "latent variable",
                    D = "data column",
                    G = "group risk",
                    P = "parameter")

# Encode leaf values with their source:
# * number -> numeric;
# * name -> left unresolved, only the model knows what declares it.
housepi_format_sources <- function(func_args) {
  # default return
  if (is.null(func_args) || all(is.na(func_args))) return(NA_character_)

  values <- as.character(func_args)
  is_num <- is.numeric(func_args) |
    !is.na(suppressWarnings(as.numeric(values))) |
    values %in% .reserved_math_consts
  for (y in values[!is_num & !is.na(values)]) {
    assert_valid_character(y, "value")
  }
  sources <- ifelse(is_num, "N", "?")
  sources[is.na(values)] <- NA_character_
  sprintf("%s|%s", values, sources)
}

# Assert the argument keys of a call respect the naming convention.
# Positional arguments carry an empty key.
assert_valid_arg_keys <- function(keys) {
  for (k in keys[nzchar(keys)]) assert_valid_character(k, "argument name")
}

# Format a call over encoded arguments as a typed housepi definition.
housepi_format_call <- function(func_name, func_args) {
  if (!func_name %in% c("+", "*"))
    assert_valid_character(func_name, "function name")

  assert_valid_arg_keys(names(func_args))
  paste0(c(func_name,
           paste0(names(func_args), collapse = ":"),
           paste0(sprintf("[%s]", func_args), collapse = ":")),
         collapse = "|")
}

# Strip redundant parentheses.
strip_useless_parens <- function(x) {
  while (rlang::is_call(x) &&
           rlang::call_name(x) == "(" &&
           length(rlang::call_args(x)) == 1) {
    x <- rlang::call_args(x)[[1]]
  }
  x
}

# Collect leaf operands by flattening chains of the same binary operator.
collect_flat_args <- function(node, op) {
  node_args <- lapply(rlang::call_args(node), strip_useless_parens)
  # a - b folds as a + (-b)
  ops <- if (op == "+") c("+", "-") else op
  if (rlang::call_name(node) == "-")
    node_args[[2]] <- rlang::call2("-", node_args[[2]])

  result <- list()
  for (a in node_args) {
    if (rlang::is_call(a) && rlang::call_name(a) %in% ops &&
          length(rlang::call_args(a)) == 2) {
      result <- c(result, collect_flat_args(a, op))
    } else {
      result <- c(result, list(a))
    }
  }
  result
}

# Assert an expression node is an admissible definition.
assert_valid_defin_node <- function(x) {
  if (is.null(x) || rlang::is_syntactic_literal(x) ||
        rlang::is_scalar_atomic(x) || rlang::is_symbol(x) ||
        rlang::is_call(x)) return(invisible(x))
  stop_h(
    "Found a definition value that is not a valid source.",
    format_type_error(x)
  )
}

# Label a leaf operand as spelled, without quotes or backticks.
# A string and a name read the same, quotes and backticks being R syntax.
housepi_label_leaf <- function(y) {
  if (is.character(y)) y
  else if (rlang::is_syntactic_literal(y) && is.numeric(y)) as.character(y)
  else if (rlang::is_symbol(y)) rlang::as_string(y)
  else rlang::as_label(y)
}

# Encode a call as a function name over encoded arguments.
housepi_func_from_call <- function(call, env = rlang::global_env()) {
  args <- lapply(rlang::call_args(call), strip_useless_parens)

  # map to function name
  # `+` and `*` stay operators, apart from their reductions sum() and prod()
  func_name <- switch(rlang::call_name(call),
                      `+` = if (length(args) == 1) NA_character_ else "+",
                      `-` = if (length(args) == 1) "minus" else "+",
                      `/` = "div",
                      `^` = "pow",
                      rlang::call_name(call))

  # flatten operator chains to avoid nested calls
  if (func_name %in% c("+", "*"))
    args <- collect_flat_args(call, func_name)

  # named operands keep their key, unnamed ones fill the free slots x, x2, ...
  nms <- rlang::names2(args)
  assert_unique_arg_names(nms, .family$func, func_name)
  slots <- c("x", paste0("x", seq_len(2 * length(args))[-1]))
  nms[!nzchar(nms)] <- utils::head(setdiff(slots, nms), sum(!nzchar(nms)))

  # list operands encoded
  func_args <- character()
  for (i in seq_along(args)) {
    y <- coerce_date_(args[[i]])
    assert_valid_defin_node(y)
    is_subcall <- rlang::is_call(y)
    folded <- if (is_subcall) try_fold_const_(y, env) else NULL
    func_args[[nms[i]]] <-
      if (!is.null(folded)) housepi_format_sources(folded)
      else if (is_subcall) housepi_func_from_call(y, env)
      else housepi_format_sources(housepi_label_leaf(y))
  }
  # canonical order: slots first, then other names as written
  func_args <- func_args[order(match(nms, slots, nomatch = length(slots) + 1))]

  if (is.na(func_name)) unname(func_args)
  else if (!length(func_args)) func_name
  else housepi_format_call(func_name, func_args)
}

# Number of '|' separators in an encoded definition.
housepi_bar_count <- function(x) {
  if (is.null(x)) 0L
  else nchar(x) - nchar(gsub("|", "", x, fixed = TRUE))
}

# Whether `x` is a sourced value, one '|'.
housepi_is_sourced_value <- function(x) isTRUE(housepi_bar_count(x) == 1)

# Whether `x` is a call:
# * default function -> no '|';
# * call -> more than one.
housepi_is_call <- function(x) isTRUE(housepi_bar_count(x) > 1)

# Combine function name and arguments.
housepi_combine_to_call <- function(func_name, func_args, is_func = FALSE) {
  func_name <- as.character(func_name)
  if (housepi_is_call(func_name) && length(func_args)) {
    stop_h("Impossible to associate arguments to a function which is a call.")
  } else if (housepi_is_sourced_value(func_name) && isTRUE(is_func)) {
    func_name <- showcaser_extract(func_name, "?")
  }
  if (!length(func_args)) return(func_name)

  assert_valid_arg_keys(rlang::names2(func_args))
  sprintf("%s|%s|%s", func_name,
          paste(rlang::names2(func_args), collapse = ":"),
          paste(sprintf("[%s]", func_args), collapse = ":"))
}

# Unpack a list() or c() container into one quosure per element.
helper_case_container <- function(x_quos) {
  if (length(x_quos) == 1) {
    x_expr <- rlang::quo_get_expr(x_quos[[1]])
    x_env <- rlang::quo_get_env(x_quos[[1]])
    # @todo other containers than list() and c()
    if (rlang::is_call(x_expr, c("list", "c"))) {
      elements <- as.list(x_expr[-1])
      # convert each element to a quosure with the same environment
      x_quos <- lapply(elements, \(elem) rlang::new_quosure(elem, x_env))
    }
  }
  x_quos
}

# Unwrap a definition argument, evaluating it when defused.
housepi_smart_unwrap <- function(x, arg_name, null_ok = FALSE) {
  if (rlang::is_quosure(x) || rlang::is_quosures(x)) x <- rlang::eval_tidy(x)
  x <- coerce_date_(x)
  if (is.null(x) && isTRUE(null_ok)) return(NULL)
  if (is.character(x) || is.numeric(x) || is.logical(x) || is.list(x)) return(x)
  stop_h("Argument {.arg {arg_name}} must be character, numeric, or logical.",
         format_type_error(x))
}

# Encode the definitions carried by housepi quosures.
housepi_extractor <- function(...) {
  .extract_node <- function(x_quo) {
    x_env <- rlang::quo_get_env(x_quo)
    x_expr <- coerce_date_(rlang::quo_get_expr(x_quo))

    # expand a formula
    if (rlang::is_formula(x_expr)) x_expr <- rlang::f_rhs(x_expr)
    x_expr <- strip_useless_parens(x_expr)

    if (is.null(x_expr)) {
      NULL
    } else if (rlang::is_call(x_expr, c("-", "+")) &&
                 length(x_expr) == 2 &&
                 rlang::is_syntactic_literal(x_expr[[2]])) {
      # unary
      housepi_format_sources(c(x = eval(x_expr)))
    } else if (rlang::is_call(x_expr)) {
      # * constant call, e.g. as.Date() -> folded to its numeric value;
      # * call over housepi nodes -> fails to evaluate, stays symbolic.
      folded <- try_fold_const_(x_expr, x_env)
      if (!is.null(folded)) housepi_format_sources(c(x = folded))
      else housepi_func_from_call(x_expr, x_env)
    } else if (rlang::is_syntactic_literal(x_expr) ||
                 rlang::is_scalar_atomic(x_expr) ||
                 rlang::is_symbol(x_expr)) {
      if (rlang::is_symbol(x_expr)) x_expr <- rlang::as_string(x_expr)
      housepi_format_sources(c(x = x_expr))
    } else {
      assert_valid_defin_node(x_expr)
    }
  }

  x_quos <- rlang::enquos(...)
  x_quos <- helper_case_container(x_quos)
  output <- unlist(lapply(x_quos, .extract_node))

  if (!is.null(names(x_quos)) && any(nzchar(names(x_quos)))) output
  else unname(output)
}

# Load a named list of quosures for definitions.
# * not `append` -> store reset first;
# * `defaults` -> named arguments for every definition not setting them,
#   e.g. a bulk-level 'status =' shared by all entries.
housepi_load_defin_quos <- function(self, implicit_name, xquos,
                                    append = FALSE, defaults = list()) {
  if (isFALSE(append)) self$initialize()
  self$store_warn_over <- character()

  # nameless entry, outside the infectivity profile -> list passed without !!!
  requires_name <- "name" %in% names(formals(self$add_quo))

  # only inject defaults the per-definition arguments did not override
  inject_defaults <- \(keys) defaults[setdiff(names(defaults), keys)]

  for (i in seq_along(xquos)) {
    nm <- names(xquos[i])
    if (requires_name && (!length(nm) || is.na(nm) || nm == "")) stop_h(
      "Name of {self$family} must be a non-empty string.",
      "Did you provide a list and forgot rlang injection operator {.code !!!}?"
    )
    def_quo <- xquos[[i]]
    def_expr <- rlang::quo_get_expr(def_quo)
    # symbol or call other than a named list() / c() -> its value, if a list
    named <- rlang::is_call(def_expr, c("list", "c")) &&
      any(nzchar(rlang::names2(rlang::call_args(def_expr))))
    if (!named && rlang::is_symbolic(def_expr)) {
      val <- tryCatch(suppressWarnings(rlang::eval_tidy(def_quo)),
                      error = \(e) NULL)
      if (rlang::is_bare_list(val)) def_expr <- rlang::call2("list", !!!val)
    }

    if (rlang::is_call(def_expr, c("list", "c"))) {
      # pass arguments as quosures to add_quo
      args_quos <- lapply(rlang::call_args(def_expr),
                          rlang::new_quosure,
                          env = rlang::quo_get_env(def_quo))

      if (!any(nzchar(names(args_quos))))
        args_quos <- rlang::set_names(list(def_quo), implicit_name)

      inj <- inject_defaults(names(args_quos))
      rlang::inject(self$add_quo(name = nm, !!!args_quos, !!!inj,
                                 .warn_over = FALSE))

    } else {
      # single implicit definition, a stray 'name' dropped by a nameless holder
      extra_arg <- rlang::set_names(list(def_quo), implicit_name)
      inj <- inject_defaults(implicit_name)
      rlang::inject(
        self$add_quo(name = nm, !!!extra_arg, !!!inj, .warn_over = FALSE)
      )
    }
  }
  definition_overwrite_warn(self$store_warn_over, self$family)
  self$store_warn_over <- character()
}

# Convert a definition list to a list of quosures.
convert_list_to_quosures <- function(input_lst, env = rlang::global_env()) {
  lapply(input_lst, function(x) {
    if (is.language(x)) rlang::new_quosure(x, env = env)
    else rlang::new_quosure(rlang::call2("list", !!!x), env = env)
  })
}

# Load a templated list of definitions.
get_template_def_list <- function(path, family) {
  templ_file <- file.path(path, template_def_files[[family]])

  # check if template file exists
  if (!file.exists(templ_file)) stop_h(
    "Impossible to use templated {family} definitions.",
    `!` = "A file with template definitions was not found."
  )
  # get template definition list and convert to quosures
  dget(templ_file)
}
