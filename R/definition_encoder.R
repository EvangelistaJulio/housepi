## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

# Reserved constant names forbidden as definition names.
.reserved_math_consts <- c("pi")

# Coerce an R Date value to its numeric day count.
coerce_date_ <- function(x) {
  if (inherits(x, "Date")) as.numeric(unclass(x))
  else x
}

# Helper to check type definitions and set defaults.
helper_default_source <- function(types, values) {
  valid_sources_lst <- list("A" = "augmented",
                            "D" = "data",
                            "G" = "group rate",
                            "N" = "numeric",
                            "P" = "parameter")

  for (i in seq_along(types)) {
    types[i] <- toupper(types[i])
    k <- types[i]
    y <- values[[i]]

    # explicit type checking
    if (is.na(k)) {
      # no explicit type provided, then format for default types
      if (is.na(y)) types[i] <- NA_character_
      else if (is.numeric(y)) types[i] <- "N"
      else if (!is.character(y)) stop_h(
        "Found an unsupported definition type value.",
        format_type_error(y)
      )
      else if (!is.na(suppressWarnings(as.numeric(y)))) types[i] <- "N"
      else if (y %in% .reserved_math_consts) types[i] <- "N"
      else types[i] <- "A"

    } else if (!(k %in% names(valid_sources_lst))) {
      stop_h(paste0("Invalid type {.val ", k, "} for value {.val ", y, "}."),
             paste(sprintf("{.val %s}", names(valid_sources_lst)),
                   sprintf("(%s)", unlist(valid_sources_lst)),
                   collapse = ", "))
    }
  }
  types
}

# Add default types if not a call.
housepi_format_sources <- function(func_args) {
  # default return
  if (is.null(func_args) || all(is.na(func_args))) return(NA_character_)

  splitted <- strsplit(as.character(func_args), "\\@")
  sapply(splitted, function(x) {
    ifelse(stringr::str_count(x[1], "\\|"), x[1],
           sprintf("%s|%s", x[1], helper_default_source(x[2], x[1])))
  }, USE.NAMES = FALSE)
}

# Format a call to match housepi standard typed definitions
housepi_format_call <- function(func_name, func_args) {
  paste0(c(func_name,
           paste0(names(func_args), collapse = ":"),
           paste0(sprintf("{%s}", housepi_format_sources(func_args)),
                  collapse = ":")),
         collapse = "|")
}

# Helper to unwrap from useless parenthesis.
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
  result <- list()
  for (a in node_args) {
    if (rlang::is_call(a) && rlang::call_name(a) == op &&
          length(rlang::call_args(a)) == 2) {
      result <- c(result, collect_flat_args(a, op))
    } else {
      result <- c(result, list(a))
    }
  }
  result
}

# Assert that expression node is an admissible definition.
assert_valid_defin_node <- function(x) {
  if (is.null(x) || rlang::is_syntactic_literal(x) ||
        rlang::is_scalar_atomic(x) || rlang::is_symbol(x) ||
        rlang::is_call(x)) return(invisible(x))
  stop_h(
    "Found a definition value that is not a valid source.",
    format_type_error(x)
  )
}

# Label a leaf operand.
housepi_label_leaf <- function(y) {
  if (rlang::is_syntactic_literal(y) && is.numeric(y)) as.character(y)
  else rlang::as_label(y)
}

# Helper to extract functin from call.
housepi_func_from_call <- function(call) {
  args <- lapply(rlang::call_args(call), strip_useless_parens)

  # map to function name
  func_name <- switch(rlang::call_name(call),
                      `+` = if (length(args) == 1) NA_character_ else "sum",
                      `-` = if (length(args) == 1) "minus" else "sub",
                      `*` = "prod",
                      `/` = "div",
                      `^` = "pow",
                      rlang::call_name(call))

  # flatten variadic operators to avoid nested chains
  if (func_name %in% c("sum", "prod")) {
    args <- collect_flat_args(call, rlang::call_name(call))
  }

  # list operands as literals
  func_args <- character()
  if (length(args) > 0) {
    for (i in seq_along(args)) {
      y <- coerce_date_(args[[i]])
      assert_valid_defin_node(y)
      func_args[[paste0("x", i)]] <- ifelse(
        rlang::is_call(y) && rlang::call_name(y) != "@",
        housepi_func_from_call(y),
        housepi_label_leaf(y)
      )
    }
    names(func_args)[1] <- "x"
  }

  if (is.na(func_name)) housepi_format_sources(func_args)
  else if (length(func_args) == 0) func_name
  else housepi_format_call(func_name, func_args)
}

# Flag to tell if object is a call or a value.
# NOTE: default function are 0 '|'
housepi_is_sourced_value <- function(x) {
  if (is.null(x)) FALSE
  else isTRUE(stringr::str_count(x, "\\|") == 1)
}

# Flag to tell if object is a call or a value.
# NOTE: default function are 0 '|'
housepi_is_call <- function(x) {
  if (is.null(x)) FALSE
  else isTRUE(stringr::str_count(x, "\\|") > 1)
}

# Combione function name and arguments.
housepi_combine_to_call <- function(func_name, func_args, is_func = FALSE) {
  func_name <- as.character(func_name)
  if (housepi_is_call(func_name) && length(func_args) > 0) {
    stop_h("Impossible to associate arguments to a function which is a call.")
  } else if (housepi_is_sourced_value(func_name) && isTRUE(is_func)) {
    func_name <- showcaser_extract(func_name, "A")
  }
  if (length(func_args) == 0) return(func_name)

  sprintf("%s|%s|%s", func_name,
          paste(names(func_args), collapse = ":"),
          paste(sprintf("{%s}", func_args), collapse = ":"))
}

# Handle case of list() or c() container.
helper_case_container <- function(x_quos) {
  if (length(x_quos) == 1) {
    x_expr <- rlang::quo_get_expr(x_quos[[1]])
    x_env <- rlang::quo_get_env(x_quos[[1]])
    # inspect for container # UPDATE is this complete?
    if (rlang::is_call(x_expr, c("list", "c"))) {
      elements <- as.list(x_expr[-1])
      # convert each element to a quosure with the same environment
      x_quos <- lapply(elements, function(elem) {
        rlang::new_quosure(elem, x_env)
      })
    }
  }
  x_quos
}

# Deal with arguments not expected to be quosure.
housepi_smart_unwrap <- function(x, arg_name = NULL) {
  # If not a quosure, return as-is
  if (!rlang::is_quosure(x) && !rlang::is_quosures(x)) return(x)

  # evaluate the quosure
  result <- coerce_date_(rlang::eval_tidy(x))

  # check if is one of the allowed types
  if (is.character(result) || is.numeric(result) ||
        is.logical(result) || is.list(result)) return(result)

  if (length(result) == 0) return(NULL)

  tryCatch({
    lapply(result, rlang::quo_get_expr)
  },
  error = function(e) {
    stop_h(
      sprintf("%s must be character, numeric, or logical.",
              ifelse(is.null(arg_name), "Defused arguments",
                     sprintf("Argument {.arg %s}", arg_name))),
      format_type_error(result)
    )
  })
}

# Deal with housepi quosures.
housepi_extractor <- function(...) {
  .extract_node <- function(x_quo) {
    x_expr <- coerce_date_(rlang::quo_get_expr(x_quo))

    # expand formula in case
    if (rlang::is_formula(x_expr)) x_expr <- rlang::f_rhs(x_expr)
    x_expr <- strip_useless_parens(x_expr)

    if (is.null(x_expr)) {
      NULL
    } else if (rlang::is_call(x_expr, c("-", "+")) &&
                 length(x_expr) == 2 &&
                 rlang::is_syntactic_literal(x_expr[[2]])) {
      # unary
      housepi_format_sources(c(x = eval(x_expr)))
    } else if (rlang::is_call(x_expr, "@")) {
      # sourced unary
      housepi_format_sources(c(x = rlang::expr_text(x_expr)))
    } else if (rlang::is_call(x_expr)) {
      # function
      housepi_func_from_call(x_expr)
    } else if (rlang::is_syntactic_literal(x_expr) ||
                 rlang::is_scalar_atomic(x_expr) ||
                 rlang::is_symbol(x_expr)) {
      housepi_format_sources(c(x = x_expr))
    } else {
      assert_valid_defin_node(x_expr)
    }
  }

  x_quos <- rlang::enquos(...)
  x_quos <- helper_case_container(x_quos)
  output <- unlist(lapply(x_quos, .extract_node))

  if (!is.null(names(x_quos)) && !all(names(x_quos) == "")) output
  else unname(output)
}

# Detect whether captured definitions declare an observation node ('obs =').
# Observation nodes are reserved to add_observation();
# this guards the augmented entry points (single and bulk forms).
housepi_quos_declare_obs <- function(xquos) {
  # single form: add_augmented(name, obs = ...)
  if ("obs" %in% names(xquos)) return(TRUE)
  # bulk form: add_augmented(name = list(obs = ...))
  any(vapply(xquos, function(q) {
    expr <- if (rlang::is_quosure(q)) rlang::quo_get_expr(q) else q
    if (rlang::is_call(expr, c("list", "c"))) {
      "obs" %in% names(rlang::call_args(expr))
    } else if (is.list(expr)) {
      "obs" %in% names(expr)
    } else {
      FALSE
    }
  }, logical(1)))
}

# Stop if augmented definitions try to declare an observation node.
assert_no_obs_declaration <- function(xquos) {
  if (housepi_quos_declare_obs(xquos)) stop_h(
    "Augmented definitions cannot declare an observation node.",
    "Did you mean {.fn add_observation} instead?"
  )
}

# Stop if no definition was supplied to the list.
assert_nonempty_defin_quos <- function(xquos, family) {
  if (length(xquos) == 0L) stop_h(
    "No {family} definition provided.",
    "Provide at least one definition."
  )
}

# Load a named list of quosures for definitions.
# Resets the store first unless 'append' is TRUE.
# 'defaults' supplies named arguments applied to every definition that does not
# set them itself (e.g. a bulk-level 'status =' shared by all entries).
housepi_load_defin_quos <- function(self, implicit_name, xquos,
                                    append = FALSE, defaults = list()) {
  assert_nonempty_defin_quos(xquos, self$family)
  if (isFALSE(append)) self$initialize()
  self$store_warn_over <- character()

  # all holder but the infectivity profile key their definitions by a name; a
  # nameless entry there is the tell-tale of a list passed without '!!!'
  requires_name <- "name" %in% names(formals(self$add_quo))

  # only inject defaults the per-definition arguments did not override
  inject_defaults <- function(keys) defaults[setdiff(names(defaults), keys)]

  for (i in seq_along(xquos)) {
    nm <- names(xquos[i])
    if (requires_name && (is.na(nm) || length(nm) == 0 || nm == "")) stop_h(
      "Name of {self$family} must be a non-empty string.",
      "Did you provide a list and forgot rlang injection operator {.code !!!}?"
    )
    def_quo <- xquos[[i]]
    def_expr <- rlang::quo_get_expr(def_quo)

    if (rlang::is_call(def_expr, c("list", "c"))) {
      # pass arguments as quosures to add_quo
      args_quos <- lapply(rlang::call_args(def_expr),
                          rlang::new_quosure,
                          env = rlang::quo_get_env(def_quo))

      if (all(names(args_quos) == "")) {
        args_quos <- rlang::set_names(list(def_quo), implicit_name)
      }
      inj <- inject_defaults(names(args_quos))
      rlang::inject(self$add_quo(name = nm, !!!args_quos, !!!inj,
                                 .warn_over = FALSE))

    } else if (is.list(def_expr)) {
      # arguments are not quosures
      inj <- inject_defaults(names(def_expr))
      rlang::inject(self$add_quo(name = nm, !!!def_expr, !!!inj,
                                 .warn_over = FALSE))
    } else {
      # single implicit definition to pass to add_quo. A stray 'name' is
      # harmless for the nameless infectivity profile (dropped by its add_quo).
      extra_arg <- rlang::set_names(list(def_quo), implicit_name)
      inj <- inject_defaults(implicit_name)
      rlang::inject(self$add_quo(name = nm, !!!extra_arg, !!!inj,
                                 .warn_over = FALSE))
    }
  }
  template_overwrite_warn(self$store_warn_over, self$family)
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
  if (family == "parameter") {
    templ_file <- file.path(path, "param.R")
  } else if (family == "augmented") {
    templ_file <- file.path(path, "augm.R")
  }

  # check if template file exists
  if (!file.exists(templ_file)) stop_h(
    "Impossible to use templated {family} definitions.",
    "A file with template definitions was not found."
  )
  # get template definition list and convert to quosures
  dget(templ_file)
}
