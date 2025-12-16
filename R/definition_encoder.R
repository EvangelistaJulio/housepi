#
#        __
#       / /_  ____  __  __ ____ __  ____  __
#      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#     / / / / /_/ / /_/ (__  )  __/ /_/ / /
#    /_/ /_/\____/\____/____/\___/ ____/_/
#                               /_/
#
# housepi - Epidemiological Analysis of Household Outbreaks.
#
# Author: Francesco Camaglia
# Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
#
# Distributed under the terms of the GNU General Public License (GPLv3).
# See the LICENSE file for details.
#
# definition_encoder.R
#

# Helper to check type definitions and set defaults.
helper_default_source <- function(types, values) {
  valid_sources_lst <- list("A" = "augmented",
                            "D" = "data",
                            "G" = "group rate",
                            "N" = "numeric",
                            "P" = "parameter")

  for (i in seq_along(types)) {
    k <- types[i]
    y <- values[[i]]

    # explicit type checking
    if (is.na(k)) {
      # no explicit type provided, then format for default types
      if (is.na(y)) types[i] <- NA_character_
      else if (is.numeric(y)) types[i] <- "N"
      else if (is.character(y)) types[i] <- ifelse(
        is.na(suppressWarnings(as.numeric(y))), "A", "N"
      )
      else types[i] <- "A"
    } else if (!(toupper(k) %in% names(valid_sources_lst))) {
      stop_h(
        paste0("Invalid type '", k, "' for value '", y, "'."),
        paste(sprintf("'%s'", names(valid_sources_lst)),
              sprintf("(%s)", unlist(valid_sources_lst)),
              collapse = ", ")
      )
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

# Helper to extract functin from call.
housepi_func_from_call <- function(call) {

  # Helper to unwrap from useless parenthesis.
  strip_useless_parens <- function(x) {
    while (rlang::is_call(x) &&
             rlang::call_name(x) == "(" &&
             length(rlang::call_args(x)) == 1) {
      inner_arg <- rlang::call_args(x)[[1]]
      if (rlang::is_symbol(inner_arg)) x <- inner_arg
      else if (rlang::call_name(inner_arg) == "@") x <- inner_arg
      else break
    }
    x
  }
  args <- lapply(rlang::call_args(call), strip_useless_parens)

  # list operands as literals
  func_args <- character()
  if (length(args) > 0) {
    for (i in seq_along(args)) {
      y <- args[[i]]
      func_args[[paste0("x", i)]] <- ifelse(
        rlang::is_call(y) && rlang::call_name(y) != "@",
        housepi_func_from_call(y),
        rlang::as_label(y)
      )
    }
    names(func_args)[1] <- "x"
  }

  # map to function name
  func_name <- switch(rlang::call_name(call),
                      `+` = if (length(args) == 1) NA_character_ else "sum",
                      `-` = if (length(args) == 1) "minus" else "sub",
                      `*` = "prod",
                      `/` = "div",
                      `^` = "pow",
                      rlang::call_name(call))

  if (is.na(func_name)) housepi_format_sources(func_args)
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
housepi_combine_to_call <- function(func_name, x) {
  if (housepi_is_call(func_name)) {
    stop_h("Impossible to associate arguments to a function which is a call.")
  } else if (housepi_is_sourced_value((func_name))) {
    # verify it's not a default sourced value mistakenly interpreted
    func_name <- showcaser_extract(func_name, "A")
    if (housepi_is_sourced_value(func_name)) stop_h(
      "Impossible to associate arguments to a sourced value."
    )
  }
  if (length(x) == 0) return(func_name)

  sprintf("%s|%s|%s", func_name,
          paste(names(x), collapse = ":"),
          paste(sprintf("{%s}", x), collapse = ":"))
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
  result <- rlang::eval_tidy(x)

  # check if is one of the allowed types
  if (is.character(result) || is.numeric(result) ||
        is.logical(result)) return(result)

  if (length(result) == 0) return(NULL)

  tryCatch({
    lapply(result, rlang::quo_get_expr)
  },
  error = function(e) {
    stop_h(sprintf(
      "%s must be character, numeric, or logical, not %s.",
      ifelse(is.null(arg_name), "Defused arguments",
             sprintf("Argument '%s'", arg_name)),
      class(result)[1]
    ))
  })
}

# Deal with housepi quosures.
housepi_extractor <- function(...) {
  x_quos <- rlang::enquos(...)
  x_quos <- helper_case_container(x_quos)

  output <- unlist(lapply(x_quos, function(x_quo) {
    x_expr <- rlang::quo_get_expr(x_quo)
    # expand formula in case
    if (rlang::is_formula(x_expr)) x_expr <- rlang::f_rhs(x_expr)

    if (is.null(x_expr)) {
      NULL
    } else if (rlang::is_call(x_expr)) {
      # apply housepi format to formulas
      if (rlang::call_name(x_expr) %in% c("-", "+") &&
            length(x_expr) == 2 &&
            rlang::is_syntactic_literal(x_expr[[2]])) {
        # unary
        housepi_format_sources(c(x = eval(x_expr)))
      } else if (rlang::call_name(x_expr) == "@") {
        # sourced unary
        housepi_format_sources(c(x = rlang::expr_text(x_expr)))
      } else {
        # function
        housepi_func_from_call(x_expr)
      }
    } else if (rlang::is_syntactic_literal(x_expr) ||
                 rlang::is_symbol(x_expr)) {
      # literals
      housepi_format_sources(c(x = x_expr))
    } else {
      NULL
    }
  }))

  if (!is.null(names(x_quos)) && !all(names(x_quos) == "")) output
  else unname(output)
}

# Load a named list of quosures for definitions.
housepi_load_defin_quos <- function(self, implicit_name, xquos) {
  for (i in seq_along(xquos)) {
    nm <- names(xquos[i])
    def_quo <- xquos[[i]]
    def_expr <- rlang::quo_get_expr(def_quo)

    if (rlang::is_call(def_expr, c("list", "c"))) {
      # pass arguments as quosures to add_quo
      args_quos <- lapply(rlang::call_args(def_expr),
                          rlang::new_quosure,
                          env = rlang::quo_get_env(def_quo))
      rlang::inject(self$add_quo(nm, !!!args_quos))
    } else if (is.list(def_expr)) {
      # arguments are not quosures
      rlang::inject(self$add_quo(nm, !!!def_expr))
    } else {
      # single implicit definition to pass to add_quo
      extra_arg <- rlang::set_names(list(def_quo), implicit_name)
      # associate name if provided
      if (is.na(nm) || length(nm) == 0 || nm == "") {
        rlang::inject(self$add_quo(!!!extra_arg))
      } else {
        rlang::inject(self$add_quo(name = nm, !!!extra_arg))
      }
    }
  }
  self
}

# Convert a definition list to a list of quosures.
convert_list_to_quosures <- function(input_lst, env = rlang::global_env()) {
  lapply(input_lst, function(x) {
    rlang::new_quosure(rlang::call2("list", !!!x), env = env)
  })
}

# Load a templated list of definitions.
housepi_load_template_list <- function(self, templ_file, family,
                                       required = NULL, append = FALSE) {
  # check if template file exists
  if (!file.exists(templ_file)) stop_h(
    sprintf("Impossible to use templated %s definition.", family),
    "A file with template definitions was not found."
  )
  # get template definition list and convert to quosures
  xquos <- convert_list_to_quosures(dget(templ_file))

  if (isTRUE(append)) {
    # check set of required augmented variables
    if (length(required) == 0) {
      warning_h(sprintf("No %s definitions required by the model.", family))
    } else {
      miss <- template_missing_find(required, self$defs[["name"]])
      if (length(miss) == 0) warning_h(
        sprintf("No template %s variable to append.", family),
        "Definitions are complete."
      )
      else if (!all(miss %in% names(xquos))) stop_h(sprintf(
        "No template found for the missing augmented definitionas.", family
      ))
      else self$load_quosures(xquos[miss], append = TRUE)
    }
  } else {
    self$load_quosures(xquos, append = FALSE)
  }
  self
}

# Assert that arguments provided to model are legitimate.
# UPDATE: to be generalized beyond distributions
housepi_assert_args <- function(hhmodel, defs, family, max_print = 5) {
  readable_defs <- housepi_defin_showcaser(defs)
  for (r in seq_len(nrow(readable_defs))) {
    dist_name <- readable_defs[r, family]
    if (!is.na(dist_name)) {
      col_args <- ifelse(family == "func", "func_args", "args")
      arg_names <- names(readable_defs[r, col_args][[1]])
      # NOTE: I am assuming `1` is a safe argument...
      fake_args <- stats::setNames(rep(1, length(arg_names)), arg_names)
      # trigger C++ errors
      tryCatch(
        expose_distribution(dist_name, fake_args, hhmodel = hhmodel),
        error = function(e) {
          stop_h(sprintf("Wrong definition found for '%s'.",
                         readable_defs[r, "name"]),
                 strsplit(e$message, "\n")[[1]][-1])
        }
      )
    }
  }
}

# Helper to load definition lists into an HousepiModel clone.
helper_load_def_list <- function(pvt_module, hhmodel, xquos, append) {
  stopifnot(inherits(hhmodel, "HousepiModel"))
  new_object <- hhmodel$clone(deep = TRUE)
  new_object$.__enclos_env__$private[[pvt_module]]$load_quosures(
    xquos, append = append
  )
  housepi_assert_args(
    new_object,
    new_object$.__enclos_env__$private[[pvt_module]]$defs,
    ifelse(pvt_module == ".augm", "model", "prior")
  )
  new_object
}
