#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# definition_encoder.R
#

#' @title Definition Encoder
#' @description
#' Helper with definitions.
#'
#' @keywords internal
#' @noRd
DefinitionEncoder <- R6::R6Class(
  classname = "DefinitionEncoder",

  private = list(
    # If it's a call to "(" with a single argument that's a symbol, unwrap it.
    .strip_useless_parens = function(x) {
      while (rlang::is_call(x) &&
               rlang::call_name(x) == "(" &&
               length(rlang::call_args(x)) == 1) {
        inner_arg <- rlang::call_args(x)[[1]]
        if (rlang::is_symbol(inner_arg)) x <- inner_arg
        else if (rlang::call_name(inner_arg) == "@") x <- inner_arg
        else break
      }
      x
    },

    # Check type definitions and set defaults.
    .default_types = function(types, values) {
      valid_types <- list("A" = "augmented",
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
        } else if (!(toupper(k) %in% names(valid_types))) {
          stop_h(
            paste0("Invalid type '", k, "' for value '", y, "'."),
            paste(sprintf("'%s'", names(valid_types)),
                  sprintf("(%s)", unlist(valid_types)),
                  collapse = ", ")
          )
        }
      }

      types
    },

    # Add default types if not a call.
    .format_types = function(func_args) {
      output <- NA_character_
      if (!(is.null(func_args) || all(is.na(func_args)))) {
        splitted <- strsplit(as.character(func_args), "\\@")
        output <- sapply(splitted, function(x) {
          ifelse(stringr::str_count(x[1], "\\|"), x[1],
                 sprintf("%s|%s", x[1], private$.default_types(x[2], x[1])))
        }, USE.NAMES = FALSE)
      }
      output
    },

    # Format a call to match housepi standard typed definitions
    .fmt_identity = function(func_args) {
      if (is.null(func_args)) NA_character_
      else private$.format_types(func_args)
    },

    # Format a call to match housepi standard typed definitions
    .format_call = function(func_name, func_args) {
      paste0(c(
        func_name,
        paste0(names(func_args), collapse = ":"),
        paste0(sprintf("{%s}", private$.format_types(func_args)),
               collapse = ":")
      ), collapse = "|")
    },

    # Helper to extract functin from call.
    .helper_func_from_call = function(call) {
      # strip meaningless parentheses from arguments
      args <- lapply(rlang::call_args(call), private$.strip_useless_parens)

      # list operands as literals
      func_args <- character()
      if (length(args) > 0) {
        for (i in seq_along(args)) {
          y <- args[[i]]
          func_args[[paste0("x", i)]] <- ifelse(
            rlang::is_call(y) && rlang::call_name(y) != "@",
            private$.helper_func_from_call(y),
            rlang::as_label(y)
          )
        }
        names(func_args)[1] <- "x"
      }

      # map to function name
      func_name <- switch(
        rlang::call_name(call),
        `+` = if (length(args) == 1) NA_character_ else "sum",
        `-` = if (length(args) == 1) "minus" else "sub",
        `*` = "prod",
        `/` = "div",
        `^` = "pow",
        rlang::call_name(call) # default
      )

      if (is.na(func_name)) private$.fmt_identity(func_args)
      else private$.format_call(func_name, func_args)
    },

    # Extract function or identiy call from rlang expression.
    .extract = function(expr) {
      # read input expression
      if (is.null(expr)) return(private$.fmt_identity(NULL))
      else if (inherits(expr, "formula")) call <- rlang::f_rhs(expr)
      else if (rlang::is_call(expr)) call <- expr
      else return(private$.fmt_identity(c(x = expr)))

      # case of a symbol actually looking as a call
      if (rlang::is_symbol(call)) private$.fmt_identity(c(x = call))
      # case of a raw typed identity
      else if (rlang::call_name(call) == "@") private$.fmt_identity(c(x = call))
      # standard call
      else private$.helper_func_from_call(call)
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field x `character(1)`\cr
    #' Encoding after initialization.
    x = NA_character_,

    #' @description Flag to tell if object is a call or a value.
    is_call = function() {
      # NOTE: default function are 0 '|'
      isTRUE(stringr::str_count(self$x, "\\|") > 1)
    },

    #' @description Return the value, if applicable.
    #' @param func_name
    #' A `character(1)` with the function name.
    #'
    #' @return `character(1)` or `NA`
    combine_to_call = function(func_name) {
      if (DefinitionEncoder$new(func_name)$is_call()) stop_h(
        "Function name is actually a call, impossible to combine."
      )
      sprintf("%s|%s|%s", func_name,
              paste(names(self$x), collapse = ":"),
              paste(sprintf("{%s}", self$x), collapse = ":"))
    },

    #' @description
    #' New object of class
    #' [HousepiAugmented][housepi::HousepiAugmented].
    initialize = function(expr) {
      if (is.list(expr)) self$x <- sapply(expr, private$.extract)
      else self$x <- private$.extract(expr)
      invisible()
    }
  )
)

############################################################################
############################################################################
############################################################################

#' @title Definition Walker
#' @description
#' Helper with type extraction from definitions.
#'
#' @keywords internal
#' @noRd
DefinitionWalker <- R6::R6Class(
  classname = "DefinitionWalker",
  private = list(
    # Helper to extract balanced {something} blocks at the top level.
    .extract_blocks = function(s) {
      blocks <- list()
      depth <- 0
      start <- NULL
      for (i in seq_along(strsplit(s, "")[[1]])) {
        ch <- substr(s, i, i)
        if (ch == "{") {
          if (depth == 0) start <- i
          depth <- depth + 1
        } else if (ch == "}") {
          depth <- depth - 1
          if (depth == 0 && !is.null(start)) {
            blocks <- c(blocks, substr(s, start, i))
            start <- NULL
          }
        }
      }
      blocks
    },

    # Recursive helper function to evaluate blocks.
    .evaluate_block = function(block, type) {
      # Remove curly braces
      block <- sub("^\\{", "", block)
      block <- sub("\\}$", "", block)

      # split the block into individual conditions
      conditions <- unlist(strsplit(block, ":"))
      output <- c()

      # Evaluate each condition
      for (cond in conditions) {
        if (grepl("\\{", cond)) {
          # recurse on nested blocks
          output <- c(output, DefinitionWalker$new(cond, type)$x)
        } else {
          pts <- unlist(strsplit(cond, "\\|"))
          if (length(pts) == 2 && pts[2] == type) output <- c(output, pts[1])
        }
      }
      output
    }
  ),

  public = list(
    #' @field x `character()`\cr
    #' Encoding after initialization.
    x = NULL,

    #' @param input
    #' A `character(1)` with encoded function string.
    #' @param type
    #' A `character(1)` with the type to extract.
    #'
    #' @return `character()`
    initialize = function(input, type) {
      self$x <- NULL
      if (!is.na(input) && !is.null(input) && input != "") {
        blocks <- private$.extract_blocks(input)
        if (length(blocks) > 0) {
          output <- c()
          for (b in blocks) output <- c(output, private$.evaluate_block(b, type))
          if (length(output) > 0) self$x <- unique(output)
        } else {
          # terminal node
          parts <- unlist(strsplit(input, "\\|"))
          if (length(parts) == 2 && parts[2] == type) self$x <- parts[1]
        }
      }
      invisible(self)
    }
  )
)

############################################################################
############################################################################
############################################################################

#' @title Definition Showcaser
#' @description
#' Helper with showcasing definitions in an accessible way.
#'
#' @keywords internal
#' @noRd
DefinitionShowcaser <-  R6::R6Class(
  classname = "DefinitionShowcaser",

  private = list(
    # helper to extract numeric values from a definition
    .extract_numeric = function(x) {
      sapply(x, function(element) {
        if (is.na(element)) return(NA_real_)
        if (!is.character(element)) return(element)
        z <- strsplit(element, "\\|")[[1]]
        if (length(z) != 2 || z[2] != "N") return(element)
        as.numeric(z[1])
      }, USE.NAMES = FALSE)
    },

    # helper to extract function name from a definition
    .extract_func = function(x) {
      sapply(x, function(element) {
        if (is.na(element)) return(NA_real_)
        if (!is.character(element)) return(element)
        z <- strsplit(element, "\\|")[[1]]
        if (length(z) <= 2) return(element)
        z[1]
      }, USE.NAMES = FALSE)
    },

    # helper to split most external semicolon separator
    .split_external = function(x) {
      chars <- strsplit(x, "")[[1]]
      depth <- 0
      split_pos <- c(0)  # start of string
      for (i in seq_along(chars)) {
        if (chars[i] == "{") depth <- depth + 1
        else if (chars[i] == "}") depth <- depth - 1
        else if (chars[i] == ":" && depth == 0) split_pos <- c(split_pos, i)
      }
      split_pos <- c(split_pos, nchar(x) + 1)

      # extract sub-strings between split positions
      result <- character(length(split_pos) - 1)
      for (i in seq_along(result)) {
        start <- split_pos[i] + 1
        end <- split_pos[i + 1] - 1
        result[i] <- substr(x, start, end)
      }

      # remove outer braces at the edges
      gsub("^\\{(.*)\\}$", "\\1", result)
    },

    # helper to extract arguments from a definition
    .extract_args = function(x) {
      lapply(x, function(element) {
        if (is.na(element)) return(NA_real_)
        if (!is.character(element)) return(c(element, NA_real_))
        z <- strsplit(element, "\\|")[[1]]
        if (length(z) <= 2) return(NA_real_)
        unlist(stats::setNames(
          private$.extract_numeric(
            private$.split_external(paste(z[3:length(z)], collapse = "|"))
          ),
          private$.split_external(z[2])
        ))
      })
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  list(
    #' @field x `data.frame()`\cr
    #' Readable definitions.
    x = NULL,

    #' @description
    #' Convert encoded definition strings into a more human-readable format.
    #'
    #' @param defs
    #' A `data.frame()` with encoded definition strings.
    #'
    #' @return `data.frame()`
    initialize = function(defs) {
      defs <- as.data.frame(defs)
      output <- as.data.frame(lapply(defs, private$.extract_numeric))
      for (col in colnames(defs)) {
        # find function definition
        if (isTRUE(any(stringr::str_count(output[[col]], "\\|") > 1))) {
          pos <- match(col, names(output))
          # format column with arguments
          # UPDATE: this is not clean at all.
          col_args <- data.frame(
            args = base::I(private$.extract_args(output[[col]]))
          )
          class(col_args$args) <- "list"
          if (col == "func")  col_args <- stats::setNames(col_args, "func_args")

          right_side_table <- cbind(
            stats::setNames(list(private$.extract_func(output[[col]])), col),
            col_args,
            output[(pos + 1):ncol(output)]
          )
          if (pos > 1) output <- cbind(output[1:(pos - 1)], right_side_table)
          else output <- right_side_table
        }
      }
      self$x <- output
      invisible(self)
    }
  )
)
