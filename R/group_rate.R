#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# group_rate.R
#

#' @title
#' Assert dataset features
#' @description
#' Assert that features are present in dataset,
#'
#' @param features
#' Features as list of dataset column names
#' @param data_headers
#' Headers of the dataset.
#'
#' @noRd
assert_dataset_features <- function(features,
                                    data_headers) {
  # assert that features are present
  missing <- c()
  for (this in features) {
    if (!(this %in% data_headers)) {
      missing <- c(missing, this)
    }
  }
  if (length(missing) > 0) {
    oss <- sprintf(
      "Unknown feature%s: \"%s\"",
      ifelse(length(missing) > 1, "s", ""),
      missing[[1]]
    )
    for (this in missing[-1]) {
      oss <- sprintf("%s, \"%s\"", oss, this)
    }
    stop_h(
      paste0(oss, "."),
      "Check your features and dataset headers."
    )
  }
}

#' @title
#' Group rate logic table
#' @description
#' Return logic table of group rates
#'
#' @param dataset
#' the dataset from which the etract group rates
#' @param by
#' Features as list of dataset column names
#'
#' @return
#' data.frame with logic table
#'
#' @noRd
group_rate_logic_table <- function(dataset, by) {
  features <- as.character(by)

  # convert features
  assert_dataset_features(features, colnames(dataset))

  # extract unique features and sort by list length
  uni_features <- lapply(features, function(f) sort(unique(dataset[[f]])))
  sorted_indices <- order(sapply(uni_features, length), decreasing = TRUE)
  features <- features[sorted_indices]
  uni_features <- uni_features[sorted_indices]
  names(uni_features) <- features

  # create a dataframe of all possible multi-features
  logic_table <- expand.grid(uni_features)

  # check that the multi-feature exists
  multifeat_count <- apply(logic_table, 1, function(row) {
    sum(apply(dataset[features], 1, function(d_row) all(d_row == row)))
  })

  # create the final truth table
  logic_table <- logic_table[multifeat_count > 0, , drop = FALSE]

  # reverse table to ease transcript
  logic_table <- rev(logic_table)

  # add counts
  logic_table$.counts <- multifeat_count[multifeat_count > 0]

  return(logic_table)
}

#' @title
#' Associate parameter names
#' @description
#' Associate parameter name to group rate logic table.
#'
#' @param this_table
#' Data.frame containing conditions for the group rate.
#' @param prefix
#' Prefix to parameter names associated to the group rate.
#' @param reference
#' Choose reference providing corresponding features.
#' Default is the group with the largest sample size.
#'
#' @return
#' input data.frame with logic table with parameter names.
#'
#' @noRd
associate_param_names <- function(this_table,
                                  prefix = NA,
                                  reference = NULL) {
  # column names to consider
  colmask <- setdiff(names(this_table),
                     c(".counts", ".param.name"))

  # extract reference index
  if (length(reference) == 0) {
    idx_refer <- which.max(this_table$.counts)

  } else {
    # find row
    idx_refer <- unname(which(
      apply(this_table[colmask], 1,
            function(row) all(row == unlist(reference[colmask])))
    ))

    if (length(idx_refer) != 1) {
      stop_h("Provided reference group is not valid.")
    }
  }

  # format prefix and associate parameter names
  prefix <- ifelse((!is.na(prefix)) && is.character(prefix),
                   paste0(prefix, "_"), "")

  # associate parameter names
  this_table$.param.name <- apply(
    this_table[colmask], 1, function(row) {
      paste0(prefix, paste(paste0(colmask, row), collapse = "_"))
    }
  )
  this_table$.param.name[idx_refer] <- "REFERENCE"

  return(this_table)
}

#' @description
#' Construct the C++ group rate emplacement.
#'
#' @param x `character(1)`\cr
#' Parameter associated to current group rate.
#' @param tab `logical(1)`\cr
#' If `TRUE` prints the tab (default).
#' @param quote `logical(1)`\cr
#' If `TRUE` prints the quotes (default).
#'
#' @return `character(1)`\cr
#'
#' @noRd
cpp_emplace <- function(x, tab = TRUE, quote = TRUE) {
  sprintf("%semplace_group_rate(r, %s);",
          ifelse(isTRUE(tab), "\t", ""),
          ifelse(isTRUE(quote), sprintf("\"%s\"", x), x))
}

#' @description
#' Construct the C++ if clauses plus emplacement.
#'
#' @param r
#' row of the group rate
#' @param is_ref
#' if `TRUE` the if clause is "REFERENCE"
#'
#' @noRd
cpp_ifclause_emplace <- function(r, is_ref = FALSE) {
  mask <- setdiff(names(r), c(".counts", ".param.name"))
  # if condition
  condition <- paste(
    paste0(
      sprintf("house.data(i, \"%s\") == ", names(r[mask])),
      r[mask]
    ),
    collapse = ") && ("
  )
  # wrap all together
  c(sprintf("%sif ((%s))", ifelse(is_ref, "", "else "), condition),
    ifelse(is_ref, 
           cpp_emplace("REFERENCE", tab = TRUE, quote = FALSE),
           cpp_emplace(r[[".param.name"]], tab = TRUE, quote = TRUE))
  )
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#' @title
#' Housepi Group Rate R6 Class
#' @description
#' Help with group rate definitions of model Rcpp object.
#'
#' @noRd
HousepiGroupRate <- R6::R6Class(
  classname = "HousepiGroupRate",
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################
  ############################################################################

  private = list(
    .rate_script = NULL,

    #' @description
    #' Generate a script to be printed within source code model.cpp
    #'
    #' @param dataset
    #' Dataset to be considered.
    .cpp_implement = function(dataset) {

      # initialize lines of the script
      rate_script <- c()
      self$logic_tables <- list()

      for (i in seq_along(self$defs)){
        rate_id <- names(self$defs[i]) # name of the rate
        this_def <- self$defs[[i]] # definition of the rate
        new_logic_table <- list()
        features <- as.list(this_def[["by"]])
        if (all(is.na(features))) {
          this_script <- c(cpp_emplace("REFERENCE", tab = FALSE, quote = FALSE))

        } else {
          features <- features[!is.na(features)]

          # initialize group rate contribution to script
          this_script <- c()

          if (isFALSE(this_def[["independent"]])) {
            features <- list(features)
          }

          # default prefix
          prefix <- as.character(this_def[["prefix"]])
          if (is.na(prefix)) {
            prefix <- format_string(tolower(rate_id))
            prefix <- substring(prefix, 1, 3)
          }

          for (this_feat in features) {
            # create logic table with group rates
            this_table <- group_rate_logic_table(dataset, this_feat)
            # associate names for parameters
            this_table <- associate_param_names(
              this_table, prefix,
              this_def[["reference"]]
            )

            # mask reference case
            mask <- this_table$.param.name == "REFERENCE"
            # concatenate new lines to the script
            this_script <- c(
              this_script,
              cpp_ifclause_emplace(this_table[mask, ], TRUE), # reference
              apply(this_table[!mask, ], 1, cpp_ifclause_emplace), # others
              "else",
              cpp_emplace("MISSING_RATE", tab = TRUE, quote = FALSE)
            )
            new_logic_table[[length(new_logic_table) + 1]] <- this_table
          }
        }

        # append the new rate chunk of script
        rate_script <- c(rate_script,
                         "",
                         sprintf("\t//----- %s -----//", rate_id),
                         sprintf("\tr = \"%s\";", rate_id),
                         sapply(this_script, function(x) paste0("\t", x)))
        # store the new rate logic table
        self$logic_tables[[rate_id]] <- new_logic_table
      }

      private$.rate_script <- private$.cpp_full_script(rate_script)
      invisible(self)
    },

    #' @description
    #' construct the C++ if clauses.
    #'
    #' @param script `list()`\cr
    #' List of cpp script lines.
    #'
    #' @return `character(1)`
    .cpp_full_script = function(script) {
      return(
        paste(c(
          "void Model::assign_group_rate_to_individual(Household &house, indiv i)",
          "{",
          ifelse(all(is.na(lapply(self$defs, `[[`, "by"))),
                 "\trateID r;\n\tUNUSED(house);\n\tUNUSED(i);",
                 "\trateID r;"),
          script,
          "}"
        ), collapse = "\n")
      )
    }
  ),

  ############################################################################
  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################
  ############################################################################

  active = list(
    #' @field param_group_rate_defs `list()`\cr
    #' List containing parameter definitions associated to group rates.
    param_group_rate_defs  = function(value) {
      if (missing(value)) {
        new_defs <- list()
        for (rate_id in names(self$logic_tables)) {
          for (this_table in self$logic_tables[[rate_id]]) {
            for (par_id in setdiff(this_table$.param.name, c("REFERENCE"))) {
              new_defs[[par_id]] <- self$param_defs[[rate_id]]
            }
          }
        }
        return(new_defs)
      }
      stop_h("'param_group_rate_defs' is read only.")
    }
  ),

  ###################
  ##### PUBLIC ######
  ###################

  public = list(
    defs = list(),
    param_defs = list(),
    logic_tables = list(),

    #' @description
    #' New instance of class
    #' [HousepiGroupRate][housepi::HousepiGroupRate].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description
    #' Find group rates required by the model.
    #'
    #' @param refer_file
    #' Path/to/.refer source file.
    require = function(refer_file) {
      # vector of regular expression to use for assertion
      regular_expr <- c(
        # group rate from : .group_rate(..., "string"
        '\\.group_rate\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      )
      # find required group rates from regular expressions
      return(template_find_requirement(refer_file, regular_expr))
    },

    #' @description
    #' Remove a given group_rate definition.
    #'
    #' @param rate_id `character(1)`\cr
    #' The group rate definition to remove.
    #' @param show_warnings `logical(1)`\cr
    #' Wheter to show warnings.
    remove = function(rate_id,
                      show_warnings = TRUE) {

      if (any(names(self$defs) == rate_id)) {
        # drop previous definition
        self$defs[rate_id] <- NULL
        self$param_defs[rate_id] <- NULL
        self$logic_tables[rate_id] <- NULL
        return(TRUE)

      } else if (isTRUE(show_warnings)) {
        warning_h(
          sprintf("Attempt to remove unknown group_rate \"%s\".", this)
        )
      }

      return(FALSE)
    },

    #' @description
    #' Add a new group rate definition.
    #'
    #' @param name `character(1)`\cr
    #' Name of the group rate.
    #' @param by `list()`\cr
    #' Feature(s) associated to group rate as dataset headers.
    #' Corresponding dataset column have to be categorical.
    #' @param prefix `character(1)`\cr
    #' Prefix to parameter names associated to the group rate.
    #' @param reference (named `list()`\cr
    #' Choose reference group providing corresponding features.
    #' Default is the group with the largest sample size.
    #' @param independent
    #' If `TRUE` features are considered separately.
    #' @param prior `character(1)`\cr
    #' Name of the prior distribution.
    #' @param args (named `list()`\cr
    #' Arguments for the prior distribution.
    #' @param rule `character(1)`\cr
    #' Name of the inference proposal rule.
    #' @param scale `double(1)`\cr
    #' Scale parameter for the proposal rule.
    #' Changed automatically during inference if `n_adpat`>0.
    #' @param min `double(1)`\cr
    #' Value of the minimum.
    #' @param max `double(1)`\cr
    #' Value of the maximum.
    #' @param show_warnings `logical(1)`\cr
    #' If `TRUE` (default) show warnings about overwriting definitions.
    #'
    #' @references
    #' Cauchemez, et al. (2009)
    #' Household Transmission of 2009 Pandemic Influenza A (H1N1) Virus in the United States.
    #' New England Journal of Medicine. 361: 2619–2627. doi:10.1056/NEJMoa0905498
    add = function(name,
                   by,
                   prefix = NA,
                   reference = NULL,
                   independent = TRUE,
                   prior = "lnorm",
                   args = list(meanlog = 0, sdlog = 3),
                   min = NA,
                   max = 10,
                   rule = "lnorm",
                   scale = 1) {

      # UPDATE: controls are missing
      if (is.na(name)) {
        stop_h("Incomplete definition: missing name.")
      }

      # remove previous definition if it exists
      if (self$remove(name, FALSE)) {
        warning_h(sprintf("Overwriting existing \"%s\" definition.", name))
      }

      # load reference
      if (!is.null(reference)) {
        if (length(by) != length(reference)) {
          stop_h(
            "Arguments 'by' and 'reference' should have the same dimension."
          )
        } else if (is.list(reference)) {
          if (length(setdiff(names(reference), by)) != 0) {
            stop_h("Unrecognized names in reference list.")
          }

        } else {
          reference <- as.list(reference)
          names(reference) <- by
        }
      }

      # add new definition
      self$defs[[name]] <- list(
        by = by, prefix = prefix,
        reference = reference, independent = isTRUE(independent)
      )

      # add parameter definitions
      self$param_defs[[name]] <- list(
        value = 1, prior = prior, args = args,
        min = min, max = max, rule = rule, scale = scale
      )

      # reorder 
      y <- sort(names(self$defs))
      self$defs <- self$defs[y]
      self$param_defs <- self$param_defs[y]

      invisible(self)
    },

    #' @description
    #' Load group rate definitions either as direct arguments
    #' or from predefined list.
    #'
    #' @param ...
    #' Named arguments with respective list of definition details.
    #' @param from (named `list()`\cr
    #' Named list of parameter definitions.
    #' @param append `logical(1)`\cr
    #' If `TRUE` appends to current definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_list = function(...,
                         from = NULL,
                         append = FALSE) {
      definitions <- load_definition_list(..., from = from)

      # remove previos definitions if not append
      if (isFALSE(append)) self$defs <- NULL

      for (i in seq_along(definitions)){
        this_def <- as.list(definitions[[i]])
        def_name <- names(definitions[i])
        if (length(names(this_def)) == 0) {
          stop_h(
            sprintf("Wrong definition for group rate \"%s\"", def_name),
            "Provide a definition list."
          )
        }
        this_def$name <- def_name
        do.call(self$add, this_def)
      }

      invisible(self)
    },

    #' @description
    #' Assert that group rate definitions match the requirement list.
    #'
    #' @param requirements `character()`\cr
    #' Vector containing required group rates.
    assert = function(requirements) {

      defined <- names(self$defs)

      # point out group rates that are not used
      unused_set <- template_find_unused(requirements, defined)
      # drop unused group rates
      if (length(unused_set) > 0) {
        for (this in unused_set) self$remove(this)
        stop_h(
          sprintf("Dropped unused group rate%s %s.",
                  ifelse(length(unused_set) > 1, "s:", ""),
                  character_to_string_list(unused_set)),
          find_closest_match(unused_set, requirements, min_dist = 3)
        )
      }

      # find missing group rates and set defaults
      missing_set <- template_find_missing(requirements, defined)
      if (length(missing_set) > 0) {
        for (this in missing_set) self$add(this, by = NA)
      }

      invisible(self)
    },

    #' @description
    #' Implement group rate definitions.
    #' Assert if implementated group rates are correct.
    #'
    #' @param dataset
    #' Dataset to be considered.
    #' @param refer_file
    #' Path/to/model.cpp source code.
    #'
    #' @return `logical(1)`\cr
    #' `TRUE` if definitions exist and are up to date.
    cpp_assert = function(dataset,
                          refer_file) {
      if (length(self$defs) > 0) {
        requirements <- self$require(refer_file)
        self$assert(requirements)

        if (length(dataset) > 0) {
          # produce associated script
          private$.cpp_implement(dataset)
          # read source file
          lines <- readLines(refer_file)
          # check if the definitions are already implemented
          return(stringr::str_detect(
            paste(lines, collapse = "\n"),
            stringr::fixed(private$.rate_script)
          ))
          # NOTE I chould also use str_squish as in:
          # > str_detect(str_squish(big), fixed(str_squish(small)))
          # but for the moment it tells me when I am not tidy
        }
      }

      return(TRUE) # default is TRUE (i.e. no definitions)
    },

    #' @description
    #' Overwrite script to source code 'model.cpp'.
    #'
    #' @param src_file `character(1)`\cr
    #' Source file to modify.
    cpp_write = function(src_file) {

      # check that source file exists
      if (!file.exists(src_file)) {
        stop_h("Impossible to find reference file.")
      }

      # read source file
      lines <- readLines(src_file)

      # find the start of the method definition
      start_line <- grep(
        "^void\\s+Model::assign_group_rate_to_individual\\s*\\(",
        lines
      )

      if (length(start_line) == 0) {
        warning_h("Model does not accommodate for group rates.")

      } else {

        # locate method definition opening brace '{'
        opening_brace_line <- NA
        for (i in start_line:length(lines)) {
          if (grepl("\\{", lines[i])) {
            opening_brace_line <- i
            break
          }
        }
        if (is.na(opening_brace_line)) {
          stop_h("Unmatched curly bracket found in reference code.")
        }

        # find method definition closing brace '}'
        brace_level <- 0
        end_line <- NA
        for (i in opening_brace_line:length(lines)) {
          # update brace level based on opening and closing braces
          brace_level <- brace_level +
            stringr::str_count(lines[i], "\\{") -
            stringr::str_count(lines[i], "\\}")
          # closing brace corresponds to level 0
          if (brace_level == 0) {
            end_line <- i
            break
          }
        }
        if (is.na(end_line)) {
          stop_h("Unmatched curly bracket found in reference code.")
        }

        # define the new implementation
        lines <- c(
          lines[1:(start_line - 1)],
          private$.rate_script,
          lines[(end_line + 1):length(lines)]
        )

        # replace the source code
        writeLines(lines, src_file)
      }

      invisible(self)
    }
  )
)
