#
# group_rate.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

# formatting options
TAB <- "    "
TAB2 <- strrep(TAB, 2)

# @description
# return logic table of group rates
# @param dataset
# the dataset from which the etract group rates
# @param features
# features as list of dataset column names
group_rate_logic_table <- function(dataset,
                                   features) {
  # convert features
  features <- as.character(features)

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
  logic_table$..counts <- multifeat_count[multifeat_count > 0]

  return(logic_table)
}

# @ description
# Associate parameter name to group rate logic table
# @param logic_table
# @param prefix
# @param reference
associate_group_names <- function(logic_table,
                                  prefix = NA,
                                  reference = NULL) {
  # extract reference index
  if (is.null(reference)) {
    idx_refer <- which.max(logic_table$..counts)
    # drop counts
    logic_table$..counts <- NULL

  } else {
    # drop counts
    logic_table$..counts <- NULL

    reference <- as.list(reference)[names(logic_table)]
    idx_refer <- unname(which(
      apply(logic_table, 1, function(row) all(row == unlist(reference)))
    ))
    if (length(idx_refer) != 1) {
      stop("Provided reference group is not valid.")
    }
  }

  # format prefix
  if (is.na(prefix) || !is.character(prefix)) {
    prefix <- ""
  } else {
    prefix <- paste0(prefix, "_")
  }

  # associate parameter names
  logic_table$..associate <- apply(
    logic_table, 1, function(row) {
      paste0(prefix, paste(paste0(names(logic_table), row), collapse = "_"))
    }
  )
  logic_table$..associate[idx_refer] <- "REFERENCE"

  return(logic_table)
}

# @ description
# construct the C++ if clauses
# @param r
# row of the group rate row
# @param start
# if `TRUE` the if clause is the first
cpp_ifclause_line <- function(r,
                              is_ref = FALSE) {
  mask <- setdiff(names(r), "..associate")
  # if condition
  condition <- paste(
    paste0(
      sprintf('house.data(i, "%s") == ', names(r[mask])),
      r[mask]), collapse = ") && ("
  )
  # wrap all together
  c(sprintf("%s%sif ((%s))", TAB, ifelse(is_ref, "", "else "), condition),
    sprintf("%semplace_group_rate(r, %s);",
            TAB2,
            ifelse(is_ref, "REFERENCE", sprintf('"%s"', r[["..associate"]])))
  )
}


# @description
# final wrapper
specify_group_rate <- function(dataset,
                               features,
                               rate_id,
                               prefix = NA,
                               reference = NULL,
                               independent = TRUE) {
  # default prefix
  if (is.na(prefix)) {
    prefix <- format_string(tolower(rate_id))
    prefix <- substring(prefix, 1, 3)
  }

  features <- as.list(features)
  if (isFALSE(independent)) {
    features <- list(features)
  }

  # initialize script
  script <- c()

  for (this_feat in features) {
    logic_table <- group_rate_logic_table(dataset, this_feat)
    logic_table <- associate_group_names(logic_table, prefix, reference)
    # mask reference case
    mask_refer <- logic_table$..associate == "REFERENCE"
    script <- c(
      script,
      # reference case
      cpp_ifclause_line(logic_table[mask_refer, ], TRUE),
      # other cases
      apply(logic_table[!mask_refer, ], 1, cpp_ifclause_line),
      # close with missing case
      sprintf("%selse", TAB),
      sprintf("%semplace_group_rate(r, MISSING_RATE);", TAB2)
    )
  }

  return(c("{", sprintf('%srateID r = "%s";', TAB, rate_id), script, "}"))
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# @title Housepi Group Rate R6 Class
# @description
# Help with group rate definitions of model Rcpp object.
# @export HousepiGroupRateClass
# @docType class
HousepiGroupRateClass <- R6::R6Class(
  classname = "HousepiGroupRateClass",
  cloneable = FALSE,

  ###################
  ##### PUBLIC ######
  ###################

  public = list(
    defs = NULL,

    # @description
    # Assign the empty parameter data.frame
    initialize = function() {
      self$drop()
      invisible(self)
    },

    # @description
    # Remove all current group rate definitions.
    drop = function() {
      self$defs <- list()
      invisible(self)
    },

    # @description
    # Add a new group rate definition.
    # @param name (`character(1)`)\cr
    # Name of the group rate.
    # @param features (`list()`)\cr
    # Feature(s) associated to group rate as dataset headers.
    # Corresponding dataset column have to be categorical.
    # @param prefix (`character(1)`)\cr
    # Prefix to parameter names associated to the group rate.
    # @param (named `list()`)\cr
    # Choose reference group providing corresponding features.
    # Default is the group with the largest sample size.
    # @param independent
    # If `TRUE` features are considered separately.
    add = function(name,
                   features,
                   prefix = NA,
                   reference = NULL,
                   independent = TRUE) {

      # UPDATE: controls are missing
      self$defs[[name]] <- list(
        features = features,
        prefix = prefix,
        reference = reference,
        independent = independent
      )
    },

    # @description
    # Assert group rate definitions from the implementation script.
    # @param refer_file
    # @param fill (`logical(1)`)\cr
    # If `TRUE` fills the missing group rate with template definitions.
    assert = function(refer_file, fill = FALSE) {

      # read the script and extract parameters
      script <- readChar(refer_file, file.info(refer_file)$size)

      # group rate from : house.group_rate(..., "string"
      regular_expr <- 'house.group_rate\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      model_rate_set <- unlist(matches[[1]][, 2])

      # re-order alphabetically
      model_rate_set <- unique(model_rate_set[order(model_rate_set)])

      # add missing parameters with standard options
      missing <- c()
      for (name in model_rate_set) {
        if (!(name %in% self$defs[["name"]])) {
          missing <- c(missing, name)
        }
      }

      # eventually fill with template definitions
      if (length(missing) > 0) {
        if (isTRUE(fill)) {
          stop(cli::format_error(c(
              "Template group rates are not yet implemented."
            )))
        } else {
          oss <- sprintf("Missing group rate(s): '%s'", missing[[1]])
        }

        for (name in missing[-1]) {
          oss <- sprintf("%s, '%s'", oss, name)
        }
        cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
      }

      # point out group rates that are not used
      unused <- c()
      for (name in names(self$defs)) {
        if (!(name %in% model_rate_set)) {
          unused <- c(unused, name)
        }
      }
      if (length(unused) > 0) {
        oss <- paste0("Group rate(s) not used in model: '", unused[[1]], "'")
        for (name in unused[-1]) {
          oss <- paste0(oss, ", '", name, "'")
        }
        cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
      }

      invisible(self)
    },

    # @description
    # Generate a script to be printed within source code model.cpp
    # @param dataset
    # Dataset to be considered.
    cpp_implement = function(dataset) {

      # initialize lines of the script
      rate_script <- c()

      for (i in seq_along(self$defs)){
        this_def <- self$defs[[i]]
        this_def[["rate_id"]] <- names(self$defs[i])
        this_def[["dataset"]] <- dataset # UPDATE : this is not great...
        rate_script <- c(rate_script,
                         do.call(specify_group_rate, this_def))
      }
      return(
        paste(c(
          "/* Assign all group rates for each individual.",
          "   Generated automatically -> do not edit by hand */",
          "void Model::assign_group_rate_to_individual(Household &house, indiv i)",
          "{",
          sapply(rate_script, function(x) { paste0(TAB, x) }),
          "}"
        ),
        collapse = "\n")
      )
    },

    # @description
    # Assert if script would change with the new implementation.
    # @param dataset
    # Dataset to be considered.
    # @param source_file
    # Path/to/model.cpp source code.
    # @return (`logical(1)`)\cr
    cpp_assert = function(dataset,
                          source_file) {
      # check if there are user defined group rates
      if (length(self$defs) > 0) {
        # produce associated script
        rate_script <- self$cpp_implement(dataset)
        # read source file
        lines <- readLines(source_file)
        # check if the definitions are already implemented
        return(!stringr::str_detect(paste(lines, collapse = "\n"),
                                    stringr::fixed(rate_script)))
      }

      # default is false (i.e. no definitions)
      return(FALSE)
    },

    # @description
    # Overwrite script to source code model.cpp
    # @param dataset
    # Dataset to be considered.
    # @return logical
    # Wheter code requires recompilation
    cpp_write = function(dataset,
                         source_file) {
      # produce associated script
      rate_script <- self$cpp_implement(dataset)
      # read source file
      lines <- readLines(source_file)

      # find the start of the method definition
      start_line <- grep("^void\\s+Model::assign_group_rate_to_individual\\s*\\(", lines)

      if (length(start_line) == 0) {
        cli::cli_alert_warning("Model does not accommodate for group rates.")

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
          stop("Unmatched curly bracket(s) found in source code.")
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
          stop("Unmatched curly bracket(s) found in source code.")
        }

        # define the new implementation
        lines <- c(
          lines[1:(start_line - 1)],
          rate_script,
          lines[(end_line + 1):length(lines)]
        )

        # replace the source code
        writeLines(lines, source_file)
      }

      invisible(self)
    }
  )
)
