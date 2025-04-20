#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# model_wrapper.R
#

#' @title
#' Housepi Wrapper R6 Class
#' @description
#' Help with calling methods from the Rcpp model object.
#' @noRd
HousepiWrapper <- R6::R6Class(
  classname = "HousepiWrapper",
  inherit = HousepiGenerator,
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################
  ############################################################################

  active = list(
    #' @field dataset `data.frame`\cr
    #' Dataset with individual level information of current dataset.
    #' If the case, it also return the associated augmented variables.
    dataset = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export dataset before model build.")
        }
        return(private$.rcpp$export_dataset())
      }
      stop_h("'dataset' is read only.")
    },

    #' @field data_required `character()`\cr
    #' Dataset headers required by the model.
    data_required = function(value) {
      if (missing(value)) {
        # vector of regular expression to use for assertion
        regular_expr <- c(
          # from : .data(, "string" or augm(, 'string'
          '.data\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
        )

        # retrieve required data headers from regular expressions
        refer_file <- self$reference
        return(template_find_requirement(refer_file, regular_expr))
      }
      stop_h("'data_required' is read only.")
    },

    #' @field summary `data.frame`\cr
    #' A summary of statuses in the dataset.
    summary = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export dataset before model build.")
        }
        return(private$.rcpp$export_data_report())
      }
      stop_h("'summary' is read only.")
    },

    #' @field summary_house_size `data.frame`\cr
    #' A summary of house sizes in the dataset.
    summary_house_size = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export dataset before model build.")
        }
        return(private$.rcpp$export_house_size_report())
      }
      stop_h("'summary_house_size' is read only.")
    },

    #' @field setup `list()`\cr
    #' Last setup used by the build process.
    setup = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export setup before model build.")
        }
        return(private$.rcpp$export_setup())
      }
      stop_h("'setup' is read only.")
    },

    #' @field `matrix()` | `list()` of `matrix()`)
    #' Export the traces of last inference call.
    traces = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export traces before inference.")

        } else {
          # load traces from last inference
          df <- private$.rcpp$export_traces()

          if (dim(df)[1] == 0) {
            # note: this should never be possible...
            stop_h("Traces are empty.")

          } else if ("CHAIN__" %in% colnames(df)) {
            df_list <- lapply(
              split(seq_len(nrow(df)), df[, "CHAIN__"]),
              function(idx) df[idx, colnames(df) != "CHAIN__", drop = FALSE]
            )
            return(df_list)

          } else {
            return(df)
          }
        }
      }
      stop_h("'traces' is read only.")
    },

    #' @field `matrix()` | `list()` of `matrix()`)
    #' Export the traces of last inference call without loglikelihood.
    traces_no_loglike = function(value) {
      if (missing(value)) {
        trs <- self$traces
        if (is.list(trs)) {
          return(lapply(trs, function(x) {
            x <- x[, colnames(x) != "loglikelihood", drop = FALSE]
          }))
        } else {
          return(trs[, colnames(trs) != "loglikelihood", drop = FALSE])
        }
      }
      stop_h("'traces_no_loglike' is read only.")
    },

    #' @field reference `matrix()`\cr
    #' Household loglikelihood traces.
    house_loglike_mtx = function(value) {
      if (missing(value)) {
        # check model has been built
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export household traces before inference.")
        }
        # extract household loglikelihood traces
        df <- private$.rcpp$export_house_loglike_mtx()
        # check
        if (dim(df)[1] == 0) {
          stop_h("Household traces are empty.")
        }
        return(df)
      }
      stop_h("'household_loglike_traces' is read only.")
    }
  ),

  ############################################################################
  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################
  ############################################################################

  public = list(
    #' @field last_data_cols `character()`\cr
    #' Headers of the last dataset imported.
    last_data_cols = NULL,

    #' @description
    #' Assert if dataset headers match implementation file.
    #'
    #' @param data_headers `character()`\cr
    #' Vector with the dataset headers.
    #' @param requirements `character()`\cr
    #' Vector containing required parameters.
    assert = function(data_headers, requirements) {
      # find missing data headers
      missing_set <- template_find_missing(requirements,
                                           data_headers)
      if (length(missing_set) > 0) {
        stop_h(
          sprintf("Missing dataset header%s for:",
                  ifelse(length(missing_set) > 1, "s", "")),
          character_to_string_list(missing_set)
        )
      }

      invisible(self)
    },

    # @description
    # Wrapper of Rcpp method to import dataset
    #
    # @param dataset `data.frame`\cr
    # Dataset with rowise individual information.
    import_dataset = function(dataset) {
      private$.rcpp$import_dataset(data.frame(lapply(dataset, as.character)))
      if (length(dataset) > 0) self$last_data_cols <- colnames(dataset)
      invisible(self)
    },

    # @description
    # Wrapper of the Rcpp import_setup method.
    #
    # @param defs `list()`\cr
    # List of setup definitions.
    import_setup = function(defs) {
      # change scientific notation locally
      withr::local_options(list(scipen = 500))
      private$.rcpp$import_setup(lapply(defs, as.character))
      invisible(self)
    },

    # @description
    # Wrapper of Rcpp method to import parameter definitions.
    #
    # @param defs `data.frame`\cr
    # Formatted parameter definitions.
    import_param = function(defs) {
      private$.rcpp$import_param(data.frame(lapply(defs, as.character)))
      invisible(self)
    },

    # @description
    # Wrapper of Rcpp method to import augmented variable definitions.
    #
    # @param defs `data.frame`\cr
    # Formatted augmented variable definitions.
    import_augm = function(defs) {
      private$.rcpp$import_augm(data.frame(lapply(defs, as.character)))
      invisible(self)
    },

    # @description
    # Wrapper of Rcpp method to import infectivity profile definitions.
    #
    # @param defs `list()`\cr
    # Infectivity profile definitions.
    import_infec_profile = function(defs) {
      private$.rcpp$import_infec_profile(lapply(defs, as.character))
      invisible(self)
    },

    # @description
    # Wrapper of Rcpp method to import multiplicity definitions.
    #
    # @param defs `data.frame`\cr
    # Multiplicity definitions.
    import_multiplicity = function(defs) {
      private$.rcpp$import_multiplicity(lapply(defs, as.character))
      invisible(self)
    },

    # @description
    # Wrapper of the Rcpp build_all_rates method.
    build_all_rates = function() {
      private$.rcpp$build_all_rates()
      invisible(self)
    },

    # @description
    # Showcase instantaneous and cumulative values of infectivity profile.
    #
    # @param opts `double()`\cr
    # Additional values for the computation.
    # @param bins `integer(1)`\cr
    # Number of bins for continuous profile (ignored if discrete).
    #
    # @return `data.frame`
    showcase_infec_profile = function(opts = list(), bins = 2000) {
      # load options
      opts <- as.list(as.numeric(opts))
      bins <- as.integer(bins)

      # reformat dataframe
      output <- private$.rcpp$showcase_infec_profile(opts, bins)
      mask <- c("gentime", "inst", "cumul")
      output$profile <- as.data.frame(output[mask])
      output[mask] <- NULL

      return(output)
    },

    # @description
    # Wrapper of the Rcpp evaluate method.
    evaluate = function() {
      return(private$.rcpp$evaluate())
    },

    # @description
    # Wrapper of the Rcpp infer method.
    infer = function() {
      private$.rcpp$infer()
      invisible(self)
    },

    # @description
    # Wrapper of the Rcpp infer method.
    #
    # @return `data.frame`
    simulate = function() {
      return(private$.rcpp$simulate())
    },

    # @description
    # Wrapper of the Rcpp method to estimate_evidence.
    #
    # @param scheme `character(1)`\cr
    # If "adaptive" uses adaptive scheme for thermodynamic ladder (default).
    # Other available scheme is "power-fractions".
    # @param n_rungs `integer(1)`\cr
    # Number of rungs in the thermodynamic ladder.
    evidence = function(scheme = "adaptive", n_rungs = 40) {
      # load options
      scheme <- as.character(scheme)
      n_rungs <- as.integer(n_rungs)

      # reformat dataframe
      output <- private$.rcpp$estimate_evidence(scheme, n_rungs)
      mask <- c("inv_temper", "mean_loglike", "var_loglike")
      output$detail <- as.data.frame(output[mask])
      output[mask] <- NULL

      return(output)
    }
  )
)
