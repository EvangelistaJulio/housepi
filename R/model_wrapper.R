#
# model_generator.R
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#

# @title Housepi Wrapper R6 Class
# @description
# Help with calling methods from the Rcpp model object.
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
    # @field dataset `data.frame`\cr
    # Dataset with individual level information of current dataset.
    # If the case, it also return the associated augmented variables.
    dataset = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export dataset before model build.")
        }
        return(private$.rcpp$export_dataset())
      }
      stop_h("'dataset' is read only.")
    },

    # @field summary `data.frame`\cr
    # A summary of statuses in the dataset.
    summary = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export dataset before model build.")
        }
        return(private$.rcpp$export_data_report())
      }
      stop_h("'summary' is read only.")
    },

    # @field summary_house_size `data.frame`\cr
    # A summary of house sizes in the dataset.
    summary_house_size = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export dataset before model build.")
        }
        return(private$.rcpp$export_house_size_report())
      }
      stop_h("'summary_house_size' is read only.")
    },

    # @field setup `list()`\cr
    # Last setup used by the build process.
    setup = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export setup before model build.")
        }
        return(private$.rcpp$export_setup())
      }
      stop_h("'setup' is read only.")
    },

    # @field `matrix()` | `list()` of `matrix()`)
    # Export the traces of last inference call.
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

    # @field `matrix()` | `list()` of `matrix()`)
    # Export the traces of last inference call without loglikelihood.
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

    # @field traces_as_coda `coda::mcmc` | `coda::mcmc.list`\cr
    # Traces of last inference as a `coda` object.
    traces_as_coda = function(value) {
      if (missing(value)) return(coerce_traces_to_coda(self$traces))
      stop_h("'traces_as_coda' is read only.")
    },

    # @field traces_results `data.frame`\cr
    # Results extracted from the traces of last inference.
    traces_results = function(value) {
      if (missing(value)) {
        return(retrieve_traces_results(self$traces_no_loglike))
      }
      stop_h("'traces_results' is read only.")
    },

    # @field reference `matrix()`\cr
    # Household loglikelihood traces.
    household_loglike_matrix = function(value) {
      if (missing(value)) {
        # check model has been built
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export household traces before inference.")
        }
        # extract household loglikelihood traces
        df <- private$.rcpp$export_household_loglike_matrix()
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
    # @description
    # Show formatted summary of loaded dataset.
    #
    # @param pct_format `character(1)`\cr
    # Formatting string for the percentage passed to `sprintf`.
    # It expects 1 double value (fraction).
    #
    # @return `flextable`
    show_summary = function(pct_format = "%.1f%%") {
      if (require_suggested_library("flextable")) {

        combined <- self$summary
        combined$pct <- sapply(combined$pct,
                               function(x) sprintf(pct_format, x))

        ft <- flextable::flextable(combined)
        ft <- flextable::separate_header(ft)
        ft <- flextable::theme_vanilla(ft)
        ft <- flextable::align(ft, align = "left", part = "all")
        ft <- flextable::merge_v(ft, j = ~ type)
        ft <- flextable::bold(ft, j = ~ type, bold = TRUE)
        ft <- flextable::autofit(ft, add_w = 0, add_h = 0)
        ft <- flextable::labelizor(
          ft, part = "all",
          labels = c(type = "", comment = "", pct = "[%]")
        )
        ft <- flextable::colformat_num(ft, big.mark = "",
                                       na_str = "", nan_str = "")
        return(ft)
      }
    },

    # @description
    # Show formatted summary of housheold sizes.
    #
    # @param pct_format `character(1)`\cr
    # Formatting string for the percentage passed to `sprintf`.
    # It expects 1 double value (fraction).
    #
    # @return `flextable`
    show_house_size_summary = function(pct_format = "%.1f%%") {
      if (require_suggested_library("flextable")) {
        combined <- self$summary_house_size
        combined$pct <- sapply(combined$pct,
                               function(x) sprintf(pct_format, x))
        ft <- flextable::flextable(combined)
        ft <- flextable::separate_header(ft)
        ft <- flextable::theme_vanilla(ft)
        ft <- flextable::align(ft, align = "left", part = "all")
        ft <- flextable::bold(ft, j = ~ comment + value, bold = TRUE)
        ft <- flextable::autofit(ft, add_w = 0, add_h = 0)
        ft <- flextable::labelizor(
          ft, part = "all",
          labels = c(comment = "", value = "n.", pct = "[%]")
        )
        ft <- flextable::colformat_num(ft, big.mark = "")
        return(ft)
      }
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Wrapper of the Rcpp import_dataset method.
    #
    # @param dataset `data.frame`\cr
    # Dataset with rowise individual information.
    import_dataset = function(dataset) {
      private$.rcpp$import_dataset(
        data.frame(lapply(dataset, as.character))
      )
      invisible(self)
    },

    # @description
    # Wrapper of the Rcpp import_setup method.
    #
    # @param setup_list `list()`\cr
    # Setup definitions.
    import_setup = function(setup_list) {
      # change scientific notation locally
      withr::local_options(list(scipen = 500))
      private$.rcpp$import_setup(
        lapply(setup_list, as.character)
      )
      invisible(self)
    },

    # @description
    # Wrapper of the Rcpp import_parameter method.
    #
    # @param defs `data.frame`\cr
    # Formatted parameter definitions.
    import_parameter = function(defs) {
      private$.rcpp$import_parameter(
        data.frame(lapply(defs, as.character))
      )
      invisible(self)
    },

    # @description
    # Wrapper of the Rcpp import_augmented method.
    #
    # @param defs `data.frame`\cr
    # Formatted augmented variable definitions.
    import_augmented = function(defs) {
      private$.rcpp$import_augmented(
        data.frame(lapply(defs, as.character))
      )
      invisible(self)
    },

    # @description
    # Wrapper of the Rcpp import_infec_profile method.
    #
    # @param defs `list()`\cr
    # Infectivity profile definitions.
    import_infectivity_profile = function(defs) {
      private$.rcpp$import_infectivity_profile(
        lapply(defs, as.character)
      )
      invisible(self)
    },

    # @description
    # Wrapper of the Rcpp build_all_rates method.
    build_all_rates = function() {
      private$.rcpp$build_all_rates()
      invisible(self)
    },

    # @description
    # Export instantaneous and cumulative values of the infectivity profile.
    #
    # @param opts `double()`\cr
    # Additional values for the computation.
    # @param n_bins `integer(1)`\cr
    # Number of bins for continuous profile (ignored if discrete).
    #
    # @return `data.frame`
    export_infectivity_profile = function(opts = list(),
                                          n_bins = 2500) {
      # load options
      opts <- as.list(as.numeric(opts))
      n_bins <- as.integer(n_bins)

      # reformat dataframe
      output <- private$.rcpp$export_infectivity_profile(opts, n_bins)
      output$profile <- as.data.frame(output[c("gentime", "inst", "cumul")])
      output[c("gentime", "inst", "cumul")] <- NULL

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
    # Wrapper of the Rcpp estimate_evidence method.
    #
    # @param scheme `character(1)`\cr
    # If "adaptive" uses adaptive scheme for thermodynamic ladder (default).
    # Other available scheme is "power-fractions".
    # @param n_rungs `integer(1)`\cr
    # Number of rungs in the thermodynamic ladder.
    estimate_evidence = function(scheme = "adaptive",
                                 n_rungs = 40) {
      return(private$.rcpp$estimate_evidence(scheme, n_rungs))
    }
  )
)
