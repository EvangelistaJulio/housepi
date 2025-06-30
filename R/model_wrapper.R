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
#' Cpp Data Refactor
#' @description
#' Restore factored quanties on dataset return from Cpp implementation.
#'
#' @param df `data.frame`\cr
#' Dataset as output from rcpp model object.
#' @param levels `list()`\cr
#' Levels of the factors in the dataset.
#'
#' @return `data.frame`
#'
#' @keywords internal
cpp_data_refactor = function(df, levels) {
  for (col in names(levels)) df[[col]] <- factor(
    df[[col]],
    levels = seq_along(levels[[col]]),
    labels = levels[[col]]
  )
  return(df)
}

############################################################################
############################################################################
############################## Housepi Wrapper #############################
############################################################################
############################################################################

#' @title
#' Housepi Wrapper R6 Class
#' @description
#' Help with calling methods from the Rcpp model object.
#'
#' @keywords internal
HousepiWrapper <- R6::R6Class(
  classname = "HousepiWrapper",
  inherit = HousepiGenerator,
  cloneable = TRUE,

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field path_to_dll `character(1)`\cr
    #' Path to dynamic library file of the model.
    path_to_dll = function(value) {
      if (missing(value)) return(private$.dll)
      stop_h("'path_to_dll' is read only.")
    },

    #' @field dataset `data.frame`\cr
    #' Dataset with individual level information of current dataset.
    dataset = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) stop_h(
          "Impossible to export dataset before model build."
        )
        return(cpp_data_refactor(private$.rcpp$export_dataset(),
                                 self$dataset_levels))
      }
      stop_h("'dataset' is read only.")
    },

    #' @field datalist `data.frame`\cr
    #' Datalist with multiple rows indexed by individual.
    datalist = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) stop_h(
          "Impossible to export datalist before model build."
        )
        return(cpp_data_refactor(private$.rcpp$export_datalist(),
                                 self$datalist_levels))
      }
      stop_h("'datalist' is read only.")
    },

    #' @field dataset_required `character()`\cr
    #' Dataset headers required by the model.
    dataset_required = function(value) {
      if (missing(value)) {
        # vector of regular expression to use for assertion
        regular_expr <- c(
          # from : .data("string")
          '.data\\s*\\(\\s*["\'](.*?)["\']'
        )

        # retrieve required dataset headers from regular expressions
        refer_file <- self$reference
        return(template_find_requirement(refer_file, regular_expr))
      }
      stop_h("'dataset_required' is read only.")
    },

    #' @field summary `data.frame`\cr
    #' A summary of statuses in the dataset.
    summary = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) stop_h(
          "Impossible to export dataset before model build."
        )
        return(private$.rcpp$export_dataset_report())
      }
      stop_h("'summary' is read only.")
    },

    #' @field summary_house_size `data.frame`\cr
    #' A summary of house sizes in the dataset.
    summary_house_size = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) stop_h(
          "Impossible to export dataset before model build."
        )
        return(private$.rcpp$export_house_size_report())
      }
      stop_h("'summary_house_size' is read only.")
    },

    #' @field setup `list()`\cr
    #' Last setup used by the build process.
    setup = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) stop_h(
          "Impossible to export setup before model build."
        )
        return(private$.rcpp$export_setup())
      }
      stop_h("'setup' is read only.")
    },

    #' @field traces `matrix()` | `list()` of `matrix()`)
    #' Export the traces of last inference call.
    traces = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export traces before inference.")

        } else {
          # load traces from last inference
          df <- private$.rcpp$export_traces()

          if (dim(df)[1] == 0) stop_h("Traces are empty.")

          else if ("CHAIN__" %in% colnames(df)) return(lapply(
            split(seq_len(nrow(df)), df[, "CHAIN__"]),
            function(idx) df[idx, colnames(df) != "CHAIN__", drop = FALSE]
          ))

          else return(df)
        }
      }
      stop_h("'traces' is read only.")
    },

    #' @field traces_no_loglike `matrix()` | `list()` of `matrix()`)
    #' Export the traces of last inference call without loglikelihood.
    traces_no_loglike = function(value) {
      if (missing(value)) {
        trs <- self$traces
        if (is.list(trs)) return(lapply(trs, function(x) {
          x <- x[, colnames(x) != "loglikelihood", drop = FALSE]
        }))
        else return(trs[, colnames(trs) != "loglikelihood", drop = FALSE])
      }
      stop_h("'traces_no_loglike' is read only.")
    },

    #' @field house_loglike_mtx `matrix()`\cr
    #' Household loglikelihood traces.
    house_loglike_mtx = function(value) {
      if (missing(value)) {
        # check model has been built
        if (is.null(private$.rcpp)) stop_h(
          "Impossible to export household traces before inference."
        )
        # extract household loglikelihood traces
        df <- private$.rcpp$export_house_loglike_mtx()
        # check
        if (dim(df)[1] == 0) stop_h("Household traces are empty.")
        return(df)
      }
      stop_h("'household_loglike_traces' is read only.")
    }

  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field last_dataset_cols `character()`\cr
    #' Headers of the last dataset imported.
    last_dataset_cols = NULL,

    #' @field last_datalist_cols `character()`\cr
    #' Headers of the last datalist imported.
    last_datalist_cols = NULL,

    #' @field dataset_levels `list()`\cr
    #' List of levels of the factors in the dataset.
    dataset_levels = NULL,

    #' @field datalist_levels `list()`\cr
    #' List of levels of the factors in the datalist.
    datalist_levels = NULL,

    #' @description
    #' Assert if dataset headers match implementation file.
    #'
    #' @param dataset_headers `character()`\cr
    #' Vector with the dataset headers.
    #' @param requirements `character()`\cr
    #' Vector containing required parameters.
    assert = function(dataset_headers, requirements) {
      # find missing dataset headers
      missing_set <- template_find_missing(requirements,
                                           dataset_headers)
      if (length(missing_set) > 0) stop_h(
        sprintf("Missing dataset header%s for:",
                ifelse(length(missing_set) > 1, "s", "")),
        character_to_string_list(missing_set)
      )
      invisible(self)
    },

    #' @description
    #' Wrapper of Rcpp method to import dataset.
    #'
    #' @param dataset `data.frame`\cr
    #' Dataset with rowise individual information.
    import_dataset = function(dataset) {
      if (length(dataset) > 0) {
        # convert factors to numeric and store levels
        self$dataset_levels <- list()
        for (col in setdiff(colnames(dataset), c("id", "house_id", "status"))) {
          x <- dataset[[col]]
          if (is.character(x) || is.factor(x) || is.logical(x)) {
            x <- factor(x)
            self$dataset_levels[[col]] <- levels(x)
            dataset[[col]] <- as.integer(x)
          }
        }
        # import data.frame in rcpp # note: it calls
        private$.rcpp$import_dataset(fmt_data_input_(dataset))
        # store column names
        self$last_dataset_cols <- colnames(dataset)
      }
      invisible(self)
    },

    #' @description
    #' Wrapper of Rcpp method to import datalist.
    #'
    #' @param datalist `data.frame`\cr
    #' Dataset with rowise individual information.
    import_datalist = function(datalist) {
      if (length(datalist) > 0) {
        # convert factors to numeric and store levels
        self$datalist_levels <- list()
        for (col in setdiff(colnames(datalist), c("id", "house_id", "status"))) {
          x <- datalist[[col]]
          if (is.character(x) || is.factor(x) || is.logical(x)) {
            x <- factor(x)
            self$datalist_levels[[col]] <- levels(x)
            datalist[[col]] <- as.integer(x)
          }
        }
        # import data.frame in rcpp # note: it calls
        private$.rcpp$import_datalist(fmt_data_input_(datalist))
        # store column names
        self$last_datalist_cols <- colnames(datalist)
      }
      invisible(self)
    },

    #' @description
    #' Wrapper of the Rcpp import_setup method.
    #'
    #' @param defs `list()`\cr
    #' List of setup definitions.
    import_setup = function(defs) {
      private$.rcpp$import_setup(fmt_setup_input_(defs))
      invisible(self)
    },

    #' @description
    #' Wrapper of Rcpp method to import parameter definitions.
    #'
    #' @param defs `data.frame`\cr
    #' Formatted parameter definitions.
    import_param = function(defs) {
      private$.rcpp$import_param(data.frame(lapply(defs, as.character)))
      invisible(self)
    },

    #' @description
    #' Wrapper of Rcpp method to import augmented variable definitions.
    #'
    #' @param defs `data.frame`\cr
    #' Formatted augmented variable definitions.
    import_augm = function(defs) {
      private$.rcpp$import_augm(data.frame(lapply(defs, as.character)))
      invisible(self)
    },

    #' @description
    #' Wrapper of Rcpp method to import infectivity profile definitions.
    #'
    #' @param defs `list()`\cr
    #' Infectivity profile definitions.
    import_infec_profile = function(defs) {
      private$.rcpp$import_infec_profile(lapply(defs, as.character))
      invisible(self)
    },

    #' @description
    #' Wrapper of Rcpp method to import multiplicity definitions.
    #'
    #' @param defs `data.frame`\cr
    #' Multiplicity definitions.
    import_multiplicity = function(defs) {
      private$.rcpp$import_multiplicity(lapply(defs, as.character))
      invisible(self)
    },

    #' @description
    #' Wrapper of Rcpp method to import group rates definitions.
    #'
    #' @param defs `data.frame`\cr
    #' Group rates definitions.
    import_group_rates = function(defs) {
      lvls <- self$dataset_levels
      for (k in intersect(colnames(defs), names(lvls))) {
        defs[[k]] <- as.integer(factor(defs[[k]], levels = lvls[[k]]))
      }
      private$.rcpp$import_group_rates(lapply(defs, as.character))
      invisible(self)
    },

    #' @description
    #' Wrapper of the Rcpp build_all_rates method.
    #' FIXME: temporary for contact rates, it will be dropped.
    build_all_rates = function() {
      private$.rcpp$build_all_rates()
      invisible(self)
    },

    #' @description
    #' Export dataset formatted for C++ standalone.
    #'
    #' @return `data.frame`
    export_dataset_to_standalone = function() {
      return(private$.rcpp$export_dataset())
    },

    #' @description
    #' Export datalist formatted for C++ standalone.
    #'
    #' @return `data.frame`
    export_datalist_to_standalone = function() {
      return(private$.rcpp$export_datalist())
    },

    #' @description
    #' Export a list containing full house and model details.
    #' 
    #' @param house_id `character(1)`\cr
    #' Household id to export.
    #'
    #' @return `list()`
    export_state_for_debug = function(house_id) {
      return(private$.rcpp$export_for_debug(house_id))
    },

    #' @description
    #' Showcase instantaneous and cumulative values of infectivity profile.
    #'
    #' @param bins `integer(1)`\cr
    #' Number of bins for continuous profile (ignored if discrete).
    #'
    #' @return `data.frame`
    showcase_infec_profile = function(bins = 2000) {
      # load options
      bins <- as.integer(bins)

      # reformat dataframe
      output <- private$.rcpp$showcase_infec_profile(bins)
      mask <- c("gentime", "inst", "cumul")
      output$profile <- as.data.frame(output[mask])
      output[mask] <- NULL

      return(output)
    },

    #' @description
    #' Wrapper of the Rcpp evaluate method.
    evaluate = function() {
      return(private$.rcpp$evaluate())
    },

    #' @description
    #' Wrapper of the Rcpp infer method.
    infer = function() {
      private$.rcpp$infer()
      invisible(self)
    },

    #' @description
    #' Wrapper of the Rcpp infer method.
    #'
    #' @return `data.frame`
    simulate = function() {
      return(cpp_data_refactor(private$.rcpp$simulate(),
                               self$dataset_levels))
    },

    #' @description
    #' Wrapper of the Rcpp method to estimate_evidence.
    #'
    #' @param scheme `character(1)`\cr
    #' If "adaptive" uses adaptive scheme for thermodynamic ladder (default).
    #' Other available scheme is "power-fractions".
    #' @param n_rungs `integer(1)`\cr
    #' Number of rungs in the thermodynamic ladder.
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
