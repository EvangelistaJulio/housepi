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
# model_class.R
#

############################################################################
############################################################################
############################## Household Model #############################
############################################################################
############################################################################

#' @title
#' Housepi Model
#' @description
#' Compile household model, manage your data
#' and variable definitions from a single class.
#'
#' @docType class
HousepiModel <- R6::R6Class(
  classname = "HousepiModel",
  inherit = HousepiGenerator,
  cloneable = TRUE,
  lock_class = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    .param = NULL,
    .augm = NULL,
    .group_rates = NULL,
    .infec_profile = NULL,
    .crash_house = NA,

    # @description
    # Augmented variables required by implementation file and other definitions.
    .get_all_required_augm = function() {
      c(private$.augm$require(private$.reference()), # by the model
        private$.infec_profile$augm_required)        # by infec profile
    },

    # @description
    # Dataset headers required by implementation file and other definitions.
    .get_all_required_data = function() {
      required <- c()
      if (!private$.is_blank()) required <- private$.data_required()
      c(required, # by the model
        private$.augm$data_required,          # by augmented
        private$.infec_profile$data_required) # by infec profile
    },

    # @description
    # Group rates required by implementation file and other definitions.
    .get_all_required_group_rates = function() {
      c(private$.group_rates$require(private$.reference()), # by the model
        private$.augm$group_rates_required,                 # by augmented
        private$.infec_profile$group_rates_required)        # by infec profile
    },

    # @description
    # Parameters required by implementation file and other definitions.
    .get_all_requried_param = function() {
      if (private$.is_blank()) required <- private$.param$defs$name
      else required <- private$.param$require(private$.reference())
      c(required,                              # by the model
        private$.group_rates$param_required,   # by group rates
        private$.augm$param_required,          # by augmented
        private$.infec_profile$param_required) # by infec profile
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Assert if augmented definitions match the required.
    # @param capital `logical(1)` to raise error.
    .assert_augm = function(capital = TRUE) {
      if (private$.is_blank()) return(NULL)
      private$.augm$assert(private$.get_all_required_augm(),
                           private$.infec_status_set,
                           capital)
    },

    # @description
    # Assert if data_wide headers match implementation file.
    .assert_data_wide = function(headers) {
      if (length(headers) == 0) stop_h(
        "Impossible to proceed. No data was provided."
      )
      required <- private$.get_all_required_data()
      template_missing_stop(required, headers, "data")
      invisible(self)
    },

    # @description
    # Assert group rate definitions from the model implementation script.
    # @param capital `logical(1)` to raise error.
    .assert_group_rates = function(capital = TRUE) {
      if (private$.is_blank()) return(NULL)
      private$.group_rates$assert(private$.get_all_required_group_rates(),
                                  capital)
    },

    # @description
    # Assert if parameter definitions match required
    # @param capital `logical(1)` to raise error.
    .assert_param = function(capital = TRUE) {
      private$.param$assert(private$.get_all_requried_param(), capital)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Generate multiplicity definitions
    # UPDATE: this will be moved inside C++ standalone in a future version.
    #
    # @return `data.frame`
    .generate_multiplicity_defs = function() {
      # check if model allows for multiplicity
      refer_file <- private$.reference()
      regexp <- ".multiplicity\\s*\\(\\s*\\)"

      if (!template_assert_expression(refer_file, regexp)) return(NULL)

      # note: this might be a modeler error
      if (length(private$.get_all_required_augm()) > 0) warning_h(
        "Model with augmentated variable called multiplicity engine."
      )

      # load data_wide
      # UPDATE: this "back & forth" of data_wide is not optimal!
      data_wide <- self$data_wide

      headers <- c(
        # data headers required by implementation
        private$.get_all_required_data(),
        # data headers required by group rates
        unlist(lapply(
          private$.group_rates$logic_tables,
          function(x) unlist(lapply(x, names))
        ))
      )
      # remove logic tables auxiliary column names
      headers <- setdiff(unique(headers), c("n", "param_name"))

      # list of splitted data by house_id
      split_data <- split(data_wide[c("status", headers)], data_wide$house_id)
      # sort individuals in each household by all columns
      sorted_data <- lapply(split_data, function(x) {
        x[do.call(order, x), , drop = FALSE]
      })
      # create hash strings from sorted rows
      hashes <- sapply(sorted_data, function(x) {
        paste(apply(x, 1, paste, collapse = "-"), collapse = "|")
      })
      # regroup households by hash
      hash_groups <- split(names(hashes), hashes)

      # create summary table with class_id and multiplicity
      mult_hhtable <- data.frame(class_id = sapply(hash_groups, `[`, 1),
                                 multiplicity = lengths(hash_groups))

      # use hashes to produce multiplicity table
      defs <- data.frame(house_id = unique(data_wide$house_id),
                         multiplicity = 0,
                         class_id = NA_integer_)
      rownames(defs) <- defs$house_id
      defs[mult_hhtable$class_id, "multiplicity"] <- mult_hhtable$multiplicity

      for (i in seq_len(nrow(mult_hhtable))) {
        defs[hash_groups[[i]], "class_id"] <- mult_hhtable$class_id[i]
      }
      defs
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Pre-build to deal with group rates.
    .build_group_rates = function() {
      # proceed only if definitions exist and data loaded
      if (length(private$.group_rates$defs) == 0) return(invisible(self))
      if (length(private$.cols_data_wide) == 0) {
        warning_h(
          "Missing data needed for building group rates.",
          "Ensure to provide all necessary data before proceeding."
        )
        return(invisible(self))
      }
      # load data_wide
      # UPDATE: this "back & forth" of data_wide is not optimal!
      # UPDATE: count in the wrapper while storing data_wide_levels
      data_wide <- self$data_wide

      # implement group rates
      private$.group_rates$implement(data_wide, private$.info_data_wide)

      # append parameter definitions if missing
      new_parameters <- private$.group_rates$parameters
      # find missing parameter definitions
      mask <- setdiff(names(new_parameters), private$.param$defs$name)
      new_parameters <- new_parameters[mask]
      if (length(new_parameters) > 0) {
        xquos <- convert_list_to_quosures(new_parameters)
        private$.param$load_quosures(xquos, append = TRUE)
      }

      invisible(self)
    },

    # @description
    # Get individual for blank model augmented draws.
    # UPDATE: not robust if data has only default columns
    .get_blank_id_draws = function(id_draws) {
      # if no default columns, assume debugging mode
      if (private$.is_blank() &&
            is.null(id_draws) &&
            !is.null(private$.augm$defs) &&
            length(setdiff(private$.cols_data_wide,
                           c("house_id", "id", "status"))) == 0) {
        id_draws <- template_blank_data_wide_()$id[1]
      }
      id_draws
    },

    # @description
    # Full model build with complete list of definitions.
    .build_full = function(steps = 0, adapt = 0,
                           burnin = 0, thin = 1,
                           chains = 1, threads = 1,
                           id_draws = NULL,
                           seed = NA_integer_, progress = FALSE) {
      # reset crash house
      private$.crash_house <- NA
      # import setup
      private$.rcpp$import_setup(char_setup_input_(list(
        steps = as.integer(steps),
        adapt = as.integer(adapt),
        burnin = as.integer(burnin),
        thin = as.integer(thin),
        chains = as.integer(chains),
        threads = as.integer(threads),
        id_draws = private$.get_blank_id_draws(id_draws),
        seed = generator_adjust_seed(seed),
        progress = as.integer(progress)
      )))

      # assert group rate definitions
      private$.assert_group_rates()
      # build group rate definitions
      private$.build_group_rates()

      # build contact rates (UPDATE: temporary)
      private$.rcpp$build_all_rates()

      # assert parameter definitions
      private$.assert_param()
      # import parameters
      private$.rcpp$import_param(private$.param$defs)

      # import group rates
      private$.rcpp$import_group_rates(
        private$.group_rates$expose_to_cpp(private$.info_data_wide)
      )

      # assert augmented variable definitions
      private$.assert_augm()

      # import infectivity profile
      private$.rcpp$import_infec_profile(private$.infec_profile$defs)

      # import augmented variables
      private$.rcpp$import_augm(private$.augm$defs)

      # associate multiplicity engine
      private$.rcpp$import_multiplicity(private$.generate_multiplicity_defs())
      invisible(self)
    },

    # UPDATE:
    # this could be designed to be more helpful
    # e. g. create a standalone folder somewhere with correct links to headers
    .export_to_cpp_standalone = function(outpath = "inst/standalone/") {
      private$.assert_group_rates()
      private$.build_group_rates()

      write.csv(private$.rcpp$export_data_wide(),
                file.path(outpath, "data_wide.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$.generate_multiplicity_defs(),
                file.path(outpath, "multiplicity.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$.rcpp$export_data_long(),
                file.path(outpath, "data_long.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$.group_rates$expose_to_cpp(private$.info_data_wide),
                file.path(outpath, "group_rates.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$.augm$defs,
                file.path(outpath, "augm.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$.param$defs,
                file.path(outpath, "param.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(char_setup_input_(self$setup),
                file.path(outpath, "setup.csv"),
                quote = FALSE, row.names = FALSE)
      write.csv(private$.infec_profile$defs,
                file.path(outpath, "infec.csv"),
                quote = FALSE, row.names = FALSE)
      invisible(self)
    },

    # @description
    # Sugar for log-likelihood computation errors in try-catch.
    .helper_log_lik_error = function(e) {
      h <- stringr::str_extract(e$message, "(?<=household ')[^']+(?=')")
      i <- stringr::str_extract(e$message, "(?<=individual ')[^']+(?=')")
      if (!is.na(i)) h <- subset(self$data_wide, id == i)$house_id
      if (!stringr::str_count(e$message, "While computing") ||
            is.na(h)) stop_h(e$message)
      private$.crash_house <- h
      # store house_id somewhere
      stop_h(
        e$message,
        "Use '$export_house_snapshot' to show last status of the hosuehold."
      )
    },

    # Helper with parameter table import.
    .helper_import_table = function(param_table) {
      # extract parameter draws from table
      if (inherits(param_table, "HousepiPosteriorDraws")) {
        return(posterior::as_draws_matrix(param_table$param_draws))
      } else if (posterior::is_draws(param_table)) {
        return(posterior::as_draws_matrix(param_table))
      }
      tryCatch({
        return(as.data.frame(param_table))
      }, error = function(e) {
        stop_h("Type of param_table not supported.")
      })
    },

    # Run measure function on posterior simulations.
    # UPDATE: this has to be improved with loops in C++
    .run_measure = function(param_table, apply_measure, progress) {
      # attach param_table
      param_table <- private$.helper_import_table(param_table)
      private$.rcpp$import_param_for_prediction(param_table)

      # define iterator, eventually with progress bar
      iterator <- seq_len(dim(param_table)[1])
      # if progress not provided, import from model
      progress <- ifelse(is.na(progress), self$setup$progress, progress)
      # FIXME: can't manually set the frequency here
      if (is.numeric(progress)) progress <- TRUE
      if (isTRUE(progress)) iterator <- cli::cli_progress_along(
        iterator, name = "simulations", type = "tasks", clear = FALSE
      )

      # static data (not returned by simulation)
      static_df <- private$.rcpp$data_wide_for_simulation()
      # run simulation and apply all measures
      # NOTE: simulate_for_expec from 0 to table_size - 1 !
      lapply(iterator, function(r) {
        apply_measure(private$.refactor_data_wide(
          cbind(static_df, private$.rcpp$simulate_for_prediction(r - 1)[[1]])
        ))
      })
    },

    # @description
    # Estimate the model logevidence (or marginal likelihood).
    # It uses the method of power posteriors with different available schemes.
    #
    # @param scheme `r rd_evi_scheme_`
    # @param rungs `r rd_evi_rungs_`
    # @param steps `r rd_smpl_steps_`
    # @param burnin `r rd_smpl_burn_`
    # @param thin `r rd_smpl_thin_`
    # @param adapt `r rd_smpl_adapt_`
    # @param chains `r rd_smpl_chain_`
    # @param threads `r rd_smpl_thread_`
    # @param seed `r rd_seed_`
    # @param progress `r rd_progress_`
    #
    # @return `numeric(1)`
    logevidence = function(scheme = "power-fractions", rungs = 40,
                           steps = 10000, burnin = 0,
                           thin = 1, adapt = 0,
                           chains = 1, threads = 1,
                           seed = NA_integer_, progress = TRUE) {
      # UPDATE: developer warning
      warning_h(
        "The method to estimate logevidence is currently under development.",
        "Its usage is deprecated."
      )

      # internal build function
      private$.build_full(steps = steps, burnin = burnin,
                          thin = thin, adapt = adapt,
                          chains = chains, threads = threads,
                          seed = seed, progress = progress)
      # load options
      scheme <- as.character(scheme)
      rungs <- as.integer(rungs)
      # evaluate marginal from the internal Rcpp object
      output <- private$.rcpp$logevidence(scheme, rungs)
      # reformat dataframe
      mask <- c("inv_temper", "mean_log_lik", "var_log_lik")
      output$detail <- as.data.frame(output[mask])
      output[mask] <- NULL
      output
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field data_wide
    #' A `data.frame` with individual level information.
    data_wide = function(value) {
      if (missing(value)) {
        return(private$.refactor_data_wide(private$.rcpp$export_data_wide()))
      }
      self$add_data_wide(value)
    },

    #' @field data_long
    #' A `data.frame` with additional data rows indexed by individual.
    data_long = function(value) {
      if (missing(value)) {
        return(private$.refactor_data_long(private$.rcpp$export_data_long()))
      }
      self$add_data_long(value)
    },

    #' @field parameters
    #' A `data.frame` with parameter formatted definitions.
    parameters = function(value) {
      if (missing(value)) {
        private$.build_group_rates()
        return(housepi_defin_showcaser(private$.param$defs))
      }
      xquos <- convert_list_to_quosures(value)
      private$.param$load_quosures(xquos, append = FALSE)
    },

    #' @field augmented
    #' A`data.frame` with augmented formatted definitions.
    augmented = function(value) {
      if (missing(value)) return(housepi_defin_showcaser(private$.augm$defs))
      xquos <- convert_list_to_quosures(value)
      private$.augm$load_quosures(xquos, append = FALSE)
    },

    #' @field infec_profile
    #' A `list` with infectivity profile definitions.
    infec_profile = function(value) {
      # UPDATE: reformatting current rcpp_export waiting for status specific
      if (missing(value)) {
        # get custom definition
        custom <- private$.infec_profile$showcase()
        if (length(custom) > 0) return(custom)

        # get template/previously implemented
        stored <- private$.rcpp$export_infec_profile()
        if (length(stored) > 0) {
          output <- data.frame(name = stored$name,
                               args = base::I(list(stored$args)),
                               min = stored$min, max = stored$max,
                               discr = stored$discr,
                               stringsAsFactors = FALSE)
          class(output$args) <- "list"
          return(output)
        }
        return(NULL)
      }
      do.call(self$add_infec_profile, value)
    },

    #' @field group_rates
    #' A `list` with group rate definitions.
    group_rates = function(value) {
      if (missing(value)) {
        private$.assert_group_rates()
        private$.build_group_rates()
        return(private$.group_rates$tidy_table)
      }
      self$drop_group_rates()
      xquos <- convert_list_to_quosures(value)
      private$.group_rates$load_quosures(xquos, append = FALSE)
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' Create a new object of class 'HousepiModel'
    #' used to run inference and simulation routines.
    #'
    #' @param choice `r rd_mdl_choice_`
    #' @param as_template `r rd_mdl_astmpl_`
    #' @param force `r rd_mdl_force_`
    #' @param quiet `r rd_mdl_quiet_`
    #' @param ... `r rd_mdl_flags_`
    initialize = function(choice = NA_character_, as_template = FALSE,
                          force = FALSE, quiet = FALSE, ...) {
      # initialize internal fields
      self$drop_all()
      # generate model wrapper
      super$initialize(choice = choice, as_template = as_template,
                       force = force, quiet = quiet, ...)
      # import blank template data_wide
      if (private$.is_blank()) self$add_data_wide(template_blank_data_wide_())

      invisible(self)
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      sms <- c(i = sprintf("Household model {.pkg %s}%s.",
                           spaced_template_choice(self$name),
                           ifelse(private$.is_template(), " (template)", "")))
      # UPDATE: to be rationalize in the future if more options are included
      if ("_HOUSEPI_DISCRETE_TIME" %in% private$.last_user_flags) sms <- c(
        sms, "Model is set to {.emph discrete time} mode."
      )
      if (length(private$.get_all_requried_param()) > 0) sms <- c(sms,
        fmt_assert_cli("Parameters", private$.assert_param(FALSE))
      )
      if (length(private$.get_all_required_augm()) > 0) sms <- c(sms,
        fmt_assert_cli("Augmented", private$.assert_augm(FALSE))
      )
      if (length(private$.get_all_required_group_rates()) > 0) sms <- c(sms,
        fmt_assert_cli("Group rates", private$.assert_group_rates(FALSE),
                       minor = TRUE)
      )
      if (length(private$.rcpp$export_infec_profile()) > 0) {
        ip_sms <- NULL
        if (length(private$.infec_profile$showcase()) == 0) {
          ip_sms <- c("Using templated infectivity profile.",
                      "A custom profile can be specified for model.")
        }
        sms <- c(sms,
                 fmt_assert_cli("Infectivity profile", ip_sms, minor = TRUE))
      }
      cli::cli_inform(sms)
      invisible(self)
    },

    #' @description
    #' Show model description in RStudio Help window.
    help = function() {
      open_model_helper_(private$.path())
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Make an internal copy of the provided household data_wide.
    #'
    #' @param data_wide `r rd_add_data_wide_`
    add_data_wide = function(data_wide) {
      if (length(data_wide) > 0) {
        # check that data_wide type is correct
        if (!is.data.frame(data_wide)) stop_h(
          "Unsupported type for data_wide.",
          "Dataset should be a data.frame."
        )

        # assert data_wide headers
        private$.assert_data_wide(colnames(data_wide))
        # convert factors
        embed <- cpp_data_factorize(data_wide)
        # store levels
        private$.info_data_wide <- embed$info
        # store column names
        private$.cols_data_wide <- embed$cols
        # store all postive status
        private$.infec_status_set <- embed$status

        # import data.frame in rcpp
        private$.rcpp$import_data_wide(embed$data)

      } else {
        # assert headers of last imported data_wide
        private$.assert_data_wide(private$.cols_data_wide)
      }

      invisible(self)
    },

    #' @description Summarise statuses in loaded data.
    #' @return`data.frame`
    summarise_cases = function() {
      private$.rcpp$summarise_cases()
    },

    #' @description Summary house sizes in loaded data.
    #' @return `data.frame`
    summarise_house_sizes = function() {
      private$.rcpp$summarise_house_sizes()
    },

    #' @description
    #' Make an internal copy of the provided data_long.
    #'
    #' @param data_long `r rd_add_data_long_`
    add_data_long = function(data_long) {
      if (!is.data.frame(data_long)) stop_h(
        "Unsupported type for data_long.",
        "Datalist should be a data.frame."
      )

      # convert factors
      embed <- cpp_data_factorize(data_long)
      # store levels
      private$.info_data_long <- embed$info
      # store column names
      private$.cols_data_long <- embed$cols
      # import data.frame in rcpp
      private$.rcpp$import_data_long(embed$data)

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add a new parameter to houspepi model definitions.
    #' If the parameter already exists, it is updated.
    #'
    #' @param name `r rd_var_name_`
    #' @param value `r rd_var_value_`
    #' @param const `r rd_var_const_`
    #' @param discr `r rd_par_discr_`
    #' @param prior `r rd_par_model_`
    #' @param ... `r rd_par_args_`
    #' @param min `r rd_var_min_`
    #' @param max `r rd_var_max_`
    #' @param proposal `r rd_par_proposal_`
    #' @param proposal_scale `r rd_var_proposal_scale_`
    add_parameter = function(name = NA_character_, value = NULL,
                             const = NA, discr = FALSE,
                             prior = NA_character_, ...,
                             min = NULL, max = NULL,
                             proposal = NA_character_,
                             proposal_scale = NA_real_) {
      private$.param$add_quo(name = name, value = rlang::enquo(value),
                             const = const, discr = discr,
                             prior = prior, rlang::enquos(...),
                             min = rlang::enquo(min), max = rlang::enquo(max),
                             proposal = proposal,
                             proposal_scale = proposal_scale)
      invisible(self)
    },

    #' @description
    #' Load parameter definitions either as arguments or from predefined list.
    #'
    #' @param ... `r rd_list_args_`
    #' @param append `r rd_list_append_`
    load_parameter_list = function(..., append = FALSE) {
      private$.param$load_quosures(rlang::enquos(...), append = append)
      invisible(self)
    },

    #' @description
    #' Load parameter definitions from model template.
    #'
    #' @param append `r rd_tmpl_append_`
    load_parameter_template = function(append = FALSE) {
      if (private$.is_blank()) warning_h(
        "No template parameters available for blank model."
      )
      # retrieve the template parameters file
      else private$.param$load_template(
        file.path(private$.path(), "template_param.R"),
        private$.param$require(private$.reference()),
        isTRUE(append)
      )
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add a new augmented variable/function to houspepi model definitions.
    #' If the augmented variable already exists, it is updated.
    #'
    #' @param name `r rd_var_name_`
    #' @param value `r rd_aug_val_`
    #' @param status `r rd_aug_status_`
    #' @param const `r rd_var_const_`
    #' @param discr `r rd_aug_discr_`
    #' @param func `r rd_aug_func_`
    #' @param func_args `r rd_aug_funargs_`
    #' @param model `r rd_aug_model_`
    #' @param ... `r rd_aug_args_`
    #' @param min `r rd_aug_min_`
    #' @param max `r rd_aug_max_`
    #' @param proposal `r rd_var_proposal_`
    #' @param proposal_scale `r rd_var_proposal_scale_`
    add_augmented = function(name = NA_character_,  value = NULL,
                             status = NA_integer_, const = NA, discr = FALSE,
                             func = NULL, func_args = NULL,
                             model = NA_character_, ...,
                             min = NULL, max = NULL,
                             proposal = NA_character_,
                             proposal_scale = NA_real_) {
      private$.augm$add_quo(name = name, value = rlang::enquo(value),
                            status = status, const = const, discr = discr,
                            func = rlang::enquo(func),
                            func_args = rlang::enquos(func_args),
                            model = model, rlang::enquos(...),
                            min = rlang::enquo(min), max = rlang::enquo(max),
                            proposal = proposal,
                            proposal_scale = proposal_scale)
      invisible(self)
    },

    #' @description
    #' Load augmented variable/function definitions
    #' either as arguments or from predefined list.
    #'
    #' @param ... `r rd_list_args_`
    #' @param append `r rd_list_append_`
    load_augmented_list = function(..., append = FALSE) {
      private$.augm$load_quosures(rlang::enquos(...), append = append)
      invisible(self)
    },

    #' @description
    #' Load augmented variable/function definitions from template.
    #'
    #' @param append `r rd_tmpl_append_`
    load_augmented_template = function(append = FALSE) {
      if (private$.is_blank()) warning_h(
        "No template augmented variables available for blank model."
      )
      # retrieve the template augmented file
      else private$.augm$load_template(
        file.path(private$.path(), "template_augm.R"),
        private$.augm$require(private$.reference()),
        isTRUE(append)
      )
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add definitions of a group rate by given data feature(s).
    #' Group rates are relative to a reference class that can be specified.
    #' If no MCMC options are provided, the following configuration is assumed
    #' - `prior`: "lognormal"
    #' - `meanlog`: 0
    #' - `sdlog`: 3
    #' - `min`: 0
    #' - `max` 10
    #'
    #' @param name `r rd_var_name_`
    #' @param by `r rd_grr_by_`
    #' @param prefix `r rd_grr_prefix_`
    #' @param reference `r rd_grr_refer_`
    #' @param independent `r rd_grr_indep_`
    #' @param prior `r rd_par_model_`
    #' @param ... `r rd_par_args_`
    #' @param min `r rd_var_min_`
    #' @param max `r rd_var_max_`
    #' @param proposal `r rd_var_proposal_`
    #' @param proposal_scale `r rd_var_proposal_scale_`
    add_group_rate = function(name, by, prefix = NA_character_,
                              reference = NULL, independent = TRUE,
                              prior = NA_character_, ...,
                              min = NULL, max = NULL,
                              proposal = NA_character_,
                              proposal_scale = NA_real_) {
      private$.group_rates$add_quo(name = name, by = by, prefix = prefix,
                                   reference = reference,
                                   independent = independent,
                                   prior = prior, rlang::enquos(...),
                                   min = rlang::enquo(min),
                                   max = rlang::enquo(max),
                                   proposal = proposal,
                                   proposal_scale = proposal_scale)
      invisible(self)
    },

    #' @description
    #' Load group rate definitions either as direct arguments
    #' or from predefined list.
    #'
    #' @param ... ` r rd_list_args_`
    #' @param append `r rd_list_append_`
    load_group_rate_list = function(..., append = FALSE) {
      if (isFALSE(append)) self$drop_group_rates()
      private$.group_rates$load_quosures(rlang::enquos(...), append = append)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Add definition of infectivity profile.
    #'
    #' @param name `r rd_var_name_`
    #' @param status `r rd_aug_status_`
    #' @param discr `r rd_inf_discr_`
    #' @param ... `r rd_inf_args_`
    #' @param max `r rd_aug_max_`
    add_infec_profile = function(name, status = NA_integer_,
                                 discr = FALSE,
                                 ..., max = NULL) {
      private$.infec_profile$add_quo(name = name,
                                     status = status, discr = discr,
                                     rlang::enquos(...),
                                     max = rlang::enquo(max))
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Evaluate the log-likelihood.
    #'
    #' @param house_id `r rd_house_id_`
    #' If `NA` (default) it computes for the whole data.
    #' @param seed `r rd_seed_`
    #'
    #' @return `numeric(1)`
    evaluate_log_lik = function(house_id = NA_character_, seed = NA_integer_) {
      # internal build function
      private$.build_full(seed = seed)
      tryCatch(output <- private$.rcpp$evaluate_log_lik(house_id),
               error = function(e) private$.helper_log_lik_error(e))
      output
    },

    #' @description
    #' Run inference with Bayesian Markov Chain Monte Carlo.
    #'
    #' @param steps `r rd_smpl_steps_`
    #' @param burnin `r rd_smpl_burn_`
    #' @param thin `r rd_smpl_thin_`
    #' @param adapt `r rd_smpl_adapt_`
    #' @param chains `r rd_smpl_chain_`
    #' @param threads `r rd_smpl_thread_`
    #' @param id_draws `r rd_smpl_iddrw_`
    #' @param seed `r rd_seed_`
    #' @param progress `r rd_progress_`
    run_inference = function(steps = 200, burnin = 0,
                             thin = 1, adapt = 0,
                             chains = 1, threads = 1,
                             id_draws = NULL,
                             seed = NA_integer_,
                             progress = TRUE) {
      # internal build function
      private$.build_full(steps = steps, burnin = burnin,
                          thin = thin, adapt = adapt,
                          chains = chains, threads = threads,
                          id_draws = id_draws, seed = seed, progress = progress)
      # infer method from cpp model object and deal with errors
      tryCatch(private$.rcpp$infer(),
               error = function(e) private$.helper_log_lik_error(e))
      HousepiPosteriorDraws$new(private$.name, private$.rcpp)
    },

    #' @description
    #' Dataset with associated augmented variables for last state.
    #'
    #' @return
    #' `data.frame`
    export_augmented_data = function() {
      private$.refactor_data_wide(
        cbind(private$.rcpp$data_wide_for_simulation(),
              private$.rcpp$export_augmented_data_wide())
      )
    },

    #' @description
    #' All variables associated to given household for last state.
    #'
    #' @param house_id `r rd_house_id_`
    #' If `NA` it calls the household which caused last error, if any.
    #'
    #' @return `list`
    export_house_snapshot = function(house_id = NA) {
      # fetch house id
      if (is.na(house_id)) {
        if (is.na(private$.crash_house)) stop_h("No 'house_id' was provided.")
        house_id <- private$.crash_house
      }
      house_id <- as.character(house_id)
      # export household snapshot
      output <- private$.rcpp$export_house_snapshot(house_id)
      output$data_wide <- private$.refactor_data_wide(output$data_wide)
      # > data long
      if (all(dim(output$data_long) == 0)) output$data_long <- NULL
      else output$data_long <- private$.refactor_data_long(output$data_long)
      # > group rates
      if (all(dim(output$group_rates) == 0)) output$group_rates <- NULL
      # > infectivity profile # UPDATE
      if (length(output$infec_profile) == 0) {
        output$infec_profile <- NULL
      } else {
        lst <- output$infec_profile
        df <- data.frame(id = lst$id, discr = lst$discr, name = lst$name,
                         min = lst$min, max = lst$max)
        df$args <- lst$args
        output$infec_profile <- df
      }
      output
    },

    #' @description
    #' Run simulation to generate a household survey.
    #' It assumes the same community structure of the loaded data.
    #'
    #' @param nsim `r rd_simul_nsim_`
    #' @param seed `r rd_seed_`
    #' @param progress `r rd_progress_`
    #' @param obs_process `r rd_simul_observ_`
    #'
    #' @return `data.frame`
    run_simulation = function(nsim = 1, seed = NA_integer_,
                              progress = FALSE, obs_process = NULL) {
      # set observation process to identity if not provided
      if (is.null(obs_process)) obs_process <- function(x) x
      else if (!is.function(obs_process)) stop_h(
        "Unrecognized function passed as observation process."
      )

      # internal build function
      private$.build_full(seed = seed, progress = progress)

      # refactor static data for simulation
      static_df <- private$.rcpp$data_wide_for_simulation()
      # produce list of simulated dataframes
      # then bind simulated entries with static part
      output <- lapply(private$.rcpp$simulate(nsim), function(x) {
        obs_process(private$.refactor_data_wide(cbind(static_df, x)))
      })
      # return a single data.frame if single simulation
      if (nsim == 1) output <- output[[1]]

      output
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Drop parameter definitions.
    drop_parameters = function() {
      private$.param <- HousepiParameter$new()
      invisible(self)
    },

    #' @description
    #' Drop augmented variable definitions.
    drop_augmented = function() {
      private$.augm <- HousepiAugmented$new()
      invisible(self)
    },

    #' @description
    #' Drop group rates definitions.
    drop_group_rates = function() {
      # remove associated parameters
      if (!is.null(private$.group_rates)) {
        for (p in names(private$.group_rates$parameters)) {
          private$.param$remove(p, show_warnings = FALSE)
        }
      }
      private$.group_rates <- HousepiGroupRate$new()
      invisible(self)
    },

    #' @description
    #' Drop all definitions.
    drop_all = function() {
      self$drop_parameters()
      self$drop_augmented()
      self$drop_group_rates()
      private$.infec_profile <- HousepiInfecProfile$new()
      invisible(self)
    }
  )
)
