#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# sampler_class.R
#

############################################################################
############################################################################
############################ Individual Sampler ############################
############################################################################
############################################################################

#' @title
#' Individual Sampler
#' @description
#' FIXME: missing.
#'
#' @examples
#' data_wide <- housepi::survey
#' data_wide$sex <- ifelse(data_wide$sex == 1, "M", "F")
#'
#' bic <- IndividualSampler$new(
#'   sex = list(levels = c("M", "F")),
#'   age = list(breaks = c(0, 10, 20, 40, 60))
#' )
#' test <- data_wide[1:5, ]
#' test
#' bic$encode(test)
#' bic$decode(bic$encode(test))
#'
#' @keywords internal
IndividualSampler <- R6::R6Class(
  classname = "IndividualSampler",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),
    .enc_dict = list(),
    .feats_pdf = data.frame(),
    .expand_pdf = data.frame(),
    .num_feats = list(),

    # @description
    # A data.frame with all distinct features considered.
    .get_distinct_features = function() {
      aux <- lapply(private$.enc_dict, function(v) v$x)
      output <- expand.grid(aux)
      colnames(output) <- names(private$.enc_dict)
      output
    },

    # @description
    # Reorder according to hierarchy
    #
    # @param this_house `data.frame`\cr
    # Input encoded hosuehold data.
    #
    # @return `data.frame`
    .reorder = function(this_house) {
      # UPDATE: can I turn this into an option?
      # NOTE: do.call(order, ... requires a list
      sort_by <- ".features"

      # sort by the two 'principal' household members
      if (private$.ascending) {
        ordering <- do.call(order, this_house[sort_by])
      } else {
        ordering <- do.call(order, lapply(this_house[sort_by], function(x) -x))
      }
      # extract the features of the two 'principal' house members
      this_house[ordering, ]
    },

    # @description
    # Check internally that encoding\decoding routine.
    .check_internally = function() {
      input <- private$.get_distinct_features()
      output <- self$decode(self$encode(input))
      all(output == input)
    },

    # @description
    # Associate a marginal probability to each individual feature.
    #
    # @param data_wide `data.frame()`\cr
    # A data.frame containing the data_wide with individual information by row.
    .compute_individual_marginals = function(data_wide) {
      df <- data_wide
      df$feature <- self$encode(df)

      # update feature pdf
      feature_counts <- table(df$feature)
      private$.feats_pdf <- data.frame(
        feature = names(feature_counts),
        freqs = feature_counts / sum(feature_counts)
      )

      # update expanded pdf
      myvec <- c("feature", attributes(private$.enc_dict)$names)
      tab <- as.data.frame(table(df[myvec]))
      tab$freqs <- tab$Freq / sum(tab$Freq)
      private$.expand_pdf <- tab[setdiff(names(tab), "Freq")]

      invisible(self)
    },

    # @description
    # Generate random individuals according to feature distribution.
    #
    # @param size `integer(1)`\cr
    # Number of individuals to generate.
    .sample_individuals = function(size) {
      if ((size > 0) && (nrow(private$.feats_pdf) > 0)) {
        # sample according to observed features
        output <- self$decode(
          sample(private$.feats_pdf$feature,
                 size, replace = TRUE,
                 prob = private$.feats_pdf$freqs)
        )
      } else {
        # nothing to sample from
        output <- matrix(ncol = length(private$.enc_dict), nrow = 0)
        output <- as.data.frame(output)
        colnames(output) <- names(private$.enc_dict)
      }

      output
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field hierarchy `data.frame`\cr
    #' Features as hierchically reordered by the sampler.
    hierarchy = function(value) {
      if (missing(value)) {
        output <- private$.get_distinct_features()
        output$.features <- self$encode(output)
        return(private$.reorder(output))
      }
      stop_h("'hierarchy' is read only.")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' New object of class
    #' [IndividualSampler][housepi::IndividualSampler].
    #'
    #' @param ...
    #' Named features with respective list of definition details.
    #' @param from `list` of `list`\cr
    #' List containing encoding dictionary definitions for each feature.
    initialize = function(..., from = NULL) {
      # assign encoding dictionary
      enc_dict <- load_definition_list(..., from = from)

      # assign features distributions
      private$.feats_pdf <- data.frame()
      private$.expand_pdf <- data.frame()

      # assign number of features
      private$.num_feats <- lapply(enc_dict, function(x) 0)

      # loop over all features of encoding dictionary
      private$.enc_dict <- list()
      for (col in names(enc_dict)) {
        # accept only lists
        if (is.list(enc_dict[[col]])) {
          if ("breaks" %in% names(enc_dict[[col]])) {
            # if breaks, consider also lower than min and larger than max
            x <- enc_dict[[col]]$breaks
            private$.num_feats[[col]] <- length(x) + 2
            private$.enc_dict[[col]]$x <- x
            private$.enc_dict[[col]]$type <- "breaks"

          } else if ("levels" %in% names(enc_dict[[col]])) {
            # if levels, consider also unknown case
            x <- enc_dict[[col]]$levels
            private$.num_feats[[col]] <- length(x) + 1
            private$.enc_dict[[col]]$x <- x
            private$.enc_dict[[col]]$type <- "levels"

          } else {
            stop_h(
              sprintf("Incomplete encoding for feature \"%s\"", col),
              "Provide one between \"breaks\" ans \"levels\"."
            )
          }

        } else {
          stop_h("The encoding dictionary must contains only lists.")
        }
      }

      # internal consistency check
      if (!private$.check_internally()) stop_h(
        "Internal inconsistency, check your encoding dictionary."
      )

      invisible(self)
    },

    #' @description
    #' Encode individuals from data_wide according to encoding dictionary.
    #'
    #' @param data_wide `data.frame()`\cr
    #' A data.frame containing the data_wide with individual information by row.
    encode = function(data_wide) {
      # initialize output
      output <- rep(0, nrow(data_wide))

      # encode features for each provided column
      multiplier <- 1
      for (i in seq_along(names(private$.enc_dict))) {

        # load locally the encode dictionary
        enc <- private$.enc_dict[[i]]
        # load locally the features to be encoded
        df <- data_wide[[names(private$.enc_dict)[i]]]

        if (enc$type == "breaks") {
          feats <- findInterval(df, c(-Inf, enc$x, Inf)) - 1

        } else if (enc$type == "levels") {
          tmp <- match(df, enc$x)
          feats <- ifelse(is.na(tmp), 0, tmp)
        }

        if (i > 1) {
          multiplier <- multiplier * private$.num_feats[[i - 1]]
        }

        output <- output + feats * multiplier
      }

      output
    },

    #' @description
    #' Decode individual featrues.
    #'
    #' @param feats `numeric()`\cr
    #' Encoded features of each individual.
    decode = function(feats) {
      output <- as.data.frame(matrix(
        ncol = length(private$.enc_dict),
        nrow = length(feats)
      ))
      colnames(output) <- names(private$.enc_dict)

      for (col in colnames(output)) {
        enc <- private$.enc_dict[[col]]

        if (enc$type == "breaks") {
          output[[col]] <- c(NaN, enc$x, NaN)[1 + feats %% private$.num_feats[[col]]]

        } else if (enc$type == "levels") {
          aux_dict <- stats::setNames(enc$x, seq_along(length(enc$x)))
          output[[col]] <- sapply(
            feats %% private$.num_feats[[col]],
            function(x) ifelse(x > 0, aux_dict[[x]], NaN)
          )
        }

        feats <- floor(feats / private$.num_feats[[col]])
      }

      output
    }
  )
)

############################################################################
############################################################################
############################# Household Sampler ############################
############################################################################
############################################################################

#' @title
#' Household Sampler
#' @description
#' FIXME: missing.
#'
#' @examples
#' bhhc <- HousepiSampler$new(enc_dict)
#' bhhc$compute_marginals(data_wide)
#'
#' # sample 5 household from model
#' set.seed(291)
#' simulated <- bhhc$sample(200)
#' head(simulated)
#'
#' @docType class
#' @keywords internal
HousepiSampler <- R6::R6Class(
  classname = "HousepiSampler",
  inherit = IndividualSampler,
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    .ascending = NA,
    .single_prob = NA,
    .joint_pdf = data.frame(),
    .hhsize_pdf = data.frame(),
    .other_pdf = data.frame()
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' New object of class
    #' [HousepiSampler][housepi::HousepiSampler].
    #'
    #' @param ...
    #' Named features with respective list of definition details.
    #' @param from `list` of `list`\cr
    #' List containing encoding dictionary definitions for each feature.
    #' @param ascending `logical(1)`\cr
    #' If `TRUE` hierarchy between feature numbers is ascending.
    #' If `FALSE` (default) is descending.
    initialize = function(..., from = NULL, ascending = FALSE) {
      # ininitialize parent class IndividualSampler
      super$initialize(..., from = from)

      private$.ascending <- isTRUE(ascending)

      # default marginals
      private$.single_prob <- 0.
      private$.joint_pdf <- data.frame()
      private$.hhsize_pdf <- data.frame()
      private$.other_pdf <- data.frame()

      invisible(self)
    },

    #' @description
    #' Compute household marginal distributions according to hierarchical model.
    #'
    #' @param data_wide `data.frame()`\cr
    #' A data.frame containing the data_wide with individual information by row.
    #' @param ascending `logical(1)`\cr
    #' If `TRUE` hierarchy between feature numbers is ascending.
    #' If `FALSE` (default) is descending.
    compute_marginals = function(data_wide, ascending = NA) {
      if (!is.na(ascending)) private$.ascending <- ascending

      # extract household size
      if (!("house_id" %in% colnames(data_wide))) stop_h(
        "No household structure provided in data_wide.",
        "Missing column 'house_id'."
      )

      df <- data_wide
      hhsize <- table(df$house_id)
      df$hhsize <- hhsize[as.character(df$house_id)]

      # add a column with encoded features
      df$.features <- self$encode(df)

      #----- single individual (size 1) households -----#

      # count for household with single individuals
      private$.compute_individual_marginals(df[df$hhsize == 1, ])
      private$.single_prob <- nrow(df[df$hhsize == 1, ]) / nrow(df)

      #----- larger households -----#

      # list with all model features
      model_features <- seq(prod(unlist(private$.num_feats))) - 1

      # all pairs of house members ordered by feature
      # note : pair (X, X) are appended manually
      combinations <- c(
        utils::combn(model_features, 2, simplify = FALSE),
        lapply(model_features, function(x) c(x, x))
      )

      # initialize frequencies
      freqs <- list()
      # joint : joint compositions of the two 'principal' house members
      freqs$joint <- data.frame(
        hhmemA = sapply(combinations, `[`, 1),
        hhmemB = sapply(combinations, `[`, 2)
      )
      freqs$joint$freqs <- 0

      # other : initialize conditioned compositions of other members
      freqs$other <- freqs$joint
      for (i in model_features) {
        freqs$other[[as.character(i)]] <- 0
      }

      # size : initialize household conditioned sizes
      freqs$hhsize <- freqs$joint
      for (i in seq(2, max(df$hhsize))) freqs$hhsize[[as.character(i)]] <- 0

      # count for household with at least 2 members
      for (h in unique(df[df$hhsize > 1, ]$house_id)) {
        this_house <- private$.reorder(df[df$house_id == h, ])
        x_1 <- this_house$.features[1]
        x_2 <- this_house$.features[2]
        # swap to match combination table convention
        if (x_1 > x_2) {
          x_1 <- this_house$.features[2]
          x_2 <- this_house$.features[1]
        }

        # count 'principal' pair
        jmask <- which(freqs$joint$hhmemA == x_1 & freqs$joint$hhmemB == x_2)
        freqs$joint[jmask, "freqs"] <- freqs$joint[jmask, "freqs"] + 1

        # count sizes
        size <- this_house$hhsize[1]
        csize <- as.character(size)
        freqs$hhsize[jmask, csize] <- freqs$hhsize[jmask, csize] + 1

        # count other members
        if (size > 2) {
          for (i in seq(3, size)) {
            cfeat <- as.character(this_house$.features[i])
            freqs$other[jmask, cfeat] <- freqs$other[jmask, cfeat] + 1
          }
        }
      }

      # normalize joint
      freqs$joint$freqs <- freqs$joint$freqs / sum(freqs$joint$freqs)

      # select size count columns
      scols <- 3:length(colnames(freqs$hhsize))
      # select other count columns
      ocols <- 3:length(colnames(freqs$other))

      # normalize hhsize and other
      for (i in seq_along(nrow(freqs$joint))) {
        norm <- max(1, sum(freqs$other[i, ocols]))
        freqs$other[i, ocols] <- freqs$other[i, ocols] / norm

        norm <- max(1, sum(freqs$hhsize[i, scols]))
        freqs$hhsize[i, scols] <- freqs$hhsize[i, scols] / norm
      }

      # assign marginals
      private$.joint_pdf <- freqs$joint
      private$.other_pdf <- freqs$other
      private$.hhsize_pdf <- freqs$hhsize

      invisible(self)
    },

    #' @description
    #' Sample households according to feature distribution.
    #'
    #' @param size `integer(1)`\cr
    #' Number of households to sample.
    sample = function(size = 1) {

      if (size < 1) stop_h(
        "Sample size must be greater than 0."
      )

      #----- single individual (size 1) households -----#

      # number of single-member households
      num_single <- stats::rbinom(1, size = size, p = private$.single_prob)
      hh_single <- NULL
      if (num_single > 0) {
        hh_single <- private$.sample_individuals(num_single)
        hh_single$house_id <- seq_len(nrow(hh_single))
      }

      #----- larger households -----#

      # number of households with size >= 2
      num_greater <- size - num_single
      hh_greater <- NULL
      if (num_greater > 0) {
        # start by sampling the two 'principal' members
        idx <- sample(seq_len(nrow(private$.joint_pdf)),
                      num_greater,
                      replace = TRUE,
                      prob = private$.joint_pdf$freqs)

        # range and respective features of hhsize
        scols <- 4:length(private$.hhsize_pdf)
        sizelist <- as.integer(colnames(private$.hhsize_pdf[scols]))
        # range and respective features of other members
        ocols <- 4:length(private$.other_pdf)
        otherlist <- as.integer(colnames(private$.other_pdf[ocols]))

        # initialize encoded dataframe
        df <- data.frame()

        for (h in seq_len(num_greater)) {
          # add 'principal' pair (ascending by convention)
          princ_feats <- as.integer(private$.joint_pdf[idx[h], c(1, 2)])

          # sample house size
          logic_a <- private$.joint_pdf$hhmemA == princ_feats[1]
          logic_b <- private$.joint_pdf$hhmemB == princ_feats[2]
          jmask <- which(logic_a & logic_b)
          hhsize <- sample(sizelist, 1,
                           prob = private$.hhsize_pdf[jmask, scols])

          # sample other members if any
          other_feats <- c()
          if (hhsize > 2) {
            for (i in seq_len(hhsize - 2)) {
              y <- sample(otherlist, 1, prob = private$.other_pdf[jmask, ocols])
              other_feats <- c(other_feats, y)
            }
          }

          if (!private$.ascending) princ_feats <- princ_feats[2:1]
          # append new household
          df <- rbind(
            df,
            data.frame(house_id = h, .features = c(princ_feats, other_feats))
          )
        }

        # create data_wide of sizes > 1 by decoding generated features
        hh_greater <- self$decode(df$.features)
        hh_greater$house_id <- df$house_id + num_single
      }

      # bind together
      output <- rbind(hh_single, hh_greater)
      # format house id
      house_id_formatter <- function(x, size) {
        sprintf(paste0("hh_%0", floor(log10(size) + 1), "d"), x)
      }
      output$house_id <- house_id_formatter(output$house_id, size)
      # associate individual index
      id_formatter <- function(x) {
        sprintf(paste0("%s_%0", floor(log10(max(sizelist)) + 1), "d"),
                x, seq_along(x))
      }
      output$id <- ave(output$house_id, output$house_id, FUN = id_formatter)
      # re-order columns
      output <- output[c("id", "house_id", names(private$.enc_dict))]

      return(output)
    }
  )
)

#' @title Household Sampler
#' @description
#' The `HousepiSampler` class allows you to define and sample synthetic household datasets
#' based on hierarchical models of covariates. It supports encoding and decoding of features,
#' computation of marginal distributions, and sampling of households.
#'
#' @examples
#' # Define a sampler with covariates
#' hhsampler <- new_household_sampler(
#'   sex = list(levels = c("M", "F")),
#'   age = list(breaks = c(0, 15, 25, 50))
#' )
#'
#' # Example household dataset
#' eg_house <- data.frame(
#'   house_id = "HH_1",
#'   age = c(43, 39, 12, 7),
#'   sex = c("M", "F", "F", "M")
#' )
#'
#' # Encode and decode the dataset
#' encoded <- hhsampler$encode(eg_house)
#' decoded <- hhsampler$decode(encoded)
#'
#' # Compute marginals from a survey dataset
#' hhsampler$compute_marginals(housepi::survey)
#'
#' # Sample 5 households
#' sampled_households <- hhsampler$sample(5)
#' head(sampled_households)
#'
#' @export
new_household_sampler <- function(...) {
  HousepiSampler$new(...)
}

#' @title Compute Marginals
#' @description
#' Computes marginal distributions for household features based on a survey dataset.
#'
#' @param object An instance of the `HousepiSampler` class.
#' @param data_wide A `data.frame` containing individual-level data with household IDs.
#'
#' @examples
#' hhsampler <- new_household_sampler(
#'   sex = list(levels = c("M", "F")),
#'   age = list(breaks = c(0, 15, 25, 50))
#'   ) |>
#'   compute_household_marginals(housepi::survey)
#' 
#' @export
compute_household_marginals <- function(object, data_wide) {
  stopifnot(inherits(object, "HousepiSampler"))
  object$clone(deep = TRUE)$compute_marginals(data_wide)
}

#' @title Sample Households
#' @description
#' Samples synthetic households based on the marginal distributions computed by the `HousepiSampler`.
#'
#' @param object An instance of the `HousepiSampler` class.
#' @param size An integer specifying the number of households to sample.
#'
#' @return A `data.frame` containing the sampled households.
#' @examples
#' sampled_households <- new_household_sampler(
#'   sex = list(levels = c("M", "F")),
#'   age = list(breaks = c(0, 15, 25, 50))
#'   ) |>
#'   compute_household_marginals(housepi::survey) |> 
#'   sample_households(5)
#' head(sampled_households)
#'
#' @export
sample_households <- function(object, size, ...) {
  stopifnot(inherits(object, "HousepiSampler"))
  object$sample(size)
}
