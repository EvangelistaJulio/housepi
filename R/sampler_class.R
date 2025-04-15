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

#' @title
#' Individual Sampler R6 Class
#' @description
#' FIXME: missing.
#'
#' @examples
#' dataset <- housepi::survey
#' dataset$sex <- ifelse(dataset$sex == 1, "M", "F")
#'
#' # model according to group.classes
#' enc_dict <- list(
#'   sex = list(x = c("M", "F"), type = "cat"),
#'   age = list(x = c(0, 10, 20, 40, 60), type = "bin")
#' )
#'
#' bic <- IndividualSampler$new(enc_dict)
#' test <- dataset[1:5, ]
#' test
#' bic$encode(test)
#' bic$decode(bic$encode(test))
#'
#' @export IndividualSampler
#' @docType class
IndividualSampler <- R6::R6Class(
  classname = "IndividualSampler",
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################
  ############################################################################

  private = list(
    .enc_dict = list(),
    .feats_pdf = data.frame(),
    .expand_pdf = data.frame(),
    .num_feats = list(),

    # @description
    # Check internally that encoding\decoding routine.
    .internal_check = function() {
      aux <- lapply(private$.enc_dict, function(v) v$x)
      input <- expand.grid(aux)
      colnames(input) <- names(private$.enc_dict)

      output <- self$decode(self$encode(input))

      return(all(output == input))
    }
  ),

  ############################################################################
  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################
  ############################################################################

  public = list(
    #' @description
    #' New instance of class
    #' [IndividualSampler][housepi::IndividualSampler].
    #'
    #' @param enc_dict `list() of list()`\cr
    #' List containing encoding dictionary definitions (FIXME).
    initialize = function(enc_dict) {
      # assign encoding dictionary
      private$.enc_dict <- enc_dict

      # assign features distributions
      private$.feats_pdf <- data.frame()
      private$.expand_pdf <- data.frame()

      # assign number of features
      private$.num_feats <- lapply(private$.enc_dict, function(x) 0)

      # loop over all features of encoding dictionary
      for (col in names(private$.enc_dict)) {
        # accept only lists
        if (is.list(private$.enc_dict[[col]])) {
          # check if type and vector are provided
          # UPDATE: find a better name for x maybe
          if (all(c("x", "type") %in% names(private$.enc_dict[[col]]))) {
            if (private$.enc_dict[[col]]$type == "bin") {
              # if bin consider lower than min and larger than max
              private$.num_feats[[col]] <- length(private$.enc_dict[[col]]$x) + 2
            } else if (private$.enc_dict[[col]]$type == "cat") {
              # if cat (i. e. categorical) consider unknown case
              private$.num_feats[[col]] <- length(private$.enc_dict[[col]]$x) + 1
            } else {
              stop_h(
                sprintf("Unrecognized type '%s' for %s.",
                        private$.enc_dict[[col]]$type, col)
              )
            }
          } else {
            stop_h(
              sprintf("Provide both 'x' and 'type' for feature %s.", col)
            )
          }
        } else {
          stop_h("The encoding dictionary must contains only vectors.")
        }
      }

      # internal consistency check
      if (!private$.internal_check()) {
        stop_h("Internal inconsistency, check your encoding dictionary.")
      }
    },

    #' @description
    #' Encode individuals from dataset according to encoding dictionary.
    #'
    #' @param dataset `data.frame()`\cr
    #' A data.frame containing the dataset with individual information by row.
    encode = function(dataset) {
      # initialize output
      output <- rep(0, nrow(dataset))

      # encode features for each provided column
      multiplier <- 1
      for (i in seq_along(names(private$.enc_dict))) {

        # load locally the encode dictionary
        enc <- private$.enc_dict[[i]]
        # load locally the features to be encoded
        df <- dataset[[names(private$.enc_dict)[i]]]

        if (enc$type == "bin") {
          feats <- findInterval(df, c(-Inf, enc$x, Inf)) - 1

        } else if (enc$type == "cat") {
          tmp <- match(df, enc$x)
          feats <- ifelse(is.na(tmp), 0, tmp)
        }

        if (i > 1) {
          multiplier <- multiplier * private$.num_feats[[i - 1]]
        }

        output <- output + feats * multiplier
      }

      return(output)
    },

    #' @description
    #' Decode individual featrues.
    #'
    #' @param feats `numeric()`\cr
    #' Encoded features of each individual.
    decode = function(feats) {
      output <- as.data.frame(
        matrix(ncol = length(private$.enc_dict),
               nrow = length(feats))
      )
      colnames(output) <- names(private$.enc_dict)

      for (col in colnames(output)) {
        enc <- private$.enc_dict[[col]]

        if (enc$type == "bin") {
          output[[col]] <- c(NaN, enc$x, NaN)[1 + feats %% private$.num_feats[[col]]]

        } else if (enc$type == "cat") {
          aux_dict <- stats::setNames(enc$x, seq_along(length(enc$x)))
          output[[col]] <- sapply(
            feats %% private$.num_feats[[col]],
            function(x) ifelse(x > 0, aux_dict[[x]], NaN)
          )
        }

        feats <- floor(feats / private$.num_feats[[col]])
      }

      return(output)
    },

    #' @description
    #' Associate a marginal probability to each individual feature.
    #'
    #' @param dataset `data.frame()`\cr
    #' A data.frame containing the dataset with individual information by row.
    compute_individual_marginals = function(dataset) {
      df <- dataset
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

    #' @description
    #' Generate random individuals according to feature distribution.
    #'
    #' @param sample_size `integer(1)`\cr
    #' Number of individuals to generate.
    generate_individuals = function(sample_size) {
      if ((sample_size > 0) && (nrow(private$.feats_pdf) > 0)) {
        # sample according to observed features
        output <- self$decode(
          sample(private$.feats_pdf$feature,
                 sample_size, replace = TRUE,
                 prob = private$.feats_pdf$freqs)
        )
      } else {
        # nothing to sample from
        output <- matrix(ncol = length(private$.enc_dict), nrow = 0)
        output <- as.data.frame(output)
        colnames(output) <- names(private$.enc_dict)
      }

      return(output)
    }
  )
)

############################################################################
############################################################################
############################################################################

#' @title
#' Household Sampler R6 Class
#' @description
#' FIXME: missing.
#'
#' @examples
#' bhhc <- HouseholdSampler$new(enc_dict)
#' bhhc$compute_household_marginals(dataset)
#'
#' # sample 5 household from model
#' set.seed(291)
#' simulated <- bhhc$generate_households(200)
#' head(simulated)
#'
#' @export HouseholdSampler
#' @docType class
HouseholdSampler <- R6::R6Class(
  classname = "HouseholdSampler",
  inherit = IndividualSampler,
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################
  ############################################################################

  private = list(
    .single_prob = NA,
    .joint_pdf = data.frame(),
    .hhsize_pdf = data.frame(),
    .other_pdf = data.frame()
  ),

  ############################################################################
  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################
  ############################################################################

  public = list(
    #' @description
    #' New instance of class
    #' [HouseholdSampler][housepi::HouseholdSampler].
    #'
    #' @param enc_dict `list() of list()`\cr
    #' List containing encoding dictionary definitions (FIXME).
    initialize = function(enc_dict) {
      # ininitialize parent class IndividualSampler
      super$initialize(enc_dict)

      # default marginals
      private$.single_prob <- 0.
      private$.joint_pdf <- data.frame()
      private$.hhsize_pdf <- data.frame()
      private$.other_pdf <- data.frame()
    },

    #' @description
    #' Compute household marginal distributions according to hierarchical model.
    #'
    #' @param dataset `data.frame()`\cr
    #' A data.frame containing the dataset with individual information by row.
    #' @param ascending `logical(1)`\cr
    #' If `TRUE` hierarchy between feature numbers is ascending.
    #' Default is `FALSE`.
    compute_household_marginals = function(dataset,
                                           ascending = FALSE) {
      # UPDATE: can I turn this into an option?
      sort_by <- "features"

      # extract household size
      hhsize <- table(dataset$house_id)
      df <- dataset
      df$hhsize <- hhsize[as.character(df$house_id)]

      # add a column with encoded features
      df$features <- self$encode(df)

      #-------------------------------------------------#
      #----- single individual (size 1) households -----#
      #-------------------------------------------------#

      # count for household with single individuals
      self$compute_individual_marginals(df[df$hhsize == 1, ])
      private$.single_prob <- nrow(df[df$hhsize == 1, ]) / nrow(df)

      #-----------------------------#
      #----- larger households -----#
      #-----------------------------#

      # list with all model features
      model_features <- seq(prod(unlist(private$.num_feats))) - 1

      # frequencies list
      freqs <- list()

      # all pairs of house members ordered by feature
      # note : pair (X, X) are appended manually
      combinations <- c(
        utils::combn(model_features, 2, simplify = FALSE),
        lapply(model_features, function(x) c(x, x))
      )

      # joint compositions of the two `main` house members
      joint_tuple <- data.frame(
        hhmemA = sapply(combinations, `[`, 1),
        hhmemB = sapply(combinations, `[`, 2)
      )

      # initialize joint freqs
      freqs$joint <- joint_tuple
      freqs$joint$freqs <- 0

      # other : initialize conditioned compositions of other members
      freqs$other <- joint_tuple
      for (i in model_features) {
        freqs$other[[as.character(i)]] <- 0
      }

      # size : initialize household conditioned sizes
      freqs$hhsize <- joint_tuple
      for (i in seq(2, max(df$hhsize))) {
        freqs$hhsize[[as.character(i)]] <- 0
      }

      # count for household with at least 2 members
      for (h in unique(df[df$hhsize > 1, ]$house_id)) {
        this_house <- df[df$house_id == h, ]

        # sort by the two `main` household members
        if (ascending == TRUE) {
          this_house <- this_house[do.call(order, this_house[sort_by]), ]

          # extract the features of the two `main` house members
          x_1 <- this_house$features[1]
          x_2 <- this_house$features[2]
        } else {
          this_house <- this_house[do.call(order, lapply(this_house[sort_by], function(x) -x)), ]

          # extract the features of the two `main` house members
          x_1 <- this_house$features[2]
          x_2 <- this_house$features[1]
        }

        jmask <- which(joint_tuple$hhmemA == x_1 & joint_tuple$hhmemB == x_2)

        # count `main` pair
        freqs$joint[jmask, "freqs"] <- freqs$joint[jmask, "freqs"] + 1

        # count sizes
        size <- this_house$hhsize[1]
        csize <- as.character(size)
        freqs$hhsize[jmask, csize] <- freqs$hhsize[jmask, csize] + 1

        # count other members
        if (size > 2) {
          for (i in seq(3, size)) {
            cfeat <- as.character(this_house$features[i])
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
    },

    #' @description
    #' Generate random households according to feature distribution.
    #'
    #' @param sample_size `integer(1)`\cr
    #' Number of households to generate.
    generate_households = function(sample_size) {

      # number of single-member households
      num_single <- stats::rbinom(1, size = sample_size, p = private$.single_prob)
      if (num_single > 0) {
        hh_single <- self$generate_individuals(num_single)
        hh_single$house_id <- seq_len(nrow(hh_single))
      } else {
        hh_single <- NULL
      }

      # number of households with size >= 2
      num_greater <- sample_size - num_single
      if (num_greater > 0) {
        # start by sampling the two `main` members
        idx <- sample(seq_len(nrow(private$.joint_pdf)),
                      num_greater,
                      replace = TRUE,
                      prob = private$.joint_pdf$freqs)

        # range and respective features of hhsize
        scols <- 3:length(private$.hhsize_pdf)
        sizelist <- as.integer(colnames(private$.hhsize_pdf[scols]))
        # range and respective features of other members
        ocols <- 3:length(private$.other_pdf)
        otherlist <- as.integer(colnames(private$.other_pdf[ocols]))

        # initialize encoded dataframe
        df <- data.frame()

        for (h in seq_len(num_greater)) {
          # add `main` pair
          df <- rbind(df, data.frame(features = private$.joint_pdf[idx[h], 1],
                                     house_id = h))
          df <- rbind(df, data.frame(features = private$.joint_pdf[idx[h], 2],
                                     house_id = h))

          # sample house size
          jmask <- which(private$.joint_pdf$hhmemA == private$.joint_pdf[idx[h], 1] &
                           private$.joint_pdf$hhmemB == private$.joint_pdf[idx[h], 2])
          size <- sample(sizelist, 1, prob = private$.hhsize_pdf[jmask, scols])

          # sample other members if any
          if (size > 2) {
            for (i in seq_len(size - 2)) {
              y <- sample(otherlist, 1, prob = private$.other_pdf[jmask, ocols])
              df <- rbind(df, data.frame(features = y,
                                         house_id = h))
            }
          }
        }

        # create dataset of sizes > 1 by decoding generated features
        hh_greater <- self$decode(df$features)
        hh_greater$house_id <- sprintf("HH_%s", df$house_id + num_single)

      } else {
        hh_greater <- NULL
      }

      # bind together
      output <- rbind(hh_single, hh_greater)
      if (nrow(output) > 0) {
        output$id <- seq_len(nrow(output))
      }

      return(output)
    }
  )
)
