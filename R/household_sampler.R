library(dplyr)

###########################################################
##########  Class : Brew Individual Composition  ##########
###########################################################
IndividualSampler <- setRefClass(
  #===== NAME =====#
  "IndividualSampler",
  
  #===== ATTRIBUTES =====#
  fields = list(
    m_enc_dict = "list",
    m_feats_pdf = "data.frame",
    m_expand_pdf = "data.frame",
    m_num_feats = "list"
  ),
  
  #===== METHODS =====#
  methods = list(
    #===== Method : initialize =====#
    initialize = function(mydict) {
      # assign encoding dictionary
      m_enc_dict <<- mydict
      
      # assign features distributions
      m_feats_pdf <<- base::data.frame()
      m_expand_pdf <<- base::data.frame()
      
      # assign number of features
      m_num_feats <<- base::lapply(m_enc_dict, function(x) 0)
      
      # loop over all features of encoding dictionary
      for (col in base::names(m_enc_dict)) {
        # accept only lists
        if (base::is.list(m_enc_dict[[col]])) {
          # check if type and vector are provided
          # UPDATE: find a better name for x maybe
          if (base::all(c("x", "type") %in% base::names(m_enc_dict[[col]]))) {
            if (m_enc_dict[[col]]$type == "bin") {
              # if bin consider lower than min and larger than max
              m_num_feats[[col]] <<- base::length(m_enc_dict[[col]]$x) + 2
            } else if (m_enc_dict[[col]]$type == "cat") {
              # if cat (i. e. categorical) consider unknown case
              m_num_feats[[col]] <<- base::length(m_enc_dict[[col]]$x) + 1
            } else{
              stop("Unrecognized type `", m_enc_dict[[col]]$type, "` for ", col, ".")
            }
          } else {
            stop("Provide both `x` and `type` for feature ", col, ".")
          }
        } else {
          stop("The encoding dictionary must contains only vectors.")
        }
      }
      
      # internal consistency check
      if (!.internal.check()) {
        stop("Internal inconsistency, check your encoding dictionary.")
      }
    },
    
    #===== Method : encode =====#
    encode = function(database) {
      # initialize output
      output <- base::rep(0, base::nrow(database))
      
      # encode features for each provided column
      multiplier <- 1
      for (i in base::seq_along(base::names(m_enc_dict))) {
        
        # load locally the encode dictionary
        enc <- m_enc_dict[[i]]
        # load locally the features to be encoded
        df <- database[[base::names(m_enc_dict)[i]]]
        
        if (enc$type == "bin") {          
          feats <- base::findInterval(df, c(-Inf, enc$x, Inf)) - 1
        }
        else if (enc$type == "cat") {
          tmp <- base::match(df, enc$x)
          feats <- base::ifelse(base::is.na(tmp), 0, tmp)
        }
        
        if (i > 1) {
          multiplier <- multiplier * m_num_feats[[i - 1]]
        }
        
        output <- output + feats * multiplier
      }
      
      return(output)
    },
    
    #===== Method : decode =====#
    decode = function(feats) {
      output <- base::as.data.frame(
        base::matrix(ncol = base::length(m_enc_dict), nrow = base::length(feats)))
      colnames(output) <- base::names(m_enc_dict)
      
      for (col in base::colnames(output)) {
        enc <- m_enc_dict[[col]]
        
        if (enc$type == "bin") {
          output[[col]] <- c(NaN, enc$x, NaN)[1 + feats %% m_num_feats[[col]]]
          
        } else if (enc$type == "cat") {
          aux_dict <- stats::setNames(enc$x, base::seq(base::length(enc$x)))
          output[[col]] <- base::sapply(
            feats %% m_num_feats[[col]], 
            function(x) base::ifelse(x > 0, aux_dict[[x]], NaN)
            )
        }
        
        feats <- base::floor(feats / m_num_feats[[col]])
      }
      
      return(output)
    },
    
    #===== Method (developer) : internal check =====#
    .internal.check = function() {
      aux <- base::lapply(m_enc_dict, function(v) v$x)
      input <- base::expand.grid(aux)
      base::colnames(input) <- base::names(m_enc_dict)
      
      output <- decode(encode(input))
      
      return(base::all(output == input))
    },
    
    #===== Method : compute individual marginals =====#
    compute_individual_marginals = function(database) {
      feats <- database %>%
        mutate(feature = encode(.))
      
      # update feature pdf
      m_feats_pdf <<- feats %>%
        group_by(feature) %>%
        tally() %>%
        mutate(freqs = n / base::sum(n)) %>%
        select(!n)
      
      # update expanded pdf
      myvec <- c("feature", base::attributes(m_enc_dict)$names)
      m_expand_pdf <<- feats %>%
        group_by(across(all_of(myvec))) %>%
        tally() %>%
        mutate(freqs = n / base::sum(n)) %>%
        select(!n)
    },
    
    #===== Method : generate individuals =====#
    generate_individuals = function(n) {
      if(base::nrow(m_feats_pdf) > 0) {
        # sample according to observed features
        output <- decode(base::sample(m_feats_pdf$feature,
                                      n,
                                      replace = TRUE,
                                      prob = m_feats_pdf$freqs))
      }
      else {
        # nothing to sample from
        output <- base::as.data.frame(
          base::matrix(ncol = base::length(m_enc_dict), nrow = 0))
        base::colnames(output) <- base::names(m_enc_dict)
      }
      
      return(output)
    }
  )
)


##########################################################
##########  Class : Brew Household Composition  ##########
##########################################################

#' FIXME
#' @export
HouseholdSampler  <- setRefClass(
  #===== NAME =====#
  "HouseholdSampler",
  
  #==== INHERITANCE =====#
  contains = "IndividualSampler",
  
  
  #===== ATTRIBUTES =====#
  fields = list(
    m_single_prob = "numeric",
    m_joint_pdf = "data.frame",
    m_hhsize_pdf = "data.frame",
    m_other_pdf = "data.frame"
  ),
  
  #===== METHODS =====#
  methods = list(
    #----- initialize -----#
    initialize = function(m_enc_dict) {
      # ininitialize parent class IndividualSampler
      callSuper(m_enc_dict)
      
      # default marginals
      m_single_prob <<- 0.
      m_joint_pdf <<- base::data.frame()
      m_hhsize_pdf <<- base::data.frame()
      m_other_pdf <<- base::data.frame()
    },
    
    #----- compute household marginals -----#
    compute_household_marginals = function(database,
                                           sortby = "features",
                                           ascending = FALSE) {
      # extract household size
      df <- database %>%
        group_by(house_index) %>%
        mutate(hhsize = n()) %>%
        ungroup()
      
      # add a column with encoded features
      df$features <- encode(df)
      
      #-------------------------------------------------#
      #----- single individual (size 1) households -----#
      #-------------------------------------------------#
      
      # count for household with single individuals
      compute_individual_marginals(df[df$hhsize == 1, ])
      m_single_prob <<- base::nrow(df[df$hhsize == 1, ]) / base::nrow(df)
      
      #-----------------------------#
      #----- larger households -----#
      #-----------------------------#
      
      # list with all model features
      model_features <- base::seq(base::prod(base::unlist(m_num_feats))) - 1
      
      # frequencies list
      freqs <- base::list()
      
      # all pairs of house members ordered by feature
      # note : pair (X, X) are appended manually
      combinations <- c(
        utils::combn(model_features, 2, simplify = FALSE),
        base::lapply(model_features, function(x) c(x, x))
        )
      
      # joint compositions of the two `main` house members
      joint_tuple <- base::data.frame(
        hhmemA = base::sapply(combinations, `[`, 1),
        hhmemB = base::sapply(combinations, `[`, 2)
        )
      joint_tuple <- joint_tuple %>%
        arrange(hhmemA, hhmemB)
      
      # initialize joint freqs
      freqs$joint <- joint_tuple %>%
        mutate(freqs = 0)
      
      # other : initialize conditioned compositions of other members
      freqs$other <- joint_tuple
      for (i in model_features) {
        freqs$other[[base::as.character(i)]] <- 0
      }
      
      # size : initialize household conditioned sizes
      freqs$hhsize <- joint_tuple
      for (i in base::seq(2, base::max(df$hhsize))) {
        freqs$hhsize[[base::as.character(i)]] <- 0
      }      
      
      # count for household with at least 2 members
      for (h in base::unique(df[df$hhsize > 1, ][["house_index"]])) {
        this_house <- df %>%
          filter(house_index == h)
        
        # sort by the two `main` household members
        if (ascending == TRUE) {
          this_house <- this_house %>%
            arrange_at(sortby)
          # extract the features of the two `main` house members
          x_1 <- this_house$features[1]
          x_2 <- this_house$features[2]
        } else {
          this_house <- this_house %>%
            arrange_at(sortby, desc)
          # extract the features of the two `main` house members
          x_1 <- this_house$features[2]
          x_2 <- this_house$features[1]
        }
        
        jmask <- base::which(joint_tuple$hhmemA == x_1 &
                               joint_tuple$hhmemB == x_2)
        
        # count `main` pair
        freqs$joint[jmask, "freqs"] <- freqs$joint[jmask, "freqs"] + 1
        
        # count sizes
        size <- this_house$hhsize[1]
        csize <- base::as.character(size)
        freqs$hhsize[jmask, csize] <- freqs$hhsize[jmask, csize] + 1
        
        # count other members
        if (size > 2) {
          for (i in base::seq(3, size)) {
            cfeat <- base::as.character(this_house$features[i])
            freqs$other[jmask, cfeat] <- freqs$other[jmask, cfeat] + 1
          }
        }
      }
      
      # normalize joint
      freqs$joint$freqs <- freqs$joint$freqs / base::sum(freqs$joint$freqs)
      
      # select size count columns
      scols <- 3:base::length(base::colnames(freqs$hhsize))
      # select other count columns
      ocols <- 3:base::length(base::colnames(freqs$other))
      
      # normalize hhsize and other
      for (i in base::seq(base::nrow(freqs$joint))) {
        norm <- base::max(1, base::sum(freqs$other[i, ocols]))
        freqs$other[i, ocols] <- freqs$other[i, ocols] / norm
        
        norm <- base::max(1, base::sum(freqs$hhsize[i, scols]))
        freqs$hhsize[i, scols] <- freqs$hhsize[i, scols] / norm
      }
      
      # assign marginals
      m_joint_pdf <<- freqs$joint
      m_other_pdf <<- freqs$other
      m_hhsize_pdf <<- freqs$hhsize
    },
    
    #----- generate households -----#
    generate_households = function(n) {
      
      num_single <- stats::rbinom(1, size = n, p = m_single_prob)
      num_greater <- n - num_single
      
      # start by sampling the two `main` members
      idx <- base::sample(base::seq(base::nrow(m_joint_pdf)),
                          num_greater,
                          replace = TRUE,
                          prob = m_joint_pdf$freqs)
      
      # range and respective features of hhsize
      scols <- 3:base::length(m_hhsize_pdf)
      sizelist <- base::as.integer(colnames(m_hhsize_pdf[scols]))
      # range and respective features of other members
      ocols <- 3:base::length(m_other_pdf)
      otherlist <- base::as.integer(colnames(m_other_pdf[ocols]))
      
      # initialize encoded dataframe
      df <- base::data.frame()
      
      for (h in base::seq_len(num_greater)) {
        # add `main` pair
        df <- base::rbind(df, 
                          base::data.frame(features = m_joint_pdf[idx[h], 1],
                                           house_index = h))
        df <- base::rbind(df,
                          base::data.frame(features = m_joint_pdf[idx[h], 2],
                                           house_index = h))
        
        # sample house size
        jmask <- base::which(m_joint_pdf$hhmemA == m_joint_pdf[idx[h], 1] &
                               m_joint_pdf$hhmemB == m_joint_pdf[idx[h], 2])
        size <- base::sample(sizelist, 1, prob = m_hhsize_pdf[jmask, scols])
        
        # sample other members if any
        if (size > 2) {
          for (i in base::seq_len(size - 2)) {
            y <- base::sample(otherlist, 1, prob = m_other_pdf[jmask, ocols])
            df <- base::rbind(df, 
                              base::data.frame(features = y, 
                                               house_index = h))
          }
        }
      }
      
      # create database of sizes == 1
      hh_single <- generate_individuals(num_single) %>%
        mutate(house_index = base::seq(n()))
      
      # create database of sizes > 1 by decoding generated features
      hh_greater <- decode(df$features) %>%
        mutate(house_index = df$house_index + num_single)
      
      output <- base::rbind(hh_single, hh_greater) %>%
        mutate(index = base::seq(n()))
      
      return(output)
    }
  )
)

#==== TESTING ====#

# model according to group.classes
# enc_dict <- list(
#   sex = list(x = c("M", "F"), type = "cat"), 
#   age = list(x = c(0, 10, 20, 40, 60), type = "bin"))
# 
# database <- read.csv("database.csv")
# database$sex <- ifelse(database$sex == 1, "F", "M")

# bic <- IndividualSampler$new(enc_dict)
# test <- database[1:5,]
# print("test: ")
# print(test)
# print("encode: ")
# print(bic$encode(test))
# print("decode: ") 
# print(bc$decode(bic$encode(test)))
# 
# bic$compute_individual_marginals(database)

#-----

# bhhc <- HouseholdSampler$new(enc_dict)
# bhhc$compute_household_marginals(database)
# 
# # sample 5 household from model
# set.seed(291)
# simulated <- bhhc$generate_households(200)
