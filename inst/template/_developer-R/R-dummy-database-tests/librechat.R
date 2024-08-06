
library(dplyr)
library(tidyr)

Pre_Process_Household_Data <- setRefClass(
  "Pre_Process_Household_Data",
  fields = list(
    encoding_dict = "list",
    num_features = "list"
  ),
  methods = list(
    initialize = function(encoding_dict) {
      .self$encoding_dict <<- encoding_dict
      .self$num_features <<- lapply(encoding_dict, function(x) 0)
      if (!.self$check_encoder_decoder()) {
        stop("Internal inconsistency, check `encoding_dict`.")
      }
    },
    
    size = function() {
      return(prod(unlist(.self$num_features)))
    },
    
    encode = function(househ_df) {
      feats <- list()
      keys_to_encode <- names(.self$encoding_dict)
      output <- rep(0, nrow(househ_df))
      
      for (i in seq_along(keys_to_encode)) {
        col <- keys_to_encode[i]
        
        if (is.list(.self$encoding_dict[[col]])) {
          feats[[col]] <- as.integer(factor(househ_df[[col]], levels = unlist(.self$encoding_dict[[col]]))) - 1
          .self$num_features[[col]] <- length(.self$encoding_dict[[col]]) - 1
        } else if (is.vector(.self$encoding_dict[[col]])) {
          feats[[col]] <- findInterval(househ_df[[col]], .self$encoding_dict[[col]]) - 1
          .self$num_features[[col]] <- length(.self$encoding_dict[[col]]) - 1
        } else {
          stop("encoding_dict values must be either list or vector.")
        }
        
        multiplier <- if (i == 1) 1 else prod(unlist(.self$num_features)[1:(i-1)])
        output <- output + feats[[col]] * multiplier
      }
      
      return(as.integer(output))
    },
    
    decode = function(features) {
      keys_to_decode <- names(.self$encoding_dict)
      output <- as.data.frame(matrix(ncol = length(keys_to_decode), nrow = length(features)))
      colnames(output) <- keys_to_decode
      
      for (col in keys_to_decode) {
        if (is.list(.self$encoding_dict[[col]])) {
          aux_dict <- setNames(names(.self$encoding_dict[[col]]), .self$encoding_dict[[col]])
          output[[col]] <- sapply(features %% .self$num_features[[col]], function(x) aux_dict[[as.character(x)]])
        } else if (is.vector(.self$encoding_dict[[col]])) {
          output[[col]] <- .self$encoding_dict[[col]][features %% .self$num_features[[col]] + 1]
        } else {
          stop("encoding_dict values must be either list or vector.")
        }
        features <- floor(features / .self$num_features[[col]])
      }
      
      return(output)
    },
    
    check_encoder_decoder = function() {
      aux <- lapply(.self$encoding_dict, function(v) if (is.vector(v)) v[-length(v)] else names(v))
      input <- expand.grid(aux)
      colnames(input) <- names(.self$encoding_dict)
      output <- .self$decode(.self$encode(input))
      
      return(all(output == input))
    }
  )
)

Nucles_Model <- setRefClass(
  "Nucles_Model",
  contains = "Pre_Process_Household_Data",
  fields = list(
    joint_pdf = "ANY",
    other_pdf = "ANY",
    house_size_pdf = "ANY"
  ),
  methods = list(
    initialize = function(encoding_dict, marginals = NULL) {
      callSuper(encoding_dict)
      .self$assign_marginals(marginals)
    },
    
    assign_marginals = function(marginals) {
      if (!is.null(marginals)) {
        .self$joint_pdf <<- marginals[['joint']]
        .self$other_pdf <<- marginals[['other']]
        .self$house_size_pdf <<- marginals[['hh_size']]
      } else {
        .self$joint_pdf <<- NULL
        .self$other_pdf <<- NULL
        .self$house_size_pdf <<- NULL
      }
    },
    
    generate_encoded_households = function(num_househ, seed) {
      set.seed(seed)
      results <- sample(.self$joint_pdf$Var1, num_househ, replace = TRUE, prob = .self$joint_pdf$Freq)
      
      for (h in seq_len(num_househ)) {
        x_1 <- results[[h]][1]
        x_2 <- results[[h]][2]
        
        size <- sample(.self$house_size_pdf %>% filter(Var1 == x_1 & Var2 == x_2)$Var3,
                       1,
                       prob = .self$house_size_pdf %>% filter(Var1 == x_1 & Var2 == x_2)$Freq)
        
        if (size > 2) {
          for (i in seq_len(size - 2)) {
            y <- sample(.self$other_pdf %>% filter(Var1 == x_1 & Var2 == x_2)$Var3,
                        1,
                        prob = .self$other_pdf %>% filter(Var1 == x_1 & Var2 == x_2)$Freq)
            results[[h]] <- c(results[[h]], y)
          }
        }
      }
      return(results)
    },
    
    generate_households = function(num_househ, seed = NULL) {
      enc_households <- .self$generate_encoded_households(num_househ, seed)
      output <- .self$decode(do.call(c, enc_households))
      
      house_index <- unlist(lapply(seq_along(enc_households), function(i) rep(i, length(enc_households[[i]]))))
      output$house_index <- as.integer(house_index)
      
      return(output)
    },
    
    compute_marginals = function(dataframe, 
                                 hh_index = 'house_index', 
                                 sort_by = 'age', 
                                 ascending = FALSE) {
      
      househ_df <- dataframe %>% arrange(.dots = lapply(sort_by, as.name), .by_group = ascending) %>%
        group_by(.dots = lapply(hh_index, as.name))
      
      househ_df$features <- .self$encode(househ_df)
      
      counts <- list()
      model_features <- seq_len(.self$size())
      
      joint_tuples <- c(combn(model_features, 2, simplify = FALSE), lapply(model_features, function(x) c(x, x)))
      joint_index <- expand.grid(joint_tuples)
      
      counts$joint <- matrix(0, nrow = nrow(joint_index))
      
      other_tuples <- expand.grid(t(joint_tuples), model_features)
      other_index <- expand.grid(other_tuples)
      
      counts$other <- matrix(0, nrow = nrow(other_index))
      
      househ_size_tuples <- expand.grid(t(joint_tuples), unique(househ_df$hh_size))
      hh_size_index <- expand.grid(househ_size_tuples)
      
      counts$hh_size <- matrix(0, nrow = nrow(hh_size_index))
      
      for (id_hh in unique(househ_df[[hh_index]])) {
        feats <- househ_df %>% filter_(.dots = lapply(hh_index, function(x) househ_df[[x]] == id_hh)) %>%
          arrange(.dots = sort_by, ascending = ascending) %>%
          .$features
        
        x_1 <- feats[1]
        x_2 <- feats[2]
        
        counts$joint[which(joint_index$Var1 == x_1 & joint_index$Var2 == x_2)] <- counts$joint[which(joint_index$Var1 == x_1 & joint_index$Var2 == x_2)] + 1
        
        for (i in seq_len(length(feats) - 2) + 2) {
          counts$other[which(other_index$Var1 == x_1 & other_index$Var2 == x_2 & other_index$Var3 == feats[i])] <- counts$other[which(other_index$Var1 == x_1 & other_index$Var2 == x_2 & other_index$Var3 == feats[i])] + 1
        }
        
        size <- unique(househ_df %>% filter_(.dots = lapply(hh_index, function(x) househ_df[[x]] == id_hh)) %>%
                         .$hh_size)
        counts$hh_size[which(hh_size_index$Var1 == x_1 & hh_size_index$Var2 == x_2 & hh_size_index$Var3 == size)] <- counts$hh_size[which(hh_size_index$Var1 == x_1 & hh_size_index$Var2 == x_2 & hh_size_index$Var3 == size)] + 1
      }
      
      freqs <- list()
      
      freqs$joint <- counts$joint / sum(counts$joint)
      
      for (what in c('other', 'hh_size')) {
        freqs[[what]] <- matrix(0, nrow = nrow(counts[[what]]))
        
        for (i in unique(other_index$Var1)) {
          for (j in unique(other_index$Var2[other_index$Var2 >= i])) {
            norm <- sum(counts[[what]][other_index$Var1 == i & other_index$Var2 == j, ])
            if (norm > 0) {
              for (x_k in unique(other_index$Var3)) {
                freqs[[what]][which(other_index$Var1 == i & other_index$Var2 == j & other_index$Var3 == x_k)] <- counts[[what]][which(other_index$Var1 == i & other_index$Var2 == j & other_index$Var3 == x_k)] / norm
              }
            }
          }
        }
      }
      
      return(freqs)
    }
    
  )
)



#----- DEBUGGIN SECTION -----#
# Load database
database <- read.csv("../data/database.csv")

# Age bins
age_max <- ceiling(max(database$age))
age_bins <- c(0, 12, 25, 45, 65, age_max)

# Sex dictionary
sex_dict <- c(M = 0, F = 1)

encoding_dict <- list(
  age = age_bins,
  sex = sex_dict
)

# Create model instance
model <- Nucles_Model$new(encoding_dict)

# Compute marginals
marginals <- model$compute_marginals(database)

# Assign marginals
model$assign_marginals(marginals)

# Setting the number of households and seed
num_households <- 100  # Example number, adjust as needed
seed_list <- c(123, 456)  # Example seed list, adjust as needed

# Generate households
simulated <- model$generate_households(num_households, seed = seed_list[1])

# Show the simulated households
print(simulated)