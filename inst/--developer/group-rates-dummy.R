
format_string <- function(name) {
  return(gsub("[[:punct:]]", "", name))
}

# options 
TAB <- "    "
TAB2 <- strrep(TAB, 2)

# @description 
# return logic table of group rates
# @param dataset
# the dataset from which the etract group rates
# @param features
# 
group_rate_logic_table <- function(dataset, features) {
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
  logic_table$..counts <- multifeat_count[multifeat_count>0]
  
  return(logic_table)
}

# @ description
# 
#
#
associate_group_names <- function(logic_table, prefix = NA, reference = NULL) {
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

##########################################
### construct the C++ switch statement ###
##########################################

# @ description
# function to produce if clauses
# @param r
# row of the group rate row
# @param start
# if `TRUE` the if clause is the first
cpp_if_clause <- function(r, is_ref = FALSE) {
  mask <- setdiff(names(r), "..associate")
  # if condition
  condition <- paste(paste0(
    sprintf('house.data(i, "%s") == ', names(r[mask])), 
    r[mask]), collapse = ") && (")
  # wrap all together
  c(sprintf('%s%sif ((%s))', TAB, ifelse(is_ref, "", "else "), condition),
    sprintf("%semplace_group_rate(r, %s);", 
            TAB2, 
            ifelse(is_ref, "REFERENCE", sprintf('"%s"', r[["..associate"]])))
  )
}


# @description 
# final wrapper
specify_group_rate <- function(dataset, 
                               features, 
                               rateID, 
                               prefix = NA, 
                               reference = NULL,
                               independent = TRUE) {
  # default prefix 
  if (is.na(prefix)) {
    prefix <-  format_string(tolower(rateID))
    prefix <- substring(prefix, 1, 3)
  }
  
  features <- as.list(features)
  if (isFALSE(independent)) {
    features <- list(features)
  }
  
  # initialize script
  script <- c()
  
  for (this_feats in features) {
    logic_table <- group_rate_logic_table(dataset, this_feats)
    logic_table <- associate_group_names(logic_table, prefix, reference)
    script <- c(
      script,
      # reference case
      cpp_if_clause(logic_table[logic_table$..associate == "REFERENCE",], TRUE),
      # other cases
      apply(logic_table[logic_table$..associate != "REFERENCE",], 1, cpp_if_clause),
      # close with missing case
      sprintf("%selse", TAB),
      sprintf('%semplace_group_rate(r, MISSING_RATE);', TAB2)
    )
  }
  
  return(c("{", sprintf('%srateID r = "%s";', TAB, rateID), script, "}"))
}

# @description
# 
assign_group_rate_to_individual <- function(dataset, definitions) {
  complete_cpp_script <- c()
  
  for (i in seq_along(definitions)){
    this_def <- definitions[[i]]
    this_def[["rateID"]] <- names(definitions[i])
    this_def[["dataset"]] <- dataset # UPDATE : this is not great...
    complete_cpp_script <- c(complete_cpp_script,
                             do.call(specify_group_rate, this_def))
  }
  cat(
    paste(c(
      "/* Assign all group rates for each individual.",
      "   Generated automatically -> do not edit by hand */",
      "void Model::assign_group_rate_to_individual(Household &house, indiv i)",
      "{",
      sapply(complete_cpp_script, function(x) { paste0(TAB, x) }),
      "}"
    ), 
    collapse='\n')
  )
}

##########################################
### construct the C++ switch statement ###
##########################################

library(MASS)

definitions <- list(
  Infectivity = list(features = list("ui", "ht")),
  Susc_Community = list(features = list("race", "smoke"), prefix = "com",
                        reference = list(smoke = 1, race = 2))
)

assign_group_rate_to_individual(birthwt, definitions)