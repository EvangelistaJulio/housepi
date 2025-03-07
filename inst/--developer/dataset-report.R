#' @title
#' Show datset report
#' @description
#' Counts typical quantities for the analysis of household surveys.
#'
#' @param df `data.frame`
#' A dataset.
#'
#' @return `flextable`
show_dataset_report <- function(df) {
  # initialize report
  report_ind <- list(total = NA, missing = NA, healthy = NA, infected = NA)
  report_hh <- list(total = NA, "no cases" = NA,
                    "at least 1 case" = NA, "at least 2 cases" = NA)
  report_size <- list(min_size = NA, max_size = NA)

  # count individuals
  report_ind$total <- ifelse("id" %in% colnames(df),
                             length(unique(df$id)), nrow(df))

  if ("status" %in% colnames(df)) {
    # format statuses to -1, 0, +1
    df$status <- as.integer(df$status)
    df[is.na(df$status), ]$status <- (-1)
    df[df$status >= 1, ]$status <- 1

    counts <- table(df$status)
    report_ind$missing <- ifelse("-1" %in% names(counts), counts["-1"], 0)
    report_ind$healthy <- ifelse("0" %in% names(counts), counts["0"], 0)
    report_ind$infected <- ifelse("1" %in% names(counts), counts["1"], 0)
  }

  if ("house_id" %in% colnames(df)) {
    # count households
    report_hh$total <- length(unique(df$house_id))

    if ("status" %in% colnames(df)) {
      # count household cases
      house_cases <- tapply(df$status, df$house_id, function(x) sum(x == 1))

      report_hh$`no cases` <- sum(house_cases == 0)
      report_hh$`at least 1 case` <- sum(house_cases >= 1)
      report_hh$`at least 2 cases` <- sum(house_cases >= 2)

      # check for household of missing statuses
      if (any(tapply(df$status, df$house_id, function(x) all(x == -1)))) {
        warning(
          "Found household with only missing status.",
          "Remove it before proceeding with analysis."
        )
      }
    }

    # count household members
    house_sizes <- table(df$house_id)

    # household size
    report_size$min_size = min(house_sizes)
    report_size$max_size = max(house_sizes)
  }

  wrap_report <- function(x, what) {
    out <- data.frame(counts = unlist(x), check.names = FALSE)
    out$what <- what
    out$comment <- row.names(out)
    out$pct <- 100 * unlist(x) / x$total
    return(out)
  }
  combined <- rbind(wrap_report(report_ind, "individual"),
                    wrap_report(report_hh, "household"))
  # re-order
  combined <- combined[, c("what", "comment", "counts", "pct")]


  combined$pct <- sapply(combined$pct, function(x) sprintf("%.0f%%", x))

  ft <- flextable::flextable(combined)
  ft <- flextable::separate_header(ft)
  ft <- flextable::theme_vanilla(ft)
  ft <- flextable::align(ft, align = "left", part = "all")
  ft <- flextable::merge_v(ft, j = ~ what)
  ft <- flextable::bold(ft, j = ~ what, bold = TRUE)
  ft <- flextable::autofit(ft, add_w = 0, add_h=0)
  ft <- flextable::labelizor(
    ft, part = "all",
    labels = c(what = "", comment = "", pct = "[%]")
  )

  # FIXME: missing information about size

  return(ft)
}
