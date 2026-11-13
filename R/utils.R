## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Unprintable separator hashing several values into one key.
key_sep_ <- "\x1F"

#' Hash the `cols` of each row of `x` into one key.
#' @noRd
row_keys_ <- function(x, cols) {
  do.call(paste, c(lapply(x[cols], as.character), list(sep = key_sep_)))
}

#' Thin draws to `n` evenly spaced elements, or rows of a `data.frame`,
#' first and last kept.
#' Spacing them out keeps autocorrelated draws as far apart as they can be.
#' * `n` NA or at least `NROW(x)` -> `x` unchanged.
#' @param x Vector or `data.frame` of draws.
#' @param n `integer(1)` Number of draws kept.
#' @noRd
helper_thin_even_ <- function(x, n) {
  if (is.na(n) || NROW(x) <= n) return(x)
  keep <- round(seq.int(1L, NROW(x), length.out = n))
  if (is.data.frame(x)) x[keep, , drop = FALSE] else x[keep]
}

#' Hash each family of a definition payload, row names of tables reset,
#' so the stored form of a table never changes its hash.
#' @param payload `list` Tables the engine imports, see `payload_()`.
#' @noRd
hash_payload_ <- function(payload) {
  canonical <- function(x) {
    if (is.data.frame(x)) `rownames<-`(x, NULL)
    else if (is.list(x)) lapply(x, canonical)
    else x
  }
  lapply(payload, \(x) rlang::hash(canonical(x)))
}

#' Whether `x` holds named parameter values, a vector or a list,
#' rather than posterior draws or a table.
#' @noRd
is_named_values_ <- function(x) {
  !posterior::is_draws(x) && is.null(dim(x)) && (is.atomic(x) || is.list(x))
}

# file.copy() with an overwrite option.
overwrite <- function(from, to, ...) {
  if (from == to) return(invisible())
  # user permission error, typically
  if (!suppressWarnings(file.copy(from, to, overwrite = TRUE, ...))) {
    raise_permissions_error()
  }
}

# File shipped by a model for each templated definition family.
template_def_files <- c(parameter = "param.R",
                        latent = "latent.R",
                        `infectivity profile` = "infecprof.R",
                        `status process` = "status.R")

#' Names of the entries of named list `a` that differ from those of `b`.
#' @noRd
differing_entries_ <- function(a, b) {
  names(a)[!mapply(identical, a, b[names(a)])]
}

#' Combine R6 `objects` into a clone of the first one.
#' Functions compare by arguments, body and enclosing environment.
#' * classes differ -> error;
#' * private field outside `pooled` differs -> error naming it,
#'   down to the entry of a named list;
#' * `secs_` pooled -> run times add up, other pooled fields left to the caller.
#' @param pooled `character()` Private fields the combination pools.
#' @param what `character(1)` Plural noun of the objects, for the error.
#' @return `list` The clone `out`, its private environment `pvt`,
#' and the private environment of each object `privs`.
#' @noRd
combine_r6_ <- function(objects, pooled, what) {
  classes <- unique(vapply(objects, \(x) class(x)[[1L]], character(1L)))
  if (length(classes) > 1L) stop_h(
    "Objects of class{?es} {.cls {classes}} cannot be combined.",
    "Combine {what} of one kind at a time."
  )
  privs <- lapply(objects, \(x) x$.__enclos_env__$private)
  first <- privs[[1L]]
  fields <- setdiff(names(Filter(Negate(is.function), as.list(first))), pooled)
  # name of the field, or of its entries, `x` and `y` differ in
  differ_in <- function(field, x, y) {
    a <- x[[field]]
    b <- y[[field]]
    entries <- is.list(a) && !is.data.frame(a) && !is.null(names(a)) &&
      identical(names(a), names(b))
    if (entries) return(differing_entries_(a, b))
    if (identical(a, b)) return()
    if (field == "name_") "model" else sub("_$", "", field)
  }
  differ <- unique(unlist(lapply(privs[-1L], function(other) {
    lapply(fields, differ_in, x = first, y = other)
  })))
  if (length(differ)) stop_h(
    "The {what} differ in {.field {differ}}.",
    "Only {what} of the same model, data and settings combine."
  )
  out <- objects[[1L]]$clone()
  pvt <- out$.__enclos_env__$private
  if ("secs_" %in% pooled)
    pvt$secs_ <- sum(vapply(privs, `[[`, numeric(1L), "secs_"))
  list(out = out, pvt = pvt, privs = privs)
}
