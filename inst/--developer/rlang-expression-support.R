expr_to_row <- function(expr) {

  # Accept either a formula (~ x - y) or a naked call (x - y)
  if (inherits(expr, "formula")) {
    call <- rlang::f_rhs(expr)
  } else if (rlang::is_call(expr)) {
    call <- expr
  } else {
    stop("expr must be a formula or a call")
  }

  # Name of the operator / function
  fn  <- rlang::call_name(call)          # e.g. "-", "+", "log", "sqrt"
  args <- rlang::call_args(call)
  n_arg <- length(args)

  # Map to label
  fn_label <- switch(fn,
    `-` = if (n_arg == 1) "unary_minus" else "subtraction",
    `+` = if (n_arg == 1) "unary_plus" else "addition",
    `*` = "multiplication",
    `/` = "division",
    `^` = "power",
    log  = "logarithm",
    sqrt = "square_root",
    exp  = "exponential",
    # default: use raw name
    fn                             
  )

  y <- c(func = fn_label,
         x1 = rlang::as_label(args[[1]]))
  if (n_arg == 2) y["x2"] <- rlang::as_label(args[[2]])
  # else
  y
}


expr_to_row(~ 1 + var_2)
expr_to_row(~ var_1@D - var_2@P)
expr_to_row(~ log(var_1@D))

# TO BE DEVELOPED
expr_to_row(1)