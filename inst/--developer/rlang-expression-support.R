
# if it's a call to "(" with a single argument that's a symbol, unwrap it
strip_meaningless_parens <- function(arg) {
  while (rlang::is_call(arg) && 
         rlang::call_name(arg) == "(" && 
         length(rlang::call_args(arg)) == 1) {
    inner_arg <- rlang::call_args(arg)[[1]]
    if (rlang::is_symbol(inner_arg)) arg <- inner_arg
    else break
  }
  return(arg)
}

expr_to_row <- function(expr) {
  # accept either a formula (~ x - y) or a naked call (x - y)
  if (inherits(expr, "formula")) call <- rlang::f_rhs(expr)
  else if (rlang::is_call(expr)) call <- expr
  else return(c(func = "eq", x = expr))
  
  # flatten nested operations of the same type
  flatten_same_op_ <- function(call, target_op) {
    if (rlang::is_call(call) && rlang::call_name(call) == target_op) {
      args <- rlang::call_args(call)
      output <- list()
      for (a in args) output <- c(output, flatten_same_op_(a, target_op))
      return(output)
    } else {
      return(list(call))
    }
  }
  
  # name of the operator / function
  fn <- rlang::call_name(call)          
  
  # handle nested operations for commutative operators
  if (fn %in% c("+", "*")) args <- flatten_same_op_(call, fn)
  else args <- rlang::call_args(call)
  
  # strip meaningless parentheses from arguments
  args <- lapply(args, strip_meaningless_parens)
  n_arg <- length(args)
  
  # map to function label
  fn_label <- switch(fn,
                     `-` = if (n_arg == 1) "minus" else "sub",
                     `+` = if (n_arg == 1) "eq" else "add",
                     `*` = "mult",
                     `/` = "div",
                     `^` = "pow",
                     log  = "log",
                     sqrt = "sqrt",
                     exp  = "exp",
                     sin = "sin",
                     cos = "cos",
                     tan = "tan",
                     fn) # default: use raw name
  
  output <- c(func = fn_label)
  
  # list operands
  for (i in seq_along(args)) {
    if (i == 1) output["x"] <- rlang::as_label(args[[i]])
    else output[paste0("x", i)] <- rlang::as_label(args[[i]])
  }
  
  return(output)
}

expr_to_row(-1)
expr_to_row(~+ a)
expr_to_row(~- a)
expr_to_row(~a + b)
expr_to_row(~a + 3 * b) # <- this kind of expression shall be developed
expr_to_row(~a + b + c)
expr_to_row(~a + b - c) # <- this kind of expression shall be developed
expr_to_row(~x * (y))
expr_to_row(~ (a) - (c))
expr_to_row(~log(x))

