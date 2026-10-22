library(testthat)

####  Helpers  ####

trapezoid <- \(x, y) 0.5 * sum(diff(x) * (y[-length(x)] + y[-1]))

# Rayleigh reference, named after the stats d/p/q convention.
drayleigh <- function(x, scale) {
  ifelse(x >= 0, x / scale^2 * exp(-x^2 / (2 * scale^2)), 0)
}
prayleigh <- function(q, scale, lower.tail = TRUE) {
  z <- pmax(q, 0)^2 / (2 * scale^2)
  if (lower.tail) -expm1(-z) else exp(-z)
}
qrayleigh <- \(p, scale) scale * sqrt(-2 * log1p(-p))

# R reference `<prefix><ref>` (e.g. "d" + "norm") bound to named arguments.
stats_ref <- function(prefix, ref, args) {
  f <- get(paste0(prefix, ref), mode = "function")
  \(q, ...) do.call(f, c(list(q), args, list(...)))
}

# Interval notation into the list the kernel reports.
interval <- function(txt) {
  ends <- as.numeric(strsplit(substr(txt, 2, nchar(txt) - 1), ", ")[[1]])
  list(min = ends[1], max = ends[2],
       closed = c(min = startsWith(txt, "["), max = endsWith(txt, "]")))
}

# Interval notation into a membership predicate.
in_interval <- function(txt) {
  bounds <- interval(txt)
  lo <- if (bounds$closed[["min"]]) `>=` else `>`
  up <- if (bounds$closed[["max"]]) `<=` else `<`
  \(x) lo(x, bounds$min) & up(x, bounds$max)
}

# Assert exposed function methods against R reference.
# A NULL reference marks a surface the kernel leaves undefined.
expect_function_matches <- function(name, args,
                                    ref_f, ref_prim = NULL,
                                    ref_log = log,
                                    domain = "(-Inf, Inf)", discrete = FALSE,
                                    x = seq(0.1, 10, 0.1)) {
  func <- do.call(browse_functions, c(list(name), args))
  func$domain |> expect_identical(interval(domain))
  func$is_discrete |> expect_identical(discrete)

  func$func(x) |> expect_equal(ref_f(x))

  if (is.null(ref_prim)) func$prim(x) |> expect_error("undefined method")
  else func$prim(x) |> expect_equal(ref_prim(x))

  if (is.null(ref_log)) func$log_func(x) |> expect_error("undefined method")
  else func$log_func(x) |> expect_equal(ref_log(func$func(x)))
}

# Assert exposed distribution methods against R reference `ref`.
# The reported support gates the values the null move is defined on.
expect_distribution_matches <- function(name, args, ref,
                                        domain = "(-Inf, Inf)",
                                        discrete = FALSE,
                                        x = seq(-20, 20, 1), h = pi) {
  distr <- do.call(browse_distributions, c(list(name), args))
  distr$domain |> expect_identical(interval(domain))
  distr$is_discrete |> expect_identical(discrete)

  distr$pdf(x) |> expect_equal(stats_ref("d", ref, args)(x))
  distr$cdf(x) |> expect_equal(stats_ref("p", ref, args)(x))
  distr$log_pdf(x) |> expect_equal(log(distr$pdf(x)))
  distr$log_ratio(x, x + h) |>
    expect_equal(log(distr$pdf(x + h) / distr$pdf(x)))
  distr$log_ratio(x + h, x) |>
    expect_equal(log(distr$pdf(x) / distr$pdf(x + h)))
  keep <- in_interval(domain)(x)
  null_move <- distr$log_ratio(x, x)
  null_move[keep] |> expect_equal(rep(0, sum(keep)))
  null_move[!keep] |> expect_all_na()
}

# Assert a truncated continuous law against its renormalized R reference.
expect_truncated_matches <- function(name, args, ref, x,
                                     min = -Inf, max = Inf) {
  bounds <- Filter(is.finite, list(min = min, max = max))
  distr <- do.call(browse_distributions, c(list(name), args, bounds))
  d <- stats_ref("d", ref, args)
  p <- stats_ref("p", ref, args)
  zz <- p(max) - p(min)
  distr$pdf(x) |> expect_equal(ifelse(x >= min & x <= max, d(x), 0) / zz)
  distr$cdf(x) |> expect_equal(pmin(pmax((p(x) - p(min)) / zz, 0), 1))
}

# Assert a discretized law stays finite far in the upper tail, where the naive
# difference of cumulatives underflows. `bin` holds the edges of the bin at x.
# Tail masses sit below the tolerance of `expect_equal`, which then compares
# them in absolute terms: only the log-scale comparison has resolution there.
expect_tail_finite <- function(name, args, ref, x, bin = c(0, 1)) {
  distr <- do.call(browse_distributions, c(list(name), args, discr = TRUE))
  p <- stats_ref("p", ref, args)
  any(p(x + bin[2]) == p(x + bin[1])) |> expect_true()
  pmf <- p(x + bin[1], lower.tail = FALSE) - p(x + bin[2], lower.tail = FALSE)
  distr$log_pdf(x) |> expect_all_finite()
  distr$pdf(x) |> expect_equal(pmf)
  distr$log_pdf(x) |> expect_equal(log(pmf))
}

# Assert a discretized law keeps its mass far in the lower tail, where the
# masses sit below the machine epsilon and `1 - exp(-z)` rounds them to 0.
expect_head_exact <- function(name, args, ref, x) {
  distr <- do.call(browse_distributions, c(list(name), args, discr = TRUE))
  p <- stats_ref("p", ref, args)
  pmf <- p(x + 1) - p(x)
  (pmf < .Machine$double.eps / 2) |> expect_all_true()
  distr$log_pdf(x) |> expect_equal(log(pmf))
}

# Rd source of a template help page:
# man/ under load_all, help db when installed
template_rd <- function(name) {
  topic <- paste0("template_", name, ".Rd")
  src <- system.file("man", topic, package = "housepi")
  if (nzchar(src)) return(paste(readLines(src), collapse = "\n"))
  paste(as.character(tools::Rd_db("housepi")[[topic]]), collapse = "")
}

test_that("Available names listing is stable", {
  # "uniform" is carried by both distribution registries and must be listed once
  browse_functions() |> expect_snapshot()
  browse_distributions() |> expect_snapshot()
  browse_infectivity_profiles() |> expect_snapshot()
})

test_that("Browsed functions do not leak memory", {
  mem_used_mb <- \() gc()["Vcells", "(Mb)"]

  func <- browse_functions("sin")
  rm(func)
  invisible(gc())
  mem_before <- mem_used_mb()

  # allocation loop
  for (i in seq_len(50)) {
    func <- browse_functions("sin")
    rm(func)
    if (i %% 10 == 0) gc()
  }
  invisible(gc())
  (mem_used_mb() - mem_before) |> expect_lt(5)
})

####  Method surface  ####

test_that("Exposed fields are read only", {
  distr <- browse_distributions("normal", mean = 0, sd = 1)
  (distr$domain <- "[0, 1]") |> expect_error("read only")
  (distr$args <- c(mean = 1)) |> expect_error("read only")
  (distr$is_discrete <- TRUE) |> expect_error("read only")
  (distr$is_positive <- TRUE) |> expect_error("read only")
  (distr$codomain <- "[0, 1]") |> expect_error("read only")
})

test_that("Exposed methods carry R formals", {
  distr <- browse_distributions("lognormal", meanlog = 0, sdlog = 1)
  distr$pdf() |>
    expect_error("[Aa]rgument.*x.*missing")
  distr$log_ratio(1) |>
    expect_error("[Aa]rgument.*xn.*missing")
  distr$log_ratio(1:2, 1) |>
    expect_error("[Ss]ize")
  distr$unknown() |>
    expect_error("attempt to apply non-function")
  # each family answers for its own methods only
  browse_functions("sin")$inst() |>
    expect_error("attempt to apply non-function")
  browse_functions("sin")$log_ratio(1, 2) |>
    expect_error("attempt to apply non-function")
})

test_that("Truncation bounds the reported domain, location stays apart", {
  # bounds read before the location shift, as the arguments do
  distr <- browse_distributions("gamma", shape = 2, scale = 1, loc = 3)
  distr$domain |> expect_equal(interval("[0, Inf)"))
  distr_2 <- browse_distributions("normal", mean = 0, sd = 1,
                                  min = -1, max = 2, loc = 10)
  distr_2$domain |> expect_equal(interval("[-1, 2]"))
  # arguments report the parametrization, location stays unreported
  distr <- browse_distributions("poisson", lambda = 3, loc = 1)
  distr$args |> expect_equal(c(lambda = 3))

  # an untruncated law reports the interval its own definition carries
  browse_distributions("poisson", lambda = 5)$domain |>
    expect_equal(interval("[0, Inf)"))
  # the upper bound tracks the 'size' argument it is aliased to
  browse_distributions("binomial", size = 10, prob = 0.4)$domain |>
    expect_equal(interval("[0, 10]"))

  # a truncated side closes, a side left at the edge of the law stays open
  browse_distributions("beta", shape1 = 2, shape2 = 2, min = 0.2)$domain |>
    expect_equal(interval("[0.2, 1)"))
  browse_functions("log", max = 5)$domain |> expect_equal(interval("(0, 5]"))
})

test_that("The codomain reports what a definition returns", {
  # a function is bounded on the values it accepts and on the values it gives
  # back, and the two are not the same interval
  func <- browse_functions("log")
  func$domain |> expect_equal(interval("(0, Inf)"))
  func$codomain |> expect_equal(interval("(-Inf, Inf)"))
  func_2 <- browse_functions("exp")
  func_2$domain |> expect_equal(interval("(-Inf, Inf)"))
  func_2$codomain |> expect_equal(interval("(0, Inf)"))
  # bounded whatever it reads
  browse_functions("sin")$codomain |> expect_equal(interval("[-1, 1]"))
  browse_functions("cos")$codomain |> expect_equal(interval("[-1, 1]"))
  # a bound "x" maps to the single value it returns
  browse_functions("minus", x = 3)$codomain |>
    expect_equal(interval("[-3, -3]"))
  browse_functions("exp", x = 0)$codomain |> expect_equal(interval("[1, 1]"))
  browse_functions("log", x = 1)$codomain |> expect_equal(interval("[0, 0]"))
  # an integer raised to a non-negative integer power stays an integer
  browse_functions("pow", x = 2, x2 = 3)$is_discrete |> expect_true()
  browse_functions("pow", x = 2, x2 = 0.5)$is_discrete |> expect_false()

  # a reduction folds its declared operands, skipping the unbound "x"
  browse_functions("sum", x2 = 5)$codomain |> expect_equal(interval("[5, 5]"))
  browse_functions("min", x2 = 5, x3 = -5)$codomain |>
    expect_equal(interval("[-5, -5]"))
  browse_functions("max", x2 = 5, x3 = -5)$codomain |>
    expect_equal(interval("[5, 5]"))

  # a distribution evaluates to a density, never to a draw: it answers here for
  # the mass it gives back, and in `domain` for the values it is defined on
  browse_distributions("gamma", shape = 2, scale = 1)$codomain |>
    expect_equal(interval("[0, Inf)"))
  # the location moves the values read in, never the density read out
  browse_distributions("gamma", shape = 2, scale = 1, loc = 3)$codomain |>
    expect_equal(interval("[0, Inf)"))
  # a mass is bounded by one, a density is not
  browse_distributions("poisson", lambda = 5)$codomain |>
    expect_equal(interval("[0, 1]"))
  browse_distributions("binomial", size = 10, prob = 0.4)$codomain |>
    expect_equal(interval("[0, 1]"))
  browse_distributions("gamma", shape = 2, scale = 1, discr = TRUE)$codomain |>
    expect_equal(interval("[0, 1]"))
})

test_that("Strict positivity follows the support", {
  # a continuous support closing at zero puts no mass on the boundary
  browse_distributions("gamma", shape = 2, scale = 1)$is_positive |>
    expect_true()
  browse_distributions("lognormal", meanlog = 0, sdlog = 1)$is_positive |>
    expect_true()
  browse_distributions("exponential", rate = 1)$is_positive |> expect_true()

  # a discrete law owning the zero bin does not, discretization included
  browse_distributions("poisson", lambda = 5)$is_positive |> expect_false()
  browse_distributions("binomial", size = 10, prob = 0.4)$is_positive |>
    expect_false()
  browse_distributions("gamma", shape = 2, scale = 1,
                       discr = TRUE)$is_positive |>
    expect_false()
  browse_distributions("uniform", min = 2, max = 6, discr = TRUE)$is_positive |>
    expect_true()

  # the floor is read in the absolute frame, after truncation and location
  browse_distributions("normal", mean = 0, sd = 1)$is_positive |> expect_false()
  browse_distributions("normal", mean = 0, sd = 1, min = 0)$is_positive |>
    expect_true()
  distr <- browse_distributions("gamma", shape = 2, scale = 1, loc = -1)
  distr$is_positive |> expect_false()

  # a function answers for the value it returns: the exponential is positive
  # whatever it reads, the others follow an "x" this surface leaves unbound
  browse_functions("exp")$is_positive |> expect_true()
  browse_functions("log")$is_positive |> expect_false()
  browse_functions("sqrt")$is_positive |> expect_false()
  browse_functions("div", x2 = 5)$is_positive |> expect_false()

  # a reduction folds over its declared operands, skipping the unbound "x"
  browse_functions("sum", x2 = 5)$is_positive |> expect_true()
  browse_functions("sum", x2 = -5)$is_positive |> expect_false()
  browse_functions("prod", x2 = 5)$is_positive |> expect_true()
  browse_functions("prod", x2 = -5)$is_positive |> expect_false()
  # the minimum needs every operand positive, the maximum just one
  browse_functions("min", x2 = 5, x3 = -5)$is_positive |> expect_false()
  browse_functions("max", x2 = 5, x3 = -5)$is_positive |> expect_true()
  browse_functions("min", x2 = 5, x3 = 1)$is_positive |> expect_true()
  browse_functions("max", x2 = -5, x3 = -1)$is_positive |> expect_false()
})

####  Functions  ####

test_that("Function arguments are validated", {
  x <- seq(0.1, 10, 0.1)
  func <- browse_functions("pow", x2 = 3)
  func$args |> expect_equal(c(x2 = 3))

  browse_functions("pow")$func(x) |>
    expect_error("Missing non default")
  browse_functions("pow", 5) |>
    expect_error("Unnamed argument")
  browse_functions("pow", x2 = 5, y2 = 3) |>
    expect_error("Unused argument.*y2.*[Dd]id you mean")
  browse_functions("pow", discr = TRUE)$func(x) |>
    expect_error("undefined method")
  # a minimum below the domain of the function
  browse_functions("sqrt", min = -1) |>
    expect_error("[Ii]nvalid min value -1.*sqrt")
  browse_functions("") |> expect_error("[Uu]nknown function")
  # an argument given twice
  browse_functions("sum", x2 = 1, x2 = 2) |>
    expect_error("[Rr]epeated argument .x2. for function .sum.")
  browse_distributions("normal", mean = 0, mean = 5) |>
    expect_error("[Rr]epeated argument .mean. for distribution .normal.")
})

test_that("Location shifts the function argument", {
  x <- seq(-10, 10, 0.1)
  loc_val <- 1.5
  browse_functions("exp", loc = loc_val)$func(x) |>
    expect_equal(exp(x - loc_val))
  browse_functions("exp", loc = loc_val)$prim(x) |>
    expect_equal(exp(x - loc_val))
  # a quarter-turn shift turns the sine into the cosine
  browse_functions("sin", loc = -0.5 * pi)$func(x) |> expect_equal(cos(x))
  # a reduction shifts its "x" operand alone
  browse_functions("sum", x2 = 1, loc = 2)$func(x) |> expect_equal(x - 1)
})

test_that("Transcendental functions match R reference", {
  expect_function_matches("pow", list(x2 = 3), \(x) x^3, \(x) x^4 / 4)
  # x2 = -1 primitive uses log(|x|)
  expect_function_matches("pow", list(x2 = -1), \(x) 1 / x, \(x) log(abs(x)))
  # x2 = 0.5 equals sqrt
  expect_function_matches("pow", list(x2 = 0.5), sqrt, \(x) x^1.5 / 1.5)
  expect_function_matches("sqrt", list(), sqrt, \(x) x^1.5 / 1.5,
                          domain = "[0, Inf)")
  expect_function_matches("exp", list(), exp, exp, x = seq(-5, 5, 0.1))
  expect_function_matches("log", list(), log, \(x) x * log(x) - x,
                          domain = "(0, Inf)", x = seq(1.1, 20, 0.1))
  expect_function_matches("sin", list(), sin, \(x) -cos(x),
                          x = seq(0.1, pi - 0.1, 0.05))
  expect_function_matches("cos", list(), cos, sin,
                          x = seq(-pi / 2 + 0.1, pi / 2 - 0.1, 0.05))
  # no primitive for the tangent
  expect_function_matches("tan", list(), tan, x = seq(0.1, pi / 2 - 0.1, 0.05))
})

test_that("Arithmetic kernels have no log surface", {
  arith <- function(name, args, ref, discrete = FALSE, x = seq(-5, 5, 0.5)) {
    expect_function_matches(name, args, ref, ref_log = NULL,
                            discrete = discrete, x = x)
  }
  fold <- c(3.3, 5.5, 6.6)
  # NOTE: avoid half-values: C++ std::round uses round-half-away-from-zero
  # while R's round() uses round-half-to-even (crazy...)
  x_step <- c(-3, -2.7, -1.3, -1, -0.3, 0, 0.4, 1.3, 2.7, 4)

  arith("minus", list(), \(x) -x)
  arith("sub", list(x2 = 3), \(x) x - 3)
  arith("sub", list(x2 = 0), identity)
  arith("div", list(x2 = 4), \(x) x / 4)
  arith("div", list(x2 = 1), identity)
  # integers are fixed points of the real-to-integer kernels
  arith("round", list(), round, discrete = TRUE, x = x_step)
  arith("ceil", list(), ceiling, discrete = TRUE, x = x_step)
  arith("floor", list(), floor, discrete = TRUE, x = x_step)
  # variadic kernels fold their constant arguments with x
  arith("sum", list(fold), \(x) x + sum(fold))
  arith("prod", list(fold), \(x) x * prod(fold))
  arith("min", list(fold), \(x) pmin(x, min(fold)))
  arith("max", list(fold), \(x) pmax(x, max(fold)))
})

test_that("Variadic kernels return NA on an NA operand", {
  # as in R without na.rm: sum(NA, 3) is NA
  for (name in c("sum", "prod", "min", "max")) {
    func <- browse_functions(name, x2 = 3)
    func$func(c(NA, 2)) |> expect_equal(c(NA, get(name)(2, 3)))
  }
})

####  Distributions  ####

test_that("Truncated discrete mass is normalised on the bins of the support", {
  # a lower bound drops the zero bin from the normalization
  x <- 1:10
  distr <- browse_distributions("binomial", size = 12, prob = 0.2, min = 1)
  shift <- stats::pbinom(0, 12, 0.2)
  distr$pdf(x) |> expect_equal(stats::dbinom(x, 12, 0.2) / (1 - shift))
  distr$cdf(x) |>
    expect_equal((stats::pbinom(x, 12, 0.2) - shift) / (1 - shift))

  # the window closes on the last bin inside the support, not one atom past it
  k <- 2:6
  distr <- browse_distributions("poisson", lambda = 3,
                                min = min(k), max = max(k))
  norm <- diff(stats::ppois(c(min(k) - 1, max(k)), 3))
  distr$pdf(k) |> expect_equal(stats::dpois(k, 3) / norm)
  sum(distr$pdf(k)) |> expect_equal(1)
  distr$cdf(max(k)) |> expect_equal(1)

  # left-discretized: bin i carries [i, i + 1), so the support holds 5 bins
  i <- 1:5
  distr <- browse_distributions("gamma", shape = 2, scale = 1,
                                min = min(i), max = max(i), discr = TRUE)
  norm <- diff(stats::pgamma(c(min(i), max(i) + 1), shape = 2, scale = 1))
  distr$pdf(i) |>
    expect_equal(diff(stats::pgamma(c(i, max(i) + 1), shape = 2,
                                    scale = 1)) / norm)
  sum(distr$pdf(i)) |> expect_equal(1)
  distr$cdf(max(i)) |> expect_equal(1)
})

test_that("Truncated discrete rand() reaches every bin of the support", {
  # both sampling branches must reproduce the mass they are scored against
  expect_draws_match_pdf <- function(distr, k, n = 2e4) {
    freq <- as.vector(table(factor(distr$rand(n), levels = k))) / n
    freq |> expect_equal(distr$pdf(k), tolerance = 0.05)
  }

  # retained mass 0.54: rejection branch
  k <- 2:6
  browse_distributions("poisson", lambda = 3, min = min(k), max = max(k)) |>
    expect_draws_match_pdf(k)

  # retained mass 0.17: inverse branch, whose window used to drop the first atom
  k <- 3:4
  browse_distributions("poisson", lambda = 3, min = min(k), max = max(k)) |>
    expect_draws_match_pdf(k)

  # left-discretized: the top bin is drawable, not only scorable
  i <- 1:5
  distr <- browse_distributions("gamma", shape = 2, scale = 1,
                                min = min(i), max = max(i), discr = TRUE)
  expect_draws_match_pdf(distr, i)
  max(distr$rand(2e4L)) |> expect_identical(max(i))
})

test_that("Discrete uniform spreads flat mass over the integers", {
  k <- 0:10
  distr <- browse_distributions("uniform", min = 0, max = 10, discr = TRUE)

  # every integer value in the support carries the same mass
  distr$pdf(k) |> expect_equal(rep(1 / length(k), length(k)))
  sum(distr$pdf(k)) |> expect_equal(1)
  distr$log_pdf(k) |> expect_equal(rep(-log(length(k)), length(k)))
  distr$pdf(c(-1, 11)) |> expect_equal(c(0, 0))

  # the cumulative closes on the last atom, and is the running mass
  distr$cdf(k) |> expect_equal(cumsum(distr$pdf(k)))
  distr$cdf(max(k)) |> expect_equal(1)

  # flat mass leaves every move even
  distr$log_ratio(c(0, 5), c(10, 6)) |> expect_equal(c(0, 0))

  # the quantile is the first atom whose cumulative reaches p
  p <- c(0.01, 0.25, 0.5, 0.75, 0.99)
  first_atom <- \(q) k[which(cumsum(distr$pdf(k)) >= q)[1]]
  distr$quantile(p) |> expect_equal(vapply(p, first_atom, 0))

  # the default support is a coin flip
  coin <- browse_distributions("uniform", discr = TRUE)
  coin$pdf(0:1) |> expect_equal(c(0.5, 0.5))
  coin$cdf(0:1) |> expect_equal(c(0.5, 1))

  # 'loc' carries the atoms without reshaping them
  browse_distributions("uniform", min = 0, max = 10, loc = 100,
                       discr = TRUE)$pdf(100 + k) |>
    expect_equal(distr$pdf(k))
})

test_that("Discretizing an already discrete law keeps its mass", {
  args <- list(size = 6, prob = 0.3)
  browse_distributions("binomial", args, discr = TRUE)$pdf(0:6) |>
    expect_equal(dbinom(0:6, 6, 0.3))

  weights <- c(a = 2, b = 1, c = 7)
  browse_distributions("categorical", weights, discr = TRUE)$pdf(0:2) |>
    expect_equal(unname(weights) / sum(weights))
})

test_that("An argument off the integers is refused where one is required", {
  # the flat lattice law owns every integer value in its support
  browse_distributions("uniform", min = 0.5, max = 3.2, discr = TRUE) |>
    expect_error("min.*uniform.*discrete")
  browse_distributions("uniform", min = 0, max = 3.2, discr = TRUE) |>
    expect_error("max.*uniform.*discrete")

  # a count of trials or successes is an integer whatever the support looks like
  browse_distributions("binomial", size = 10.5, prob = 0.5) |>
    expect_error("size.*binomial.*discrete")
  browse_distributions("negative_binomial", size = 5.5, prob = 0.5) |>
    expect_error("size.*negative_binomial.*discrete")
})

test_that("Options outside the support of a law are rejected", {
  browse_distributions("beta", shape1 = 2, shape2 = 3, discr = TRUE) |>
    expect_error("invalid.*discrete")
  browse_distributions("beta", shape1 = 2, shape2 = 3, min = -0.5) |>
    expect_error("min")
  browse_distributions("beta", shape1 = 2, shape2 = 3, max = 1.5) |>
    expect_error("max")
  browse_distributions("binomial", size = 10, prob = 0.5, min = -1) |>
    expect_error("min")
  browse_distributions("binomial", size = 10, prob = 0.5, min = 11) |>
    expect_error("max")
  # bounds crossing each other, reported in the located frame
  browse_distributions("normal", mean = 0, sd = 1, min = 2, max = 1, loc = 1) |>
    expect_error("[Ii]nconsistent domain \\[3, 2\\]")
})

test_that("Standard distributions match R reference d/p/log_pdf/log_ratio", {
  expect_distribution_matches("normal", list(mean = 2, sd = 3), "norm")
  expect_distribution_matches("uniform", list(min = -1, max = 5), "unif",
                              domain = "[-1, 5]")
  expect_distribution_matches("lognormal", list(meanlog = 1, sdlog = 0.5),
                              "lnorm", domain = "[0, Inf)")
  expect_distribution_matches("beta", list(shape1 = 2, shape2 = 1.5), "beta",
                              domain = "(0, 1)", x = seq(-10, 10, 0.1), h = 0.1)
  expect_distribution_matches("gamma", list(shape = 2, scale = 1.2), "gamma",
                              domain = "[0, Inf)")
  expect_distribution_matches("exponential", list(rate = 0.25), "exp",
                              domain = "[0, Inf)")
  expect_distribution_matches("poisson", list(lambda = 5), "pois",
                              domain = "[0, Inf)", discrete = TRUE)
  expect_distribution_matches("binomial", list(size = 10, prob = 0.5), "binom",
                              domain = "[0, 10]", discrete = TRUE)
  expect_distribution_matches("negative_binomial", list(size = 10, prob = 0.5),
                              "nbinom", domain = "[0, Inf)", discrete = TRUE)
  expect_distribution_matches("geometric", list(prob = 0.5), "geom",
                              domain = "[0, Inf)", discrete = TRUE)
  expect_distribution_matches("weibull", list(shape = 2, scale = 1), "weibull",
                              domain = "[0, Inf)")
  expect_distribution_matches("rayleigh", list(scale = 2), "rayleigh",
                              domain = "[0, Inf)", x = seq(0, 10, 0.5))
})

test_that("Categorical distribution is correct", {
  # unnormalized weights still give a normalized pmf
  args <- list(a = 2, b = 1, c = 7)
  distr <- browse_distributions("categorical", args)
  pmf <- unname(unlist(args)) / sum(unlist(args))
  distr$pdf(0:2) |> expect_equal(pmf)
  distr$pdf(c(-1, 3)) |> expect_equal(c(0, 0))
  distr$log_pdf(0:2) |> expect_equal(log(pmf))

  # cumulative pmf
  distr$cdf(0:2) |> expect_equal(cumsum(pmf))
  distr$cdf(c(-1, 5)) |> expect_equal(c(0, 1))

  # log-ratio
  distr$log_ratio(c(0, 1), c(1, 2)) |>
    expect_equal(log(pmf[c(2, 3)] / pmf[c(1, 2)]))

  # the quantile is the first label whose cumulative weight reaches p
  distr$quantile(c(0.1, 0.25, 0.5, 0.9)) |> expect_equal(c(0, 1, 2, 2))

  browse_distributions("categorical") |> expect_error("no weights")
  browse_distributions("categorical", a = 1, b = -2) |>
    expect_error("non-strictly positive")
})

test_that("Default parameters match R reference for d", {
  x <- seq(-20, 20, 1)
  browse_distributions("normal")$pdf(x) |> expect_equal(stats::dnorm(x))
  browse_distributions("uniform")$pdf(x) |> expect_equal(stats::dunif(x))
  browse_distributions("lognormal")$pdf(x) |> expect_equal(stats::dlnorm(x))
  browse_distributions("exponential")$pdf(x) |> expect_equal(stats::dexp(x))
})

test_that("Truncated continuous laws renormalize their R reference", {
  expect_truncated_matches("normal", list(mean = 2, sd = 3), "norm",
                           seq(-20, 20, 1), min = 0, max = 5)
  expect_truncated_matches("normal", list(mean = 2, sd = 1.5), "norm",
                           seq(-3, 3, 0.5), max = 4)
  expect_truncated_matches("lognormal", list(meanlog = 1, sdlog = 0.5), "lnorm",
                           seq(0, 8, 0.5), min = 1, max = 5)
  expect_truncated_matches("exponential", list(rate = 0.3), "exp",
                           seq(0, 8, 0.5), max = 5)
  expect_truncated_matches("beta", list(shape1 = 2, shape2 = 5), "beta",
                           seq(0, 1, 0.05), min = 0.1, max = 0.8)
})

test_that("Gamma and lognormal accept a mean/sd parametrization", {
  browse_distributions("gamma", mean = 2, scale = 1.1) |>
    expect_error("mix")
  browse_distributions("gamma", shape = 2, sd = 1.1) |>
    expect_error("mix")
  browse_distributions("lognormal", mean = 2, sdlog = 1.1) |>
    expect_error("mix")
  browse_distributions("lognormal", meanlog = 2, sd = 1.1) |>
    expect_error("mix")

  m <- 2.5
  s <- 1.3
  sdlog <- sqrt(log1p((s * s) / (m * m)))
  canonical <- list(
    gamma = list(shape = m * m / (s * s), scale = s * s / m),
    lognormal = list(meanlog = log(m) - 0.5 * sdlog * sdlog, sdlog = sdlog)
  )
  x <- seq(0.1, 12, 0.3)
  for (name in names(canonical)) {
    d_can <- do.call(browse_distributions, c(list(name), canonical[[name]]))
    d_alt <- browse_distributions(name, mean = m, sd = s)
    d_alt$pdf(x) |> expect_equal(d_can$pdf(x), info = name)
    d_alt$cdf(x) |> expect_equal(d_can$cdf(x), info = name)
    d_alt$log_pdf(x) |> expect_equal(d_can$log_pdf(x), info = name)
  }
})

test_that("Non-integer x maps to floor for Discrete::Defined distributions", {
  # pdf and cdf at k and any value in [k, k+1) are identical
  distr <- browse_distributions("poisson", lambda = 3)
  k <- seq(0, 9, 1)
  distr$pdf(k + 0.7) |> expect_equal(distr$pdf(k))
  distr$cdf(k + 0.7) |> expect_equal(distr$cdf(k))
})

test_that("Discretized exponential and gamma use left-discretization", {
  # left discretization: pdf(k) = cdf(k + 1) - cdf(k)
  x <- seq(0, 5, 1)
  distr <- browse_distributions("exponential", rate = 0.4, discr = TRUE)
  expected_pmf <- stats::pexp(x + 1, 0.4) - stats::pexp(x, 0.4)
  distr$pdf(x) |> expect_equal(expected_pmf)
  distr$cdf(x) |> expect_equal(cumsum(expected_pmf), tolerance = 1e-10)

  x <- seq(0, 10, 1)
  gcdf <- stats_ref("p", "gamma", list(shape = 2, scale = 1.5))
  browse_distributions("gamma", shape = 2, scale = 1.5, discr = TRUE)$pdf(x) |>
    expect_equal(gcdf(x + 1) - gcdf(x))
})

test_that("Discretized normal methods go through center-discretization", {
  args <- list(mean = 5, sd = 2)
  distr <- browse_distributions("normal", args, discr = TRUE)
  normcdf <- stats_ref("p", "norm", args)

  # center discretization: pdf(k) = pnorm(k+0.5) - pnorm(k-0.5)
  x <- seq(0, 10, 1)
  expected_pmf <- normcdf(x + 0.5) - normcdf(x - 0.5)
  distr$pdf(x) |> expect_equal(expected_pmf)

  # log_func_ Center branch
  distr$log_pdf(x) |> expect_equal(log(expected_pmf))
  # primitive_ Center branch (cdf is the primitive)
  distr$cdf(x) |> expect_equal(normcdf(x + 0.5))
  # log_ratio_ Center branch with distinct rounded values
  xn <- x + 1
  distr$log_ratio(x, xn) |>
    expect_equal(log(distr$pdf(xn) / distr$pdf(x)))
  # log_ratio when rounded values coincide should be 0
  distr$log_ratio(x, x + 0.1) |> expect_equal(rep(0, length(x)))

  # truncation at max renormalizes over the bins below it
  args <- list(mean = 2, sd = 1.5)
  normcdf <- stats_ref("p", "norm", args)
  distr <- browse_distributions("normal", args, max = 4, discr = TRUE)
  x <- seq(-3, 5, 1)
  expec_pmf <- ifelse(x > 4, 0, normcdf(x + 0.5) - normcdf(x - 0.5)) /
    normcdf(4.5)
  distr$pdf(x) |> expect_equal(expec_pmf)
})

test_that("Discretized lognormal keeps the mass of the origin bin", {
  args <- list(meanlog = 0, sdlog = 1)
  x <- 0:5000
  distr <- browse_distributions("lognormal", args, discr = TRUE)

  lncdf <- stats_ref("p", "lnorm", args)
  distr$pdf(x) |> expect_equal(lncdf(x + 1) - lncdf(x))
  sum(distr$pdf(x)) |> expect_equal(1)

  # the pmf and the cdf agree on the origin bin
  distr$pdf(0) |> expect_gt(0)
  distr$cdf(0:3) |> expect_equal(cumsum(distr$pdf(0:3)))

  # in continuous time the origin is a zero-density point, not a NaN
  contin <- browse_distributions("lognormal", args)
  contin$pdf(0) |> expect_equal(0)
  contin$log_pdf(0) |> expect_equal(-Inf)
  contin$log_ratio(0, 0) |> expect_equal(0)
  contin$log_ratio(0, 2) |> expect_equal(Inf)
  contin$log_ratio(2, 0) |> expect_equal(-Inf)
})

test_that("Discretized laws stay finite far in the upper tail", {
  # complementary-CDF branches: center bins for the normal, gamma/erfc
  # complement for the lognormal, boost gamma_q complement for the gamma,
  # closed forms for the weibull and rayleigh
  expect_tail_finite("normal", list(mean = 0, sd = 1), "norm", seq(9, 30, 1),
                     bin = c(-0.5, 0.5))
  expect_tail_finite("lognormal", list(meanlog = 0, sdlog = 0.5), "lnorm",
                     seq(60, 200, 5))
  expect_tail_finite("gamma", list(shape = 2, scale = 1.5), "gamma",
                     seq(55, 120, 1))
  expect_tail_finite("weibull", list(shape = 2, scale = 1), "weibull",
                     seq(7, 12, 1))
  expect_tail_finite("rayleigh", list(scale = 1), "rayleigh", seq(9, 30, 1))

  # symmetry: reflecting about the mean matches the lower tail
  distr <- browse_distributions("normal", mean = 0, sd = 1, discr = TRUE)
  distr$pdf(9:30) |> expect_equal(distr$pdf(-(9:30)))
})

test_that("Discretized laws keep their mass far in the lower tail", {
  expect_head_exact("exponential", list(rate = 1e-17), "exp", 0:5)
  expect_head_exact("weibull", list(shape = 0.5, scale = 1e40), "weibull", 0:5)
  expect_head_exact("rayleigh", list(scale = 1e9), "rayleigh", 0:5)
})

test_that("Geometric cumulative keeps its mass below the machine epsilon", {
  # `1 - (1 - prob)^(x + 1)` rounds to 0
  distr <- browse_distributions("geometric", prob = 1e-17)
  distr$cdf(0:5) |> expect_equal(pgeom(0:5, 1e-17))
})

test_that("Geometric quantile is left-continuous at the cumulative atoms", {
  # atoms up to `1 - cdf` near 1e-14, where `p` resolves single masses
  distr <- browse_distributions("geometric", prob = 0.01)
  distr$quantile(distr$cdf(0:3000)) |> expect_identical(as.double(0:3000))
  distr$quantile(pgeom(0:3000, 0.01)) |> expect_identical(as.double(0:3000))
  # a tiny probability, where the cumulative loses relative precision
  distr_2 <- browse_distributions("geometric", prob = 1e-6)
  distr_2$quantile(distr_2$cdf(0:2000)) |> expect_identical(as.double(0:2000))
})

test_that("Poisson cumulative keeps its mass below the machine epsilon", {
  # `1 - P(X > x)` rounds to 0
  distr <- browse_distributions("poisson", lambda = 100)
  distr$cdf(0:5) |> expect_equal(ppois(0:5, 100))
})

test_that("rand() draws are reproducible and statistically plausible", {
  d <- browse_distributions("normal", mean = 2, sd = 3)
  set.seed(42L)
  x1 <- d$rand(50L)
  set.seed(42L)
  x2 <- d$rand(50L)
  x1 |> expect_identical(x2) # same R seed -> same draws
  x1 |> expect_length(50L)
  # different seed -> different draws
  set.seed(43L)
  identical(x1, d$rand(50L)) |> expect_false()

  # mean / sd of 20k draws should be close to the truth
  big <- d$rand(20000L)
  mean(big) |> expect_equal(2, tolerance = 0.05)
  sd(big) |> expect_equal(3, tolerance = 0.05)

  # without set.seed, successive calls are non-deterministic
  identical(d$rand(10L), d$rand(10L)) |> expect_false()

  # invalid n
  d$rand(-1L) |>
    expect_error("non-negative")
})

test_that("rand() honours location shift and truncation", {
  d <- browse_distributions("normal", mean = 0, sd = 1, loc = 5)
  d$pdf(5) |> expect_equal(stats::dnorm(0))
  mean(d$rand(20000L)) |> expect_equal(5, tolerance = 0.05)

  x <- browse_distributions("normal", mean = 0, sd = 1,
                            min = -1, max = 2)$rand(5000L)
  (x >= -1 & x <= 2) |> expect_all_true()
})

test_that("cdf honours location shift above the upper bound", {
  args <- list(min = 2, max = 4, loc = 10)
  shifted <- browse_distributions("uniform", args)
  args <- as.list(unlist(args) + args$loc)
  args$loc <- NULL
  plain <- browse_distributions("uniform", args)
  x <- seq(9, 17, 0.5)
  shifted$pdf(x) |> expect_equal(plain$pdf(x))
  shifted$cdf(x) |> expect_equal(plain$cdf(x))

  # shifting by 'loc' does not touch arguments
  args <- list(shape = 2, scale = 1, loc = 3, max = 4)
  shifted <- browse_distributions("gamma", args)
  args$loc <- NULL
  plain <- browse_distributions("gamma", args)
  shifted$cdf(x) |> expect_equal(plain$cdf(x - 3))
  shifted$cdf(20) |> expect_equal(1)
})

test_that("Discrete draws from rand() have type integer", {
  browse_distributions("normal", mean = 0, sd = 1)$rand(5L) |>
    expect_type("double")
  browse_distributions("normal", mean = 0, sd = 1, discr = TRUE)$rand(5L) |>
    expect_type("integer")

  browse_distributions("poisson", lambda = 3)$rand(5L) |>
    expect_type("integer")
  # an integer location keeps the draws on the lattice
  browse_distributions("poisson", lambda = 3, loc = 2)$rand(5L) |>
    expect_type("integer")
  # past machine integer max, draw back to double
  browse_distributions("poisson", lambda = 3e9)$rand(2L) |>
    expect_type("double")
})

test_that("Location off the integers is refused when discrete", {
  browse_distributions("poisson", lambda = 3, loc = 0.5) |>
    expect_error("[Ss]et.*loc.*integer")
  browse_distributions("gamma", shape = 2, discr = TRUE, loc = -1.5) |>
    expect_error("[Ss]et.*loc.*integer")
  # an integer shift, and a fractional bound that only truncates, are fine
  browse_distributions("poisson", lambda = 3, loc = -2)$pdf(-2) |>
    expect_equal(stats::dpois(0, 3))
  browse_distributions("poisson", lambda = 3, min = 0.5)$domain |>
    expect_equal(interval("[0.5, Inf)"))
  # a continuous definition is free to sit anywhere
  browse_distributions("gamma", shape = 2, scale = 1, loc = 0.5)$pdf(1.5) |>
    expect_equal(stats::dgamma(1, shape = 2, scale = 1))
})

test_that("quantile returns NA outside [0, 1]", {
  distr <- browse_distributions("normal", mean = 0, sd = 1)
  distr$quantile(c(-2, -0.1, -0.001, 1.001, 1.5, 7.7)) |> expect_all_na()
})

test_that("quantile matches R reference for standard distributions", {
  # the grid closes on 0 and 1, where an unbounded support gives -Inf and Inf
  p <- seq(0, 1, 0.1)
  cases <- list(
    list("normal", list(mean = 2, sd = 3), "norm"),
    list("uniform", list(min = -3, max = 2), "unif"),
    list("lognormal", list(meanlog = 1, sdlog = 0.5), "lnorm"),
    list("gamma", list(shape = 2, scale = 1.5), "gamma"),
    list("gamma", list(mean = 2.2, sd = sqrt(1.1)), "gamma",
         list(shape = 4.4, scale = 0.5)),
    list("exponential", list(rate = 0.3), "exp"),
    list("beta", list(shape1 = 2, shape2 = 3), "beta"),
    list("weibull", list(shape = 2, scale = 2), "weibull"),
    list("rayleigh", list(scale = 2), "rayleigh"),
    list("poisson", list(lambda = 5), "pois"),
    list("binomial", list(size = 10, prob = 0.4), "binom"),
    list("geometric", list(prob = 0.3), "geom"),
    list("negative_binomial", list(size = 5, prob = 0.4), "nbinom")
  )
  for (case in cases) {
    distr <- do.call(browse_distributions, c(case[1], case[[2]]))
    ref <- stats_ref("q", case[[3]], case[[if (length(case) == 4) 4 else 2]])
    distr$quantile(p) |> expect_equal(ref(p), info = case[[1]])
  }
})

test_that("quantile is truncation-aware", {
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  cases <- list(
    list(args = list("normal", mean = 0, sd = 1, min = 0),
         ref = \(p) qnorm(pnorm(0) + p * (1 - pnorm(0)))),
    list(args = list("normal", mean = 3, sd = 2, min = 1, max = 5),
         ref = function(p) {
           qnorm(pnorm(1, 3, 2) + p * (pnorm(5, 3, 2) - pnorm(1, 3, 2)), 3, 2)
         }),
    list(args = list("exponential", rate = 1, max = 3),
         ref = \(p) qexp(p * pexp(3))),
    list(args = list("gamma", shape = 2, scale = 1.5, min = 1, max = 8),
         ref = function(p) {
           qgamma(pgamma(1, shape = 2, scale = 1.5) +
                    p * (pgamma(8, 2, scale = 1.5) - pgamma(1, 2, scale = 1.5)),
                  shape = 2, scale = 1.5)
         })
  )
  for (case in cases) {
    distr <- do.call(browse_distributions, case$args)
    distr$quantile(p) |>
      expect_equal(case$ref(p), tolerance = 1e-10, info = case$args[[1]])
    all(distr$quantile(p) >= max(-Inf, case$args$min)) |>
      expect_true(info = paste(case$args[[1]], "above min"))
    all(distr$quantile(p) <= min(Inf, case$args$max)) |>
      expect_true(info = paste(case$args[[1]], "below max"))
  }
})

test_that("Gamma pdf and log_pdf handle x == 0 across shape regimes", {
  scale <- 2
  # test case: shape < 1
  d1 <- browse_distributions("gamma", shape = 0.5, scale = scale)
  d1$pdf(0) |> expect_equal(Inf)
  d1$log_pdf(0) |> expect_equal(Inf)

  # test case: shape == 1
  d2 <- browse_distributions("gamma", shape = 1, scale = scale)
  d2$pdf(0) |> expect_equal(1 / scale)
  d2$log_pdf(0) |> expect_equal(-log(scale))

  # test case: shape > 1
  d3 <- browse_distributions("gamma", shape = 2, scale = scale)
  d3$pdf(0) |> expect_equal(0)
  d3$log_pdf(0) |> expect_equal(-Inf)
})

test_that("rand() is reproducible and plausible for all distributions", {
  cases <- list(
    list(args = list("rayleigh", scale = 2),
         mean = 2 * sqrt(pi / 2), int = FALSE),
    list(args = list("exponential", rate = 0.5),
         mean = 2, int = FALSE),
    list(args = list("beta", shape1 = 2, shape2 = 5),
         mean = 2 / 7, int = FALSE, max = 1),
    list(args = list("poisson", lambda = 5),
         mean = 5, int = TRUE),
    list(args = list("binomial", size = 20, prob = 0.3),
         mean = 6, int = TRUE, max = 20),
    list(args = list("negative_binomial", size = 5, prob = 0.4),
         mean = 5 * 0.6 / 0.4, int = TRUE),
    list(args = list("geometric", prob = 0.3),
         mean = 0.7 / 0.3, int = TRUE),
    list(args = list("categorical", a = 2, b = 1, c = 7),
         mean = 1.5, int = TRUE, max = 2)
  )
  for (cc in cases) {
    d <- do.call(browse_distributions, cc$args)
    nm <- cc$args[[1]]
    set.seed(42L)
    x1 <- d$rand(50L)
    set.seed(42L)
    d$rand(50L) |> expect_identical(x1, info = nm)
    big <- d$rand(2000L)
    all(big >= 0) |> expect_true(info = nm)
    if (cc$int) all(big == floor(big)) |> expect_true(info = nm)
    if (!is.null(cc$max)) all(big <= cc$max) |> expect_true(info = nm)
    mean(big) |> expect_equal(cc$mean, tolerance = 0.15, info = nm)
  }
})

####  Infectivity profiles  ####

test_that("Ferretti20 infectivity profile integrates and discretizes", {
  x <- seq(0, 15, 0.001)
  profile <- browse_infectivity_profiles("ferretti20", incub_period = 5)
  trapezoid(x, profile$inst(x)) |> expect_equal(profile$cumul(x)[length(x)])
  profile$is_discrete |> expect_false()
  profile$domain |> expect_equal(interval("[0, Inf)"))

  # left discretization works as expected
  x <- seq(1.0, 1.9, 0.1)
  profile <- browse_infectivity_profiles("ferretti20",
                                         incub_period = 5, discr = TRUE)
  profile$inst(x) |> expect_all_equal(profile$cumul(1))
  profile$is_discrete |> expect_true()
  profile$domain |> expect_equal(interval("[1, Inf)"))
})

test_that("Discretized Ferretti20 stays finite far in the upper tail", {
  # default shape, scale, c_loc and mean_incub, incubation period of 5 days
  a <- 5.85
  s <- 1.85
  m <- -4
  i_t <- 5 / 5.42
  pema <- (1 + exp(m / s))^-a
  norm <- 1 + (i_t - 1) * pema - i_t * (1 + exp((m + 5.42) / s))^-a
  # complementary cumulative after onset, 1 - (1 + exp((m - t) / s))^-a
  ccdf <- function(x) -expm1(-a * log1p(exp((m - x + 5) / s))) / norm

  distr <- browse_infectivity_profiles("ferretti20", incub_period = 5)
  (1 - distr$cumul(seq(5, 20, 0.5))) |> expect_equal(ccdf(seq(5, 20, 0.5)))

  # bins whose cumulatives round to the same double,
  # masses renormalized over the discretized domain [1, Inf)
  x <- seq(70, 120, 5)
  distr_2 <- browse_infectivity_profiles("ferretti20",
                                         incub_period = 5, discr = TRUE)
  (distr$cumul(x + 1) == distr$cumul(x)) |> expect_all_true()
  log(distr_2$inst(x)) |>
    expect_equal(log(ccdf(x) - ccdf(x + 1)) - log1p(-distr$cumul(1)))
})

####  Print method  ####

test_that("Exposed object print behavior is stable", {
  print(browse_functions("pow", x2 = 3)) |>
    expect_snapshot()
  print(browse_distributions("normal", mean = 2, sd = 3)) |>
    expect_snapshot()
  print(browse_infectivity_profiles("gamma", shape = 2, scale = 1)) |>
    expect_snapshot()
  # no-argument object prints the header and the support
  print(browse_functions("sin")) |>
    expect_snapshot()
  # truncating to the same bound leaves the values on the integers
  print(browse_distributions("poisson", lambda = 3, min = 0.5)) |>
    expect_snapshot()
  # an integer shift reports the argument and no grid
  print(browse_distributions("poisson", lambda = 3, loc = 2)) |>
    expect_snapshot()
  # a continuous definition reports the argument and no grid
  print(browse_distributions("gamma", shape = 2, scale = 1, loc = 3)) |>
    expect_snapshot()
})

####  Template models  ####

test_that("Template name normalization", {
  misstypes <- c("Final-Size", "final_size", " final siZe ", "FINAL SIZE")
  for (f in misstypes) spaced_choice(f) |> expect_identical("final size")

  find_closest_template_name("final siz") |>
    expect_match("[Dd]id you mean", all = FALSE)
  find_closest_template_name("") |>
    expect_match("browse_template", all = FALSE)
})

test_that("browse_templates branches inform correctly", {
  browse_templates(42) |> expect_error("not a character")
  browse_templates("unknown") |> expect_error("not a known template")
  browse_templates("final siz") |> expect_error("[Dd]id you mean")

  # print existing templates
  browse_templates() |> expect_snapshot()
})

test_that("browse_templates triggers help preview when path is NA", {
  # mock the internal helper so the test never spawns a real previewer
  called_with <- NULL
  with_mocked_bindings(
    open_model_helper_ = function(path) {
      called_with <<- path
      invisible()
    },
    code = browse_templates("final-size")
  )
  called_with |> expect_identical("final_size")
})

test_that("Every template has a help page", {
  for (name in list.files(template_models_dir())) {
    template_rd(name) |> expect_match("\\\\name\\{template_")
  }
})

test_that("browse_templates copies template to a new path", {
  dest <- tempfile("template_")
  on.exit(unlink(dest, recursive = TRUE), add = TRUE)
  model_cpp <- file.path(dest, "model.cpp")

  # first copy: creates the model sources
  browse_templates("final-size", path = dest)
  file.exists(model_cpp) |> expect_true()
  file.exists(file.path(dest, "param.R")) |> expect_true()
  # final_size has no latent file: nothing should be copied
  file.exists(file.path(dest, "latent.R")) |> expect_false()
  shipped <- readLines(model_cpp)

  # second copy without replace: warns and does not overwrite
  writeLines("edited", model_cpp)
  browse_templates("final-size", path = dest) |>
    expect_warning("already exists")
  readLines(model_cpp) |> expect_identical("edited")

  # with replace = TRUE: overwrites silently
  browse_templates("final-size", path = dest, replace = TRUE)
  readLines(model_cpp) |> expect_identical(shipped)
})

test_that("browse_templates rejects '.' as destination path", {
  sandbox <- tempfile("sandbox_")
  dir.create(sandbox)
  on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)
  withr::with_dir(sandbox, {
    browse_templates("final-size", path = ".") |>
      expect_warning("current directory")
    # helper redirects "." to a folder named after the template
    dir.exists("final-size") |> expect_true()
  })
})

test_that("Model help pages name the keys the engine registers", {
  modules <- c(Data = "Data", Param = "Parameters",
               Latent = "Latent", GroupRisk = "Group Risks")
  for (name in c("final_size", "case_ascertained")) {
    cpp <- system.file("models", name, "model.cpp", package = "housepi") |>
      readLines() |>
      paste(collapse = "\n")
    rd <- template_rd(name)
    # loop over all modules
    for (k in names(modules)) {
      # EMPTY_KEYS declarations match nothing and register no key
      fmt <- "(?<=REGISTER_KEYS\\(%s,)[^)]*"
      decl <- regmatches(cpp, regexpr(sprintf(fmt, k), cpp, perl = TRUE))
      keys <- trimws(unlist(strsplit(decl, ",")))
      # section body up to the next section: its \code{} items are the keys
      fmt <- "^.*\\\\section\\{%s\\}\\{(.*?)\\\\section.*$"
      section <- sub(sprintf(fmt, modules[[k]]), "\\1", rd)
      if (identical(section, rd)) section <- ""
      pattern <- "(?<=\\\\item \\\\code\\{)[^}]+"
      listed <- regmatches(section, gregexpr(pattern, section, perl = TRUE))
      listed[[1]] |> expect_setequal(keys)
    }
  }
})
