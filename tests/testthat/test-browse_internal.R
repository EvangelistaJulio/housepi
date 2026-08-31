library(testthat)

################################################################################
###################################  HELPERS  ##################################
################################################################################

trapezoid <- function(x, y) {
  0.5 * sum(diff(x) * (y[-length(x)] + y[-1]))
}

rayleigh_d <- function(x, s) {
  y <- rep(0, length(x))
  m <- x >= 0
  y[m] <- x[m] / s^2 * exp(-x[m]^2 / (2 * s^2))
  y
}
rayleigh_p <- function(x, s) {
  y <- rep(0, length(x))
  m <- x >= 0
  y[m] <- 1 - exp(-x[m]^2 / (2 * s^2))
  y
}

# Interval notation, as the kernel reports it, into a support predicate.
in_range <- function(range) {
  ends <- as.numeric(strsplit(substr(range, 2, nchar(range) - 1), ", ")[[1]])
  lo <- if (startsWith(range, "[")) `>=` else `>`
  up <- if (endsWith(range, "]")) `<=` else `<`
  function(x) lo(x, ends[1]) & up(x, ends[2])
}

# Assert exposed function methods against R reference.
# A NULL reference marks a surface the kernel leaves undefined.
expect_function_matches <- function(name, args,
                                    ref_f, ref_prim = NULL,
                                    ref_log = log,
                                    range = "(-Inf, Inf)", discrete = FALSE,
                                    x = seq(0.1, 10, 0.1), h = 0.1) {
  fn <- do.call(browse_functions, c(list(name), args))
  expect_identical(fn$range, range)
  expect_identical(fn$is_discrete, discrete)

  expect_equal(fn$func(x), ref_f(x))

  if (is.null(ref_prim)) expect_error(fn$prim(x), "undefined method")
  else expect_equal(fn$prim(x), ref_prim(x))

  keep <- in_range(range)(x)
  null_move <- fn$log_ratio(x, x)
  expect_equal(null_move[keep], rep(0, sum(keep)))
  if (is.null(ref_log)) {
    expect_error(fn$log_func(x), "undefined method")
    expect_equal(fn$log_ratio(x, x + h), rep(0, length(x)))
  } else {
    expect_equal(fn$log_func(x), ref_log(fn$func(x)))
    expect_equal(fn$log_ratio(x, x + h), log(ref_f(x + h) / ref_f(x)))
  }
}

# Assert exposed distribution methods against R reference.
# The reported support gates the values the null move is defined on.
expect_distribution_matches <- function(name, args,
                                        ref_d, ref_p,
                                        range = "(-Inf, Inf)", discrete = FALSE,
                                        x = seq(-20, 20, 1), h = pi) {
  distr <- do.call(browse_distributions, c(list(name), args))
  expect_identical(distr$range, range)
  expect_identical(distr$is_discrete, discrete)

  expect_equal(distr$pdf(x), ref_d(x))
  expect_equal(distr$cdf(x), ref_p(x))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + h), log(distr$pdf(x + h) / distr$pdf(x)))
  expect_equal(distr$log_ratio(x + h, x), log(distr$pdf(x) / distr$pdf(x + h)))
  keep <- in_range(range)(x)
  null_move <- distr$log_ratio(x, x)
  expect_equal(null_move[keep], rep(0, sum(keep)))
  expect_true(all(is.na(null_move[!keep])))
}

################################################################################
####################################  TESTS  ###################################
################################################################################

test_that("Available names listing is stable", {
  # "uniform" is carried by both distribution registries and must be listed once
  expect_snapshot(browse_functions())
  expect_snapshot(browse_distributions())
  expect_snapshot(browse_infectivity_profiles())
})

test_that("Expose function has no leaks", {
  mem_used_mb <- function() gc()["Vcells", "(Mb)"]

  fun <- browse_functions("sin")
  rm(fun)
  invisible(gc())
  mem_before <- mem_used_mb()

  # heavy allocation loop
  for (i in seq_len(200)) {
    fun <- browse_functions("sin")
    rm(fun)
    if (i %% 20 == 0) gc()
  }
  invisible(gc())
  expect_lt(mem_used_mb() - mem_before, 5)
})

######################################
## BROWSE INTERNAL : METHOD SURFACE ##
######################################

test_that("Exposed fields are read only", {
  distr <- browse_distributions("normal", mean = 0, sd = 1)
  expect_error(distr$range <- "[0, 1]", "read only")
  expect_error(distr$args <- c(mean = 1), "read only")
  expect_error(distr$is_discrete <- TRUE, "read only")
})

test_that("Exposed methods carry R formals", {
  distr <- browse_distributions("lognormal", meanlog = 0, sdlog = 1)
  distr$pdf() |>
    expect_error('argument "x" is missing')
  distr$log_ratio(1) |>
    expect_error('argument "xn" is missing')
  distr$snake() |>
    expect_error("attempt to apply non-function")
  # each family answers for its own methods only
  browse_functions("sin")$inst() |>
    expect_error("attempt to apply non-function")
})

test_that("Truncation and location shift the reported support", {
  # bounds read in the same frame as the values compared against them
  ocelot <- browse_distributions("gamma", shape = 2, scale = 1, loc = 3)
  expect_equal(ocelot$range, "[3, Inf)")
  snake <- browse_distributions("normal", mean = 0, sd = 1,
                                min = -1, max = 2, loc = 10)
  expect_equal(snake$range, "[9, 12]")
  # arguments report the parametrization, location stays unreported
  distr <- browse_distributions("poisson", lambda = 3, loc = 1)
  expect_equal(distr$args, c(lambda = 3))
})

#################################
## BROWSE INTERNAL : FUNCTIONS ##
#################################

test_that("Function arguments are validated", {
  x <- seq(0.1, 10, 0.1)
  fn <- browse_functions("pow", x2 = 3)
  expect_equal(fn$args, c(x2 = 3))
  fn$log_ratio(x, x[-1]) |>
    expect_error("[Ss]ize")

  browse_functions("pow")$func(x) |>
    expect_error("Missing non default")
  browse_functions("pow", 5) |>
    expect_error("Unnamed argument")
  browse_functions("pow", x2 = 5, y2 = 3) |>
    expect_error("Unused argument.*y2.*[Dd]id you mean")
  browse_functions("pow", discr = TRUE)$func(x) |>
    expect_error("undefined method")
})

test_that("Location shifts the function argument", {
  x <- seq(-10, 10, 0.1)
  loc_val <- 1.5
  expect_equal(browse_functions("exp", loc = loc_val)$func(x), exp(x - loc_val))
  expect_equal(browse_functions("exp", loc = loc_val)$prim(x), exp(x - loc_val))
  # a quarter-turn shift turns the sine into the cosine
  expect_equal(browse_functions("sin", loc = -0.5 * pi)$func(x), cos(x))
})

test_that("Transcendental functions match R reference", {
  expect_function_matches(
    "pow", list(x2 = 3),
    function(x) x^3,
    function(x) x^4 / 4
  )
  # x2 = -1 primitive uses log(|x|)
  expect_function_matches(
    "pow", list(x2 = -1),
    function(x) 1 / x,
    function(x) log(abs(x))
  )
  # x2 = 0.5 equals sqrt
  expect_function_matches(
    "pow", list(x2 = 0.5),
    sqrt,
    function(x) x^1.5 / 1.5
  )
  expect_function_matches(
    "sqrt", list(),
    sqrt,
    function(x) x^1.5 / 1.5,
    range = "[0, Inf)"
  )
  expect_function_matches(
    "exp", list(),
    exp,
    exp,
    x = seq(-5, 5, 0.1), h = pi
  )
  expect_function_matches(
    "log", list(),
    log,
    function(x) x * log(x) - x,
    range = "(0, Inf)",
    x = seq(1.1, 20, 0.1), h = pi
  )
  expect_function_matches(
    "sin", list(),
    sin,
    function(x) -cos(x),
    x = seq(0.1, pi - 0.1, 0.05)
  )
  expect_function_matches(
    "cos", list(),
    cos,
    sin,
    x = seq(-pi / 2 + 0.1, pi / 2 - 0.1, 0.05)
  )
  expect_function_matches(
    "tan", list(),
    tan, # NOTE: not supporting a primitive for the tangent
    x = seq(0.1, pi / 2 - 0.1, 0.05)
  )
})

test_that("Arithmetic kernels have no log surface", {
  x <- seq(-5, 5, 0.5)
  fold <- c(3.3, 5.5, 6.6)
  # NOTE: avoid half-values: C++ std::round uses round-half-away-from-zero
  # while R's round() uses round-half-to-even (crazy...)
  x_step <- c(-3, -2.7, -1.3, -1, -0.3, 0, 0.4, 1.3, 2.7, 4)

  expect_function_matches(
    "minus", list(),
    function(x) -x,
    ref_log = NULL, x = x, h = pi
  )
  expect_function_matches(
    "sub", list(x2 = 3),
    function(x) x - 3,
    ref_log = NULL, x = x, h = pi
  )
  expect_function_matches(
    "sub", list(x2 = 0),
    function(x) x,
    ref_log = NULL, x = x, h = pi
  )
  expect_function_matches(
    "div", list(x2 = 4),
    function(x) x / 4,
    ref_log = NULL, x = x, h = pi
  )
  expect_function_matches(
    "div", list(x2 = 1),
    function(x) x,
    ref_log = NULL, x = x, h = pi
  )
  # integers are fixed points of the real-to-integer kernels
  expect_function_matches(
    "round", list(),
    round,
    ref_log = NULL, discrete = TRUE, x = x_step, h = pi
  )
  expect_function_matches(
    "ceil", list(),
    ceiling,
    ref_log = NULL, discrete = TRUE, x = x_step, h = pi
  )
  expect_function_matches(
    "floor", list(),
    floor,
    ref_log = NULL, discrete = TRUE, x = x_step, h = pi
  )
  # variadic kernels fold their constant arguments with x
  expect_function_matches(
    "sum", list(fold),
    function(x) x + sum(fold),
    ref_log = NULL, x = x, h = pi
  )
  expect_function_matches(
    "prod", list(fold),
    function(x) x * prod(fold),
    ref_log = NULL, x = x, h = pi
  )
  expect_function_matches(
    "min", list(fold),
    function(x) pmin(x, min(fold)),
    ref_log = NULL, x = x, h = pi
  )
  expect_function_matches(
    "max", list(fold),
    function(x) pmax(x, max(fold)),
    ref_log = NULL, x = x, h = pi
  )
})

#####################################
## BROWSE INTERNAL : DISTRIBUTIONS ##
#####################################

test_that("Truncated binomial is correct", {
  x <- seq(1, 10, 1)
  args <- list(size = 12, prob = 0.2)
  distr <- browse_distributions("binomial", args, min = 1)
  shift <- stats::pbinom(0, args$size, args$prob)
  distr$pdf(x) |>
    expect_equal(stats::dbinom(x, args$size, args$prob) / (1 - shift))
  distr$cdf(x) |>
    expect_equal((stats::pbinom(x, args$size, args$prob) - shift) / (1 - shift))
})

test_that("Truncated discrete mass is normalised on the bins inside the range", {
  # the window closes on the last bin inside the range, not one atom past it
  k <- 2:6
  distr <- browse_distributions("poisson", lambda = 3, min = min(k), max = max(k))
  norm <- diff(stats::ppois(c(min(k) - 1, max(k)), 3))
  expect_equal(distr$pdf(k), stats::dpois(k, 3) / norm)
  expect_equal(sum(distr$pdf(k)), 1)
  expect_equal(distr$cdf(max(k)), 1)

  # left-discretized: bin i carries [i, i + 1), so the range holds 5 bins
  i <- 1:5
  distr <- browse_distributions("gamma", shape = 2, scale = 1,
                                min = min(i), max = max(i), discr = TRUE)
  norm <- diff(stats::pgamma(c(min(i), max(i) + 1), shape = 2, scale = 1))
  expect_equal(distr$pdf(i),
               diff(stats::pgamma(c(i, max(i) + 1), shape = 2, scale = 1)) / norm)
  expect_equal(sum(distr$pdf(i)), 1)
  expect_equal(distr$cdf(max(i)), 1)
})

test_that("Truncated discrete rand() reaches every bin of the support", {
  # both sampling branches must reproduce the mass they are scored against
  expect_draws_match_pdf <- function(distr, k, n = 2e4) {
    freq <- as.vector(table(factor(distr$rand(n), levels = k))) / n
    expect_equal(freq, distr$pdf(k), tolerance = 0.05)
  }

  set.seed(37L)
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
  expect_identical(max(distr$rand(2e4L)), max(i))
})

test_that("Discrete uniform spreads flat mass over the whole values", {
  k <- 0:10
  distr <- browse_distributions("uniform", min = min(k), max = max(k), discr = TRUE)

  # every whole value in the range carries the same mass
  expect_equal(distr$pdf(k), rep(1 / length(k), length(k)))
  expect_equal(sum(distr$pdf(k)), 1)
  expect_equal(distr$log_pdf(k), rep(-log(length(k)), length(k)))
  expect_equal(distr$pdf(c(-1, 11)), c(0, 0))

  # the cumulative closes on the last atom, and is the running mass
  expect_equal(distr$cdf(k), cumsum(distr$pdf(k)))
  expect_equal(distr$cdf(max(k)), 1)

  # flat mass leaves every move even
  expect_equal(distr$log_ratio(c(0, 5), c(10, 6)), c(0, 0))

  # the quantile is the first atom whose cumulative reaches p
  p <- c(0.01, 0.25, 0.5, 0.75, 0.99)
  expect_equal(distr$quantile(p), vapply(p, function(q) k[which(cumsum(distr$pdf(k)) >= q)[1]], 0))

  # the default range is a coin flip
  coin <- browse_distributions("uniform", discr = TRUE)
  expect_equal(coin$pdf(0:1), c(0.5, 0.5))
  expect_equal(coin$cdf(0:1), c(0.5, 1))

  # 'loc' carries the atoms without reshaping them
  browse_distributions("uniform", min = 0, max = 10, loc = 100, discr = TRUE)$pdf(100 + k) |>
    expect_equal(distr$pdf(k))
})

test_that("An argument off the whole values is refused where one is required", {
  # the flat lattice law owns every whole value in its range
  browse_distributions("uniform", min = 0.5, max = 3.2, discr = TRUE) |>
    expect_error("min.*uniform.*discrete")
  browse_distributions("uniform", min = 0, max = 3.2, discr = TRUE) |>
    expect_error("max.*uniform.*discrete")

  # a count of trials or successes is whole whatever the support looks like
  browse_distributions("binomial", size = 10.5, prob = 0.5) |>
    expect_error("size.*binomial.*discrete")
  browse_distributions("negative_binomial", size = 5.5, prob = 0.5) |>
    expect_error("size.*negative_binomial.*discrete")
})

test_that("Standard distributions match R reference d/p/log_pdf/log_ratio", {
  expect_distribution_matches(
    "normal", list(mean = 2, sd = 3),
    function(x) stats::dnorm(x, 2, 3),
    function(x) stats::pnorm(x, 2, 3)
  )
  expect_distribution_matches(
    "uniform", list(min = -1, max = 5),
    function(x) stats::dunif(x, -1, 5),
    function(x) stats::punif(x, -1, 5),
    range = "[-1, 5]"
  )
  expect_distribution_matches(
    "lognormal", list(meanlog = 1, sdlog = 0.5),
    function(x) stats::dlnorm(x, 1, 0.5),
    function(x) stats::plnorm(x, 1, 0.5),
    range = "[0, Inf)"
  )
  expect_distribution_matches(
    "beta", list(shape1 = 2, shape2 = 1.5),
    function(x) stats::dbeta(x, 2, 1.5),
    function(x) stats::pbeta(x, 2, 1.5),
    range = "(0, 1)",
    x = seq(-10, 10, 0.1), h = 0.1
  )
  expect_distribution_matches(
    "gamma", list(shape = 2, scale = 1.2),
    function(x) stats::dgamma(x, shape = 2, scale = 1.2),
    function(x) stats::pgamma(x, shape = 2, scale = 1.2),
    range = "[0, Inf)"
  )
  expect_distribution_matches(
    "exponential", list(rate = 0.25),
    function(x) stats::dexp(x, 0.25),
    function(x) stats::pexp(x, 0.25),
    range = "[0, Inf)"
  )
  expect_distribution_matches(
    "poisson", list(lambda = 5),
    function(x) stats::dpois(x, 5),
    function(x) stats::ppois(x, 5),
    range = "[0, Inf)", discrete = TRUE
  )
  expect_distribution_matches(
    "binomial", list(size = 10, prob = 0.5),
    function(x) stats::dbinom(x, 10, 0.5),
    function(x) stats::pbinom(x, 10, 0.5),
    range = "[0, 10]", discrete = TRUE
  )
  expect_distribution_matches(
    "negative_binomial", list(size = 10, prob = 0.5),
    function(x) stats::dnbinom(x, 10, 0.5),
    function(x) stats::pnbinom(x, 10, 0.5),
    range = "[0, Inf)", discrete = TRUE
  )
  expect_distribution_matches(
    "geometric", list(prob = 0.5),
    function(x) stats::dgeom(x, 0.5),
    function(x) stats::pgeom(x, 0.5),
    range = "[0, Inf)", discrete = TRUE
  )
  expect_distribution_matches(
    "weibull", list(shape = 2, scale = 1),
    function(x) stats::dweibull(x, 2, 1),
    function(x) stats::pweibull(x, 2, 1),
    range = "[0, Inf)"
  )
  expect_distribution_matches(
    "rayleigh", list(scale = 2),
    function(x) rayleigh_d(x, 2),
    function(x) rayleigh_p(x, 2),
    range = "[0, Inf)",
    x = seq(0, 10, 0.5)
  )
})

test_that("Categorical distribution is correct", {
  # unnormalized weights still give a normalized pmf
  args <- list(fortune = 2, fatman = 1, vamp = 7)
  distr <- browse_distributions("categorical", args)
  pmf <- unname(unlist(args)) / sum(unlist(args))
  expect_equal(distr$pdf(0:2), pmf)
  expect_equal(distr$pdf(c(-1, 3)), c(0, 0))
  expect_equal(distr$log_pdf(0:2), log(pmf))

  # cumulative pmf
  expect_equal(distr$cdf(0:2), cumsum(pmf))
  expect_equal(distr$cdf(c(-1, 5)), c(0, 1))

  # log-ratio
  distr$log_ratio(c(0, 1), c(1, 2)) |>
    expect_equal(log(pmf[c(2, 3)] / pmf[c(1, 2)]))

  browse_distributions("categorical") |>
    expect_error("no weights")
  browse_distributions("categorical", peace = 1, walker = -2) |>
    expect_error("non-strictly positive")
})

test_that("Default parameters match R reference for d", {
  x <- seq(-20, 20, 1)
  expect_equal(browse_distributions("normal")$pdf(x), stats::dnorm(x))
  expect_equal(browse_distributions("uniform")$pdf(x), stats::dunif(x))
  expect_equal(browse_distributions("lognormal")$pdf(x), stats::dlnorm(x))
  expect_equal(browse_distributions("exponential")$pdf(x), stats::dexp(x))
})

test_that("Normal both-sided truncation is correct", {
  low <- 0
  up <- 5
  x <- seq(-20, 20, 1)
  distr <- browse_distributions("normal", mean = 2, sd = 3, min = low, max = up)
  ref_d <- function(x) stats::dnorm(x, 2, 3)
  ref_p <- function(x) stats::pnorm(x, 2, 3)
  zz <- ref_p(up) - ref_p(low)

  expect_equal(distr$pdf(x), ifelse(x >= low & x <= up, ref_d(x), 0) / zz)
  expect_equal(
    distr$cdf(x),
    ifelse(x < low, 0,
      ifelse(x > up, 1, (ref_p(x) - ref_p(low)) / zz)
    )
  )
})

test_that("Gamma accepts shape/scale and mean/sd parametrizations", {
  browse_distributions("gamma", mean = 2, scale = 1.1) |>
    expect_error("mix")
  browse_distributions("gamma", shape = 2, sd = 1.1) |>
    expect_error("mix")

  m <- 2.5
  s <- 1.3
  shape <- m * m / (s * s)
  scale <- s * s / m
  x <- seq(0.1, 12, 0.3)
  d_can <- browse_distributions("gamma", shape = shape, scale = scale)
  d_alt <- browse_distributions("gamma", mean = m, sd = s)
  expect_equal(d_can$pdf(x), d_alt$pdf(x))
  expect_equal(d_can$cdf(x), d_alt$cdf(x))
  expect_equal(d_can$log_pdf(x), d_alt$log_pdf(x))

  # mean/sd matches stats::dgamma with derived shape/scale
  expect_equal(d_alt$pdf(x), stats::dgamma(x, shape = shape, scale = scale))
})

test_that("Lognormal accepts meanlog/sdlog and mean/sd parametrizations", {
  browse_distributions("lognormal", mean = 2, sdlog = 1.1) |>
    expect_error("mix")
  browse_distributions("lognormal", meanlog = 2, sd = 1.1) |>
    expect_error("mix")

  m <- 2.5
  s <- 1.3
  sdlog <- sqrt(log1p((s * s) / (m * m)))
  meanlog <- log(m) - 0.5 * sdlog * sdlog
  x <- seq(0.1, 12, 0.3)
  d_can <- browse_distributions("lognormal", meanlog = meanlog, sdlog = sdlog)
  d_alt <- browse_distributions("lognormal", mean = m, sd = s)
  expect_equal(d_can$pdf(x), d_alt$pdf(x))
  expect_equal(d_can$cdf(x), d_alt$cdf(x))
  expect_equal(d_can$log_pdf(x), d_alt$log_pdf(x))

  # mean/sd matches stats::dlnorm with derived meanlog/sdlog
  expect_equal(d_alt$pdf(x), stats::dlnorm(x, meanlog = meanlog, sdlog = sdlog))
})

test_that("CDF and quantile make sense", {
  x <- seq(-10, 20, 0.1)
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  distrs <- list(
    list(name = "uniform", min = -3, max = 2),
    list(name = "normal", mean = 3, sd = 2),
    list(name = "gamma", shape = 2, scale = 1.5),
    list(name = "gamma", mean = 2.2, sd = sqrt(1.1)),
    list(name = "exponential", rate = 0.4),
    list(name = "lognormal", meanlog = 1, sdlog = 0.5),
    list(name = "weibull", shape = 2, scale = 2),
    list(name = "rayleigh", scale = 2)
  )
  for (d in distrs) {
    distr <- do.call(browse_distributions, d)
    cdf_vals <- distr$cdf(x)
    pdf_vals <- distr$pdf(x)
    expect_true(all(diff(cdf_vals) >= -1e-12))
    expect_true(all(cdf_vals >= 0 & cdf_vals <= 1))
    expect_true(all(pdf_vals >= 0))
    expect_equal(distr$cdf(distr$quantile(p)), p)
  }
})

test_that("CDF is the integral of PDF for continuous distributions", {
  # gamma: integrate pdf from 0 to 15, expect cdf(15)
  x <- seq(0.001, 15, length.out = 3000)
  distr <- browse_distributions("gamma", shape = 2, scale = 1.5)
  expect_equal(trapezoid(x, distr$pdf(x)), distr$cdf(tail(x, 1)),
    tolerance = 1e-3
  )

  # lognormal: integrate pdf from 0 to 20, expect cdf(20)
  x <- seq(0.001, 20, length.out = 3000)
  distr <- browse_distributions("lognormal", meanlog = 1, sdlog = 0.5)
  expect_equal(trapezoid(x, distr$pdf(x)), distr$cdf(tail(x, 1)),
    tolerance = 1e-3
  )

  # rayleigh: integrate pdf from 0 to 12, expect cdf(12)
  x <- seq(0.001, 12, length.out = 3000)
  distr <- browse_distributions("rayleigh", scale = 2)
  expect_equal(trapezoid(x, distr$pdf(x)), distr$cdf(tail(x, 1)),
    tolerance = 1e-3
  )
})

test_that("CDF equals cumsum of PDF for discrete distributions", {
  x <- seq(0, 20, 1)

  # poisson: CDF(k) = sum of PMF from 0 to k
  distr <- browse_distributions("poisson", lambda = 5)
  expect_equal(distr$cdf(x), cumsum(distr$pdf(x)), tolerance = 1e-10)

  # geometric: CDF(k) = sum of PMF from 0 to k
  distr <- browse_distributions("geometric", prob = 0.3)
  expect_equal(distr$cdf(x), cumsum(distr$pdf(x)), tolerance = 1e-10)
})

test_that("Non-integer x maps to floor for Discrete::Defined distributions", {
  lambda <- 3
  distr <- browse_distributions("poisson", lambda = lambda)
  x_floor <- seq(0, 9, 1)
  for (k in x_floor) {
    # pdf at k and any value in [k, k+1) should be identical
    expect_equal(distr$pdf(k), distr$pdf(k + 0.7))
    expect_equal(distr$cdf(k), distr$cdf(k + 0.7))
  }
})

test_that("Lognormal truncation is correct", {
  min_val <- 1
  max_val <- 5
  x <- seq(0, 8, 0.5)
  distr <- browse_distributions("lognormal",
    meanlog = 1, sdlog = 0.5,
    min = min_val, max = max_val
  )

  lnpdf <- function(x) stats::dlnorm(x, meanlog = 1, sdlog = 0.5)
  lncdf <- function(x) stats::plnorm(x, meanlog = 1, sdlog = 0.5)
  zz <- lncdf(max_val) - lncdf(min_val)

  expec_pdf <- ifelse(x >= min_val & x <= max_val, lnpdf(x), 0) / zz
  expec_cdf <- ifelse(
    x < min_val, 0,
    ifelse(x > max_val, 1, (lncdf(x) - lncdf(min_val)) / zz)
  )

  expect_equal(distr$pdf(x), expec_pdf)
  expect_equal(distr$cdf(x), expec_cdf)
})

test_that("Exponential truncation at max is correct", {
  rate <- 0.3
  max_val <- 5
  x <- seq(0, 8, 0.5)
  distr <- browse_distributions("exponential", rate = rate, max = max_val)

  exppdf <- function(x) stats::dexp(x, rate = rate)
  expcdf <- function(x) stats::pexp(x, rate = rate)
  zz <- expcdf(max_val)

  expec_pdf <- ifelse(x >= 0 & x <= max_val, exppdf(x), 0) / zz
  expec_cdf <- ifelse(x < 0, 0, ifelse(x > max_val, 1, expcdf(x) / zz))

  expect_equal(distr$pdf(x), expec_pdf)
  expect_equal(distr$cdf(x), expec_cdf)
})

test_that("Discretized exponential uses left-discretization", {
  rate <- 0.4
  x <- seq(0, 5, 1)
  distr <- browse_distributions("exponential", rate = rate, discr = TRUE)

  # left discretization: pdf(k) = pexp(k+1) - pexp(k)
  expcdf <- function(x) stats::pexp(x, rate = rate)
  expected_pmf <- expcdf(x + 1) - expcdf(x)
  expect_equal(distr$pdf(x), expected_pmf)
  expect_equal(distr$cdf(x), cumsum(expected_pmf), tolerance = 1e-10)
})

test_that("Discretized gamma uses left-discretization", {
  args <- list(shape = 2, scale = 1.5)
  x <- seq(0, 10, 1)
  distr <- browse_distributions("gamma", args, discr = TRUE)

  gammcdf <- function(x) do.call(stats::pgamma, c(list(q = x), args))
  expect_equal(distr$pdf(x), gammcdf(x + 1) - gammcdf(x))
})

test_that("Discretized normal uses center-discretization", {
  mu <- 5
  sd <- 2
  x <- seq(0, 10, 1)
  distr <- browse_distributions("normal", mean = mu, sd = sd, discr = TRUE)

  # center discretization: pdf(k) = pnorm(k+0.5) - pnorm(k-0.5)
  normcdf <- function(x) stats::pnorm(x, mean = mu, sd = sd)
  expect_equal(distr$pdf(x), normcdf(x + 0.5) - normcdf(x - 0.5))
})

test_that("Discretized normal stays finite far in the upper tail", {
  # far above the mean both CDFs round to 1 and the naive difference cancels
  # to 0 (-Inf log-pmf); the complementary-CDF branch keeps it finite.
  x <- seq(9, 30, 1)
  distr <- browse_distributions("normal", mean = 0, sd = 1, discr = TRUE)

  # the bug this guards against: naive CDF difference underflows to exactly 0
  naive <- stats::pnorm(x + 0.5) - stats::pnorm(x - 0.5)
  expect_true(any(naive == 0))

  # stable reference via upper-tail (complementary) CDF
  ccdf <- function(x) stats::pnorm(x, lower.tail = FALSE)
  expected_pmf <- ccdf(x - 0.5) - ccdf(x + 0.5)

  expect_true(all(is.finite(distr$log_pdf(x))))
  expect_equal(distr$pdf(x), expected_pmf)
  expect_equal(distr$log_pdf(x), log(expected_pmf))

  # symmetry: reflecting about the mean matches the lower tail
  expect_equal(distr$pdf(x), distr$pdf(-x))
})

test_that("Discretized lognormal keeps the mass of the origin bin", {
  args <- list(meanlog = 0, sdlog = 1)
  x <- 0:5000
  distr <- browse_distributions("lognormal", args, discr = TRUE)

  lncdf <- function(x) do.call(stats::plnorm, c(list(q = x), args))
  expect_equal(distr$pdf(x), lncdf(x + 1) - lncdf(x))
  expect_equal(sum(distr$pdf(x)), 1)

  # the pmf and the cdf agree on the origin bin
  expect_gt(distr$pdf(0), 0)
  expect_equal(distr$cdf(0:3), cumsum(distr$pdf(0:3)))

  # in continuous time the origin is a zero-density point, not a NaN
  contin <- browse_distributions("lognormal", args)
  expect_equal(contin$pdf(0), 0)
  expect_equal(contin$log_pdf(0), -Inf)
  expect_equal(contin$log_ratio(0, 0), 0)
  expect_equal(contin$log_ratio(0, 2), Inf)
  expect_equal(contin$log_ratio(2, 0), -Inf)
})

test_that("Discretized lognormal stays finite far in the upper tail", {
  # left-discretized, heavy upper tail: exercises the gamma/erfc complement path
  args <- list(meanlog = 0, sdlog = 0.5)
  x <- seq(60, 200, 5)
  distr <- browse_distributions("lognormal", args, discr = TRUE)

  lncdf <- function(x) do.call(stats::plnorm, c(list(q = x), args))
  expect_true(any(lncdf(x + 1) - lncdf(x) == 0)) # naive underflows

  ccdf <- function(x) {
    do.call(stats::plnorm, c(list(q = x, lower.tail = FALSE), args))
  }
  expected_pmf <- ccdf(x) - ccdf(x + 1)

  expect_true(all(is.finite(distr$log_pdf(x))))
  expect_equal(distr$pdf(x), expected_pmf)
})

test_that("Discretized gamma stays finite far in the upper tail", {
  # left-discretized; exercises the boost gamma_q complement override
  args <- list(shape = 2, scale = 1.5)
  x <- seq(60, 120, 5)
  distr <- browse_distributions("gamma", args, discr = TRUE)

  gcdf <- function(x) do.call(stats::pgamma, c(list(q = x), args))
  expect_true(any(gcdf(x + 1) - gcdf(x) == 0)) # naive underflows

  gccdf <- function(x) {
    do.call(stats::pgamma, c(list(q = x, lower.tail = FALSE), args))
  }
  expected_pmf <- gccdf(x) - gccdf(x + 1)

  expect_true(all(is.finite(distr$log_pdf(x))))
  expect_equal(distr$pdf(x), expected_pmf)
})

test_that("rand() draws are reproducible and statistically plausible", {
  d <- browse_distributions("normal", mean = 2, sd = 3)
  set.seed(42L)
  x1 <- d$rand(50L)
  set.seed(42L)
  x2 <- d$rand(50L)
  expect_identical(x1, x2) # same R seed -> same draws
  expect_length(x1, 50L)
  # different seed -> different draws
  set.seed(43L)
  expect_false(identical(x1, d$rand(50L)))

  # mean / sd of 20k draws should be close to the truth
  set.seed(7L)
  big <- d$rand(20000L)
  expect_equal(mean(big), 2, tolerance = 0.05)
  expect_equal(sd(big), 3, tolerance = 0.05)

  # without set.seed, successive calls are non-deterministic
  expect_false(identical(d$rand(10L), d$rand(10L)))

  # invalid n
  d$rand(-1L) |>
    expect_error("non-negative")
})

test_that("rand() honours location shift (relocate_ Active branch)", {
  loc_val <- 5
  d <- browse_distributions("normal", mean = 0, sd = 1, loc = loc_val)
  expect_equal(d$pdf(loc_val), stats::dnorm(0))
  expect_equal(mean(d$rand(20000L)), loc_val, tolerance = 0.05)
})

test_that("cdf honours location shift above the upper bound", {
  args <- list(min = 2, max = 4, loc = 10)
  shifted <- browse_distributions("uniform", args)
  args <- as.list(unlist(args) + args$loc)
  args$loc <- NULL
  plain <- browse_distributions("uniform", args)
  x <- seq(9, 17, 0.5)
  expect_equal(shifted$pdf(x), plain$pdf(x))
  expect_equal(shifted$cdf(x), plain$cdf(x))

  # shifting by 'loc' does not tocuh arguments
  args <- list(shape = 2, scale = 1, loc = 3, max = 4)
  shifted <- browse_distributions("gamma", args)
  args$loc <- NULL
  plain <- browse_distributions("gamma", args)
  expect_equal(shifted$cdf(x), plain$cdf(x - 3))
  expect_equal(shifted$cdf(20), 1)
})

test_that("rand() on truncated distribution lands inside the bounds", {
  d <- browse_distributions("normal", mean = 0, sd = 1, min = -1, max = 2)
  x <- d$rand(5000L)
  expect_true(all(x >= -1 & x <= 2))
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
    expect_error("[Ss]et.*`loc`.*integer")
  browse_distributions("gamma", shape = 2, discr = TRUE, loc = -1.5) |>
    expect_error("[Ss]et.*`loc`.*integer")
  # an integer shift, and a fractional bound that only truncates, are fine
  browse_distributions("poisson", lambda = 3, loc = -2) |> expect_no_error()
  browse_distributions("poisson", lambda = 3, min = 0.5) |> expect_no_error()
  # a continuous definition is free to sit anywhere
  browse_distributions("gamma", shape = 2, scale = 1, loc = 0.5) |>
    expect_no_error()
})

test_that("quantile returns NA outside [0, 1]", {
  distr <- browse_distributions("normal", mean = 0, sd = 1)
  expect_true(is.na(distr$quantile(-0.1)))
  expect_true(is.na(distr$quantile(1.5)))
  expect_true(all(is.na(distr$quantile(c(-2, -0.001, 1.001, 7.7)))))
})

test_that("quantile at boundaries returns +-Inf when domain is unbounded", {
  distr <- browse_distributions("normal", mean = 0, sd = 1)
  expect_equal(distr$quantile(0), -Inf)
  expect_equal(distr$quantile(1), Inf)
})

test_that("Discretized normal methods go through center-discretization", {
  args <- list(mean = 5, sd = 2)
  distr <- browse_distributions("normal", args, discr = TRUE)
  normcdf <- function(x) do.call(stats::pnorm, c(list(q = x), args))

  x <- seq(0, 10, 1)
  expected_pmf <- normcdf(x + 0.5) - normcdf(x - 0.5)

  # log_func_ Center branch
  expect_equal(distr$log_pdf(x), log(expected_pmf))
  # primitive_ Center branch (cdf is the primitive)
  expect_equal(distr$cdf(x), normcdf(x + 0.5))
  # log_ratio_ Center branch with distinct rounded values
  xn <- x + 1
  distr$log_ratio(x, xn) |>
    expect_equal(log(distr$pdf(xn) / distr$pdf(x)))
  # log_ratio when rounded values coincide should be 0
  expect_equal(distr$log_ratio(x, x + 0.1), rep(0, length(x)))
})

test_that("Upper truncatation of normal distribution is correct", {
  args <- list(mean = 2, sd = 1.5)
  max_val <- 4
  distr <- browse_distributions("normal", args, max = max_val)

  normpdf <- function(x) do.call(stats::dnorm, c(list(x = x), args))
  normcdf <- function(x) do.call(stats::pnorm, c(list(q = x), args))
  zz <- normcdf(max_val)

  x <- seq(-3, 3, 0.5)
  expec_pdf <- ifelse(x > max_val, 0, normpdf(x)) / zz
  expec_cdf <- ifelse(x > max_val, 1, normcdf(x) / zz)

  expect_equal(distr$pdf(x), expec_pdf)
  expect_equal(distr$cdf(x), expec_cdf)

  # case of center discretization with truncation at max
  distr <- browse_distributions("normal", args, max = max_val, discr = TRUE)
  zz <- normcdf(floor(max_val) + 0.5)

  x <- seq(-3, 5, 1)
  expec_pmf <- ifelse(x > max_val, 0, normcdf(x + 0.5) - normcdf(x - 0.5)) / zz

  expect_equal(distr$pdf(x), expec_pmf)
})

test_that("quantile matches R reference for standard distributions", {
  p <- seq(0, 1, 0.1)
  cases <- list(
    list(args = list("normal", mean = 2, sd = 3),
         ref = function(p) stats::qnorm(p, 2, 3)),
    list(args = list("lognormal", meanlog = 1, sdlog = 0.5),
         ref = function(p) stats::qlnorm(p, 1, 0.5)),
    list(args = list("gamma", shape = 2, scale = 1.5),
         ref = function(p) stats::qgamma(p, shape = 2, scale = 1.5)),
    list(args = list("exponential", rate = 0.3),
         ref = function(p) stats::qexp(p, 0.3)),
    list(args = list("beta", shape1 = 2, shape2 = 3),
         ref = function(p) stats::qbeta(p, 2, 3)),
    list(args = list("poisson", lambda = 5),
         ref = function(p) stats::qpois(p, 5)),
    list(args = list("binomial", size = 10, prob = 0.4),
         ref = function(p) stats::qbinom(p, 10, 0.4)),
    list(args = list("geometric", prob = 0.3),
         ref = function(p) stats::qgeom(p, 0.3)),
    list(args = list("negative_binomial", size = 5, prob = 0.4),
         ref = function(p) stats::qnbinom(p, 5, 0.4))
  )
  for (case in cases) {
    distr <- do.call(browse_distributions, case$args)
    expect_equal(distr$quantile(p), case$ref(p), info = case$args[[1]])
  }
})

test_that("quantile is truncation-aware", {
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  cases <- list(
    list(args = list("normal", mean = 0, sd = 1, min = 0),
         ref = function(p) qnorm(pnorm(0) + p * (1 - pnorm(0)))),
    list(args = list("normal", mean = 3, sd = 2, min = 1, max = 5),
         ref = function(p) {
           qnorm(pnorm(1, 3, 2) + p * (pnorm(5, 3, 2) - pnorm(1, 3, 2)), 3, 2)
         }),
    list(args = list("exponential", rate = 1, max = 3),
         ref = function(p) qexp(p * pexp(3))),
    list(args = list("gamma", shape = 2, scale = 1.5, min = 1, max = 8),
         ref = function(p) {
           qgamma(pgamma(1, shape = 2, scale = 1.5) +
                    p * (pgamma(8, 2, scale = 1.5) - pgamma(1, 2, scale = 1.5)),
                  shape = 2, scale = 1.5)
         })
  )
  for (case in cases) {
    distr <- do.call(browse_distributions, case$args)
    expect_equal(distr$quantile(p), case$ref(p),
      tolerance = 1e-10, info = case$args[[1]]
    )
    expect_true(all(distr$quantile(p) >= case$args$min %||% -Inf),
      info = paste(case$args[[1]], "above min")
    )
    expect_true(all(distr$quantile(p) <= case$args$max %||% Inf),
      info = paste(case$args[[1]], "below max")
    )
  }
})

test_that("Beta rejects discrete option", {
  browse_distributions("beta", shape1 = 2, shape2 = 3, discr = TRUE) |>
    expect_error("invalid.*discrete")
})

test_that("Beta rejects out-of-domain truncation bounds", {
  browse_distributions("beta", shape1 = 2, shape2 = 3, min = -0.5) |>
    expect_error("min")
  browse_distributions("beta", shape1 = 2, shape2 = 3, max = 1.5) |>
    expect_error("max")
})

test_that("Beta truncation normalises correctly", {
  low <- 0.1
  up <- 0.8
  args <- list(shape1 = 2, shape2 = 5)
  x <- seq(0, 1, 0.05)
  d <- browse_distributions("beta", args, min = low, max = up)
  ref_d <- function(x) do.call(stats::dbeta, c(list(x = x), args))
  ref_p <- function(x) do.call(stats::pbeta, c(list(q = x), args))
  zz <- ref_p(up) - ref_p(low)
  expect_equal(d$pdf(x), ifelse(x >= low & x <= up, ref_d(x), 0) / zz)
})

test_that("Beta rand() produces plausible draws", {
  args <- list(shape1 = 2, shape2 = 5)
  d <- browse_distributions("beta", args)
  set.seed(99L)
  x1 <- d$rand(50L)
  set.seed(99L)
  x2 <- d$rand(50L)
  expect_identical(x1, x2)
  big <- d$rand(2000L)
  expect_true(all(big >= 0 & big <= 1))
  expect_equal(mean(big), args$shape1 / sum(unlist(args)), tolerance = 0.02)
})

test_that("Gamma pdf and log_pdf handle x == 0 across shape regimes", {
  scale <- 2
  # test case: shape < 1
  d1 <- browse_distributions("gamma", shape = 0.5, scale = scale)
  expect_equal(d1$pdf(0), Inf)
  expect_equal(d1$log_pdf(0), Inf)

  # test case: shape == 1
  d2 <- browse_distributions("gamma", shape = 1, scale = scale)
  expect_equal(d2$pdf(0), 1 / scale)
  expect_equal(d2$log_pdf(0), -log(scale))

  # test case: shape > 1
  d3 <- browse_distributions("gamma", shape = 2, scale = scale)
  expect_equal(d3$pdf(0), 0)
  expect_equal(d3$log_pdf(0), -Inf)
})

test_that("Binomial rejects out-of-domain truncation bounds", {
  browse_distributions("binomial", size = 10, prob = 0.5, min = -1) |>
    expect_error("min")
  browse_distributions("binomial", size = 10, prob = 0.5, min = 11) |>
    expect_error("max")
})

test_that("rand() is reproducible and plausible for all distributions", {
  cases <- list(
    list(args = list("rayleigh", scale = 2),
         mean = 2 * sqrt(pi / 2), int = FALSE),
    list(args = list("exponential", rate = 0.5),
         mean = 2, int = FALSE),
    list(args = list("poisson", lambda = 5),
         mean = 5, int = TRUE),
    list(args = list("binomial", size = 20, prob = 0.3),
         mean = 6, int = TRUE, max = 20),
    list(args = list("negative_binomial", size = 5, prob = 0.4),
         mean = 5 * 0.6 / 0.4, int = TRUE),
    list(args = list("geometric", prob = 0.3),
         mean = 0.7 / 0.3, int = TRUE)
  )
  for (cc in cases) {
    d <- do.call(browse_distributions, cc$args)
    nm <- cc$args[[1]]
    set.seed(42L)
    x1 <- d$rand(50L)
    set.seed(42L)
    expect_identical(d$rand(50L), x1, info = nm)
    set.seed(1L)
    big <- d$rand(2000L)
    expect_true(all(big >= 0), info = nm)
    if (cc$int) expect_true(all(big == floor(big)), info = nm)
    if (!is.null(cc$max)) expect_true(all(big <= cc$max), info = nm)
    expect_equal(mean(big), cc$mean, tolerance = 0.15, info = nm)
  }
})

############################################
## BROWSE INTERNAL : INFECTIVITY PROFILES ##
############################################

test_that("Ferretti20 infectivity profile behave has expected", {
  x <- seq(0, 15, 0.001)
  profile <- browse_infectivity_profiles("ferretti20", incub_period = 5)
  expect_equal(trapezoid(x, profile$inst(x)), profile$cumul(x)[length(x)])
  expect_false(profile$is_discrete)
  expect_equal(profile$range, "[0, Inf)")

  # left discretization works as expected
  x <- seq(1.0, 1.9, 0.1)
  profile <- browse_infectivity_profiles("ferretti20",
                                         incub_period = 5, discr = TRUE)
  expect_true(all(profile$inst(x) == profile$cumul(1)))
  expect_true(profile$is_discrete)
  expect_equal(profile$range, "[1, Inf)")
})

####################################
## BROWSE INTERNAL : PRINT METHOD ##
####################################

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
