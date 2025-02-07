
hmc <- housepi::HouseholdModelClass$new()

################################################################
################################################################
################################################################
################################################################
################################################################

# test binomial cdf #

size = 8
prob = 0.25
hmc$add_infectivity_profile(
  name = "binom", args = list(size = size, prob = prob)
)

out <- hmc$export_infectivity_profile()
plot(out[[1]], out[[3]], type="l")
abline(h = 1, lty =2, col = "red")
abline(h = 0, lty =2, col = "red")

x <- seq(min(out[[1]]), max(out[[1]]), out[[1]][2] - out[[1]][1])
y <- pbinom(x, size = size, prob = prob)
lines(x, y, lty = 3, col = "green")

################################################################
################################################################
################################################################
################################################################
################################################################


hmc$add_infectivity_profile(
  name = "lnorm",
  discrete = TRUE
)
out <- hmc$export_infectivity_profile(opts = list(2))

plot(out[[1]], out[[2]], type="l")
abline(h = 1, lty =2, col = "red")
abline(h = 0, lty =2, col = "red")

################################################################
################################################################
################################################################
################################################################
################################################################

# test regularised incomplete beta function #

shape1 = 0.02
shape2 = 105.6
out <- hmc$add_infectivity_profile(
  name = "beta", args = list(shape1 = shape1, shape2 = shape2)
  )$export_infectivity_profile()
plot(out[[1]], out[[3]], type="l")
abline(h = 1, lty =2, col = "red")
abline(h = 0, lty =2, col = "red")

x <- seq(min(out[[1]]), max(out[[1]]), 0.001)
y <- pbeta(x, shape1 = shape1, shape2 = shape2)
lines(x, y, lty =3, col = "green")



# test regularised incomplete gamma function #

shape = 0.5
scale = 10.3
out <- hmc$add_infectivity_profile(
  name = "gamma", args = list(shape = shape, scale = scale)
)$export_infectivity_profile(xlim = c(NA, 50))
plot(out[[1]], out[[3]], type="l")
abline(h = 1, lty =2, col = "red")
abline(h = 0, lty =2, col = "red")

x <- seq(min(out[[1]]), max(out[[1]]), 0.001)
y <- pgamma(x, shape = shape, scale = scale)
lines(x, y, lty =3, col = "green")

