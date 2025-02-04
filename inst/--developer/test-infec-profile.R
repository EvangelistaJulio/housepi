
hmc <- housepi::HouseholdModelClass$new("blank")

hmc$add_infectivity_profile(
  name = "lnorm", 
)

# FIXME: the string-ization of the values is still to be fixed...
out <- hmc$export_infectivity_profile(
  xlim = list(max = "5")
)

plot(out[[1]], out[[2]], type="l")

x <- seq(0, 5, 0.01)
y <- dlnorm(x)
lines(x, y, type="l", lty=2, col = "red")

# FIXME: changing profile requires to generate a new one...

#lines(out[[1]], out[[3]], type="l", col = "red")
# FIXME: beta cdf is still very wrong

