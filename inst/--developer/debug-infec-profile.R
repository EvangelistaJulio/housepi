library(housepi)
library(ggplot2)

# load model
hhmodel <- HouseholdModel$new("augm-onset-single-intro")
# parameter definitions
hhmodel$load_param_template()

# CHECKED: this works
#output <- hhmodel$export_infectivity_profile(c(3))

# CHECKED: this works
# hhmodel$add_infectivity_profile(name = "reparam_gamma")
# output <- hhmodel$export_infectivity_profile()

# FIXME: this combination breaks R
hhmodel$add_infectivity_profile(name = "reparam_gamma", discrete = TRUE)
output <- hhmodel$export_infectivity_profile()


# FIXME: this is not normalised
# hhmodel$add_infectivity_profile(name = "pois", min = 2,
#                                 args = c(lambda = 4))
# FIXME: this breaks R
#hhmodel$add_infectivity_profile(name = "nbinom", min =2)

#############################
# DRAFT for a show function #
#############################

if (output$discrete) {
  # FLAG
  print(paste0("norm:", sum(output$profile$inst)))
  # FLAG-END

  ggplot(output$profile, aes(x = gentime)) +
    geom_bar(aes(y = inst), stat = "identity", fill = "coral") +
    labs(title = sprintf("%s [mean=%.1f, sd=%.1f]",
                         output$name, output$mean, sqrt(output$var)),
         x = "Generation Delay",
         y = "Daily Infectivity Profile")

} else {
  # FLAG
  dx <- output$profile$gentime[2] - output$profile$gentime[1]
  f <- output$profile$inst
  print(paste0("norm:", dx * (sum(f[2:(length(f)-1)]) + 0.5 * (f[1] + f[length(f)]))))
  # FLAG-END

  coeff <- max(output$profile$cumul) / max(output$profile$inst)
  ggplot(output$profile, aes(x = gentime)) +
    theme(axis.ticks.y.right = element_line(color = "coral"),
          axis.text.y.right = element_text(color = "coral"),
          axis.title.y.right = element_text(color = "coral")) +
    geom_line(aes(y = cumul / coeff), color = "coral") +
    geom_line(aes(y = inst)) +
    scale_y_continuous(name = "Instantaneous",
                       sec.axis = sec_axis(~. * coeff, name = "Cumulative")) +
    labs(x = "Generation Time", title = "Infectivity Profile")
}


