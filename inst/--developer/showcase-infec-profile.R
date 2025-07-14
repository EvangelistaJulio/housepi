

hhsimul <- HouseholdModel$new("augm-infection-time")

# extract current loaded profile
infec <- hhsimul$infec_profile_defs

# plot instantaneous and cumulative profile with ggplot2
library(ggplot2)

# coefficient to rescale cumulative
housepi_infec_profile(infec$name, x, discr = infec$discr, args = infec$args)

rescale <- max(infec$profile$cumul) / max(infec$profile$inst)

ggplot(infec$profile, aes(x = gentime)) +
  # plot cumulative profile
  theme(axis.ticks.y.right = element_line(color = "coral"),
        axis.text.y.right = element_text(color = "coral"),
        axis.title.y.right = element_text(color = "coral")) +
  geom_line(aes(y = cumul / rescale), color = "coral") +
  # plot instantaneous
  geom_line(aes(y = inst)) +
  # apply rescaling
  scale_y_continuous(name = "Instantaneous",
                     sec.axis = sec_axis(~. * rescale, name = "Cumulative")) +
  labs(x = "Generation Time")