library(housepi)
library(ggplot2)

# load model
hhmodel <- HouseholdModel$new("augm-onset-single-intro")
# parameter definitions
hhmodel$load_param_template()

hhmodel$add_infectivity_profile(name = "ferretti20")
output <- hhmodel$export_infectivity_profile(c(3))

#############################
# DRAFT for a show function #
#############################

df <- as.data.frame(output[c("gentime", "inst", "cumul")])
title <- sprintf("%s [mean=%.1f, sd=%.1f]",
                 stringr::str_to_title(output$name),
                 output$mean, sqrt(output$var))

if (output$discrete) {
  # FLAG
  print(paste0("norm:", sum(df$inst)))
  # FLAG-END

  ggplot(df, aes(x = gentime)) +
    geom_bar(aes(y = inst), stat = "identity", fill = "coral") +
    labs(title = title,
         x = "Generation Delay",
         y = "Daily Infectivity Profile")

} else {
  # FLAG
  Dx <- df$gentime[2] - df$gentime[1]
  f <- df$inst
  print(paste0("norm:", Dx * (sum(f[2:(length(f)-1)]) + 0.5 * (f[1] + f[length(f)]))))
  # FLAG-END

  coeff <- max(df$cumul) / max(df$inst)
  ggplot(df, aes(x = gentime)) +
    theme(axis.ticks.y.right = element_line(color = "coral"),
          axis.text.y.right = element_text(color = "coral"),
          axis.title.y.right = element_text(color = "coral")) +
    geom_line(aes(y = cumul / coeff), color = "coral") +
    geom_line(aes(y = inst)) +
    scale_y_continuous(name = "Instantaneous",
                       sec.axis = sec_axis(~. * coeff, name = "Cumulative")) +
    labs(x = "Generation Time", title = title)
}

