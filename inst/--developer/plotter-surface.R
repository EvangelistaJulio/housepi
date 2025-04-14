######################
### SAR by SURFACE ###
######################

library(dplyr)
library(ggplot2)

dataset <- A_H1N1_copanflu_B1_T80 %>%
  rename(house_id = house_index, surface = HH_surface)

breaks <- c(min(dataset$surface), 40, 60, 80, 100, 120, max(dataset$surface))

house_df <- dataset %>%
  group_by(house_id) %>%
  add_count() %>%
  summarise(size = unique(n),
            cases = sum(status),
            surface = unique(surface)) %>%
  mutate(bin = cut(surface, breaks = breaks))

# binomial mean error function
err_fun <- function(x) {
  sqrt(mean(x) * (1 - mean(x) / length(x))) / sqrt(length(x))
}

# secondary attack rate tibble
sar_df <- house_df %>%
  filter(size > 1, cases > 0) %>%
  group_by(size, bin) %>%
  add_count() %>%
  filter(n > 4) %>%
  summarise(sar = (mean(cases - 1)) / (unique(size) - 1), 
            sar_sd = (err_fun(cases - 1)) / (unique(size) - 1),
            counts = unique(n)) %>%
  filter(sar_sd > 0) %>%
  mutate_at(vars(size), factor)

# plot theme
mytheme <- theme_bw() +
  theme(text = element_text(size = 12, color = "black"),
        axis.text.x = element_text(angle = 30, vjust = 0.5),
        axis.text = element_text(size = 12, color = "black"))

# plot SAR
ggplot(sar_df, aes(x = bin, y = sar, color = size, group = 1, fill = size)) +
  geom_ribbon(aes(ymin = sar - sar_sd, ymax = sar + sar_sd), alpha = 0.2, size = NA) +
  # labels with household counts
  geom_label(aes(x = bin, y = sar - 0.025, label = sprintf("[%s]", counts)),
             label.size = NA, # avoid showing padding contour
             size = mytheme$text$size / 3, family = mytheme$text$family,
             fill = NA, show.legend = FALSE) +
  geom_point() +
  facet_wrap(~ size, ncol = 2) +
  mytheme +
  labs(x = "Surface [m^2]", y = "Secondary Attack Rate")
