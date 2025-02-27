library(dplyr)
library(ggplot2)

# house level tibble
house_df <- housepi::survey %>% 
  group_by(house_id) %>% 
  add_count() %>% 
  summarise(size = unique(n), cases = sum(status, na.rm = TRUE))
#  summarise(size = unique(n), cases = sum(status), surface = unique(HH_surface))

# binomial mean error function
err_fun <- function(x) { 
  sqrt(mean(x) * (1 - mean(x) / length(x))) / sqrt(length(x)) 
}

# secondary attack rate tibble
sar_df <- house_df %>% 
  filter(size > 1, cases > 0) %>%
  group_by(size) %>% 
  add_count() %>% 
  filter(n > 2) %>%
  summarise(sar = (mean(cases - 1)) / (unique(size) - 1), 
            sar_sd = (err_fun(cases - 1)) / (unique(size) - 1),
            counts = unique(n))

# plot theme
mytheme <- theme_bw() + 
  theme(text = element_text(size = 12, color = "black"),
        axis.text = element_text(size = 12, color = "black"))

# plot SAR
col <- "#F8766D"
ggplot(sar_df, aes(x = size, y = sar)) +
  # horizontal line with mean SAR
  geom_hline(aes(yintercept = mean(sar)), linetype="dashed") +
  # labels with household counts
  geom_label(aes(x = size + 0.35, y = 0.95 * sar, label = sprintf("[%s]", counts)),
             label.size = NA, # avoid showing padding contour
             size = mytheme$text$size / 3, family = mytheme$text$family) +
  # error bars
  geom_errorbar(aes(ymin = sar - sar_sd, ymax = sar + sar_sd), 
                width = .1, position = position_dodge(), color = col) +
  # scatter plot
  geom_point(size = 2, color = col) +
  # labels and ticks
  mytheme + 
  scale_x_continuous("Household Size", 
                     breaks = seq(min(sar_df$size), max(sar_df$size)), 
                     minor_breaks = NULL) +
  scale_y_continuous("Secondary Attack Rate", breaks = waiver())
# save plot
#ggsave("test.pdf", width = 11, height = 8, units = "cm", dpi = 300)

  