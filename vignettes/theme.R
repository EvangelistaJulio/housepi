## General Plot Theme for 'housepi'
##

# viggnettes theme
knitr::opts_chunk$set(echo = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 5, fig.fullwidth = TRUE,
                      fig.align = "center", out.height = "auto")
options(knitr.kable.NA = "")

# ggplot theme
my_element_text <- ggplot2::element_text(size = 12, face = "plain")
ggplot2::theme_set(
  ggplot2::theme_gray() +
    ggplot2::theme(
      text = my_element_text,
      axis.text = my_element_text,
      strip.text = my_element_text,
      legend.text = my_element_text,
      legend.title = ggplot2::element_blank()
    )
)

# bayesplot theme
bayesplot::bayesplot_theme_set(ggplot2::theme_get())
bayesplot::color_scheme_set(scheme = "brightblue")
