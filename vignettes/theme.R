## General Plot Theme for 'housepi'
##

# vignettes theme
knitr::opts_chunk$set(echo = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 5, fig.fullwidth = TRUE,
                      fig.align = "center", out.height = "auto")
options(knitr.kable.NA = "")

# ggplot theme
my_element_text <- ggplot2::element_text(size = 12, face = "plain")
ggplot2::theme_set(
  ggplot2::theme_gray() +
    ggplot2::theme(text = my_element_text,
                   axis.text = my_element_text,
                   strip.text = my_element_text,
                   legend.text = my_element_text,
                   legend.title = ggplot2::element_blank())
)
# mimic element_geom for older ggplot2 versions
ggplot2::update_geom_defaults("line", list(color = "coral", size = 1))
ggplot2::update_geom_defaults("point", list(color = "coral", size = 3))
ggplot2::update_geom_defaults("vline", list(color = "coral", size = 1))
ggplot2::update_geom_defaults("hline", list(color = "coral", size = 1))
ggplot2::update_geom_defaults("abline", list(color = "coral", size = 1))
ggplot2::update_geom_defaults("segment", list(color = "coral", size = 1))
ggplot2::update_geom_defaults("errorbar", list(color = "coral", size = 1))

# bayesplot theme
bayesplot::bayesplot_theme_set(ggplot2::theme_get())
bayesplot::color_scheme_set(scheme = "brightblue")

# session info
show_session_info <- function(sampling = TRUE) {
  output <- ""

  if (sampling) output <- paste0(
    cli::symbol$warning,
    " The number of steps were reduced for demonstration.",
    " Increase them in realistic applications.<br>"
  )

  output <- paste0(
    output,
    cli::symbol$info,
    paste0(utils::capture.output(sessionInfo())[1:3], collapse = ". "),
    "."
  )

  cat(sprintf(
"<div style='
  border:1px solid #ccc;
  padding:0.6em 1em;
  border-radius:6px;
  background:#f9f9f9;
  font-size:0.85em;
  color:#444;
  text-align:left;
'>
%s
</div>",
    output
  ))
}
