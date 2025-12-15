## General Plot Theme for housepi
##

# vignettes theme
knitr::opts_chunk$set(echo = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 5, fig.fullwidth = TRUE,
                      fig.align = "center", out.height = "auto")
options(knitr.kable.NA = "")

# ggplot theme
my_elem_text <- ggplot2::element_text(size = 12, face = "plain")
ggplot2::theme_set(
  ggplot2::theme_gray() +
    ggplot2::theme(text = my_elem_text,
                   axis.text = my_elem_text,
                   strip.text = my_elem_text,
                   legend.text = my_elem_text,
                   legend.title = ggplot2::element_blank())
)
# mimic element_geom for older ggplot2 versions
ggplot2::update_geom_defaults("line", list(color = "#FF7F50", size = 1))
ggplot2::update_geom_defaults("point", list(color = "#FF7F50", size = 3))
ggplot2::update_geom_defaults("vline", list(color = "#FF7F50", size = 1))
ggplot2::update_geom_defaults("hline", list(color = "#FF7F50", size = 1))
ggplot2::update_geom_defaults("abline", list(color = "#FF7F50", size = 1))
ggplot2::update_geom_defaults("segment", list(color = "#FF7F50", size = 1))
ggplot2::update_geom_defaults("errorbar", list(color = "#FF7F50", size = 1,
                                               width = 0.1))
# bayesplot theme
bayesplot::bayesplot_theme_set(ggplot2::theme_get())
bayesplot::color_scheme_set(c("#6ec5ff", "#3db1ff",
                              "#0d8fff", "#0070dd",
                              "#0055aa", "#003d77"))

# flextable theme
flextable::set_flextable_defaults(
  font.family = "Roboto", font.size = 13, font.color = "#111111",
  background.color = "white",
  theme_fun = "theme_vanilla"
)

#' @title Show session info
#' @param sampling `logical(1)`
#' @keywords internal
#' @noRd
show_session_info <- function(sampling = TRUE) {
  if (sampling) output <- paste0(
    cli::symbol$warning,
    " Number of steps were reduced for demonstration.",
    " Increase them in realistic applications.<br>"
  )
  else output <- ""

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
