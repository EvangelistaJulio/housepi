## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' @title Set a warm default theme for housepi visualizations
#' @param base_size `numeric(1)` Base font size, in points.
#' @param base_family `character(1)` Base font family, `""` for the device's.
#' @param transparent `logical(1)`
#' Whether to leave the plot background transparent.
#' @param colour_strips `logical(1)`
#' Whether to colour the facet strips and centre their labels,
#' as in the housepi vignettes.
#' Otherwise strips are white and labels align to the bottom of the strip,
#' next to their panel.
#' @description
#' Set the housepi themes of `ggplot2`, `bayesplot` and `flextable`,
#' for each package installed.
#' @return `NULL`, invisibly.
#' The function sets the themes.
#' @examples
#' \donttest{
#' # apply housepi defaults to ggplot2, bayesplot and flextable
#' set_housepi_theme()
#' }
#' @export
set_housepi_theme <- function(base_size = 11,
                              base_family = "",
                              transparent = FALSE,
                              colour_strips = FALSE) {
  # ggplot theme
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    my_elem_text <- ggplot2::element_text(size = base_size,
                                          color = "black",
                                          face = "plain")
    my_elem_tag <- ggplot2::element_text(size = base_size + 1,
                                         color = "black",
                                         face = "bold")
    my_panel <- ggplot2::element_rect(fill = "#FEF9EC", color = NA)
    my_strip <- ggplot2::element_rect(
      fill = if (isTRUE(colour_strips)) "#FBE09B" else "white", color = NA
    )
    my_strip_text <- ggplot2::element_text(
      size = base_size, color = "black", face = "plain",
      vjust = if (isTRUE(colour_strips)) NULL else 0
    )
    my_grid <- ggplot2::element_line(color = "#FDFEFE", linewidth = 1)
    my_plot <- ggplot2::element_rect(
      color = NA, fill = if (isTRUE(transparent)) "transparent" else "white"
    )
    ggplot2::theme_set(
      ggplot2::theme_minimal(base_size = base_size,
                             base_family = base_family) +
        ggplot2::theme(text = my_elem_text,
                       axis.text = my_elem_text,
                       strip.text = my_strip_text,
                       plot.tag = my_elem_tag,
                       panel.background = my_panel,
                       strip.background = my_strip,
                       plot.background = my_plot,
                       legend.position = "bottom",
                       legend.text = my_elem_text,
                       panel.grid.major = my_grid,
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.title = ggplot2::element_blank())
    )

    # mimic element_geom for older ggplot2 versions
    my_text_feats <- list(size = (base_size - 1) / ggplot2::.pt)
    my_point_feats <- list(color = "#FF7F50", size = 2)
    my_line_feats <- list(color = "#FF7F50", size = 1)
    my_errbar_feats <- c(my_line_feats, width = 0.1)
    # labels on a borderless box of the panel colour, hiding the grid lines
    my_label_feats <- c(my_text_feats, fill = my_panel$fill, linewidth = 0)
    ggplot2::update_geom_defaults("text", my_text_feats)
    ggplot2::update_geom_defaults("label", my_label_feats)
    ggplot2::update_geom_defaults("line", my_line_feats)
    ggplot2::update_geom_defaults("point", my_point_feats)
    ggplot2::update_geom_defaults("vline", my_line_feats)
    ggplot2::update_geom_defaults("hline", my_line_feats)
    ggplot2::update_geom_defaults("abline", my_line_feats)
    ggplot2::update_geom_defaults("segment", my_line_feats)
    ggplot2::update_geom_defaults("errorbar", my_errbar_feats)

    # bayesplot theme
    if (requireNamespace("bayesplot", quietly = TRUE)) {
      bayesplot::bayesplot_theme_set(ggplot2::theme_get())
      bayesplot::color_scheme_set(c("#6ec5ff", "#3db1ff", "#0d8fff",
                                    "#0070dd", "#0055aa", "#003d77"))
    }
  }

  # flextable theme
  if (requireNamespace("flextable", quietly = TRUE))
    flextable::set_flextable_defaults(font.family = base_family,
                                      font.size = base_size,
                                      font.color = "black",
                                      background.color = "white",
                                      padding.left = 2,
                                      padding.right = 2,
                                      theme_fun = "theme_vanilla")
  invisible()
}
