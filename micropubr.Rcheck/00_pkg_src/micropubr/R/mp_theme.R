# Publication themes ------------------------------------------------------

#' Default publication theme for micropubr figures
#'
#' A clean, minimal `ggplot2` theme with no panel grid clutter, framed axes
#' and balanced typography. It is applied automatically by every `mp_*()`
#' plotting function and can be customised or replaced via the shared `theme`
#' argument.
#'
#' @param font.size Base font size in points. Defaults to the
#'   `micropubr.font.size` option (11).
#' @param base.family Font family. Defaults to "" (device default).
#' @param grid Which major grid lines to draw: `"none"` (default), `"x"`,
#'   `"y"` or `"both"`.
#' @param legend.position Legend position passed to [ggplot2::theme()].
#'
#' @return A `ggplot2` theme object.
#' @export
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   mp_theme_publication()
mp_theme_publication <- function(font.size = getOption("micropubr.font.size", 11),
                                 base.family = "",
                                 grid = c("none", "x", "y", "both"),
                                 legend.position = "right") {
  grid <- match.arg(grid)
  half <- font.size / 2

  thm <- ggplot2::theme_bw(base_size = font.size, base_family = base.family) +
    ggplot2::theme(
      panel.border     = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour = "black",
                                               linewidth = 0.4),
      axis.ticks       = ggplot2::element_line(colour = "black",
                                               linewidth = 0.4),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold",
                                               size = font.size),
      plot.title       = ggplot2::element_text(face = "bold",
                                               size = font.size + 2,
                                               hjust = 0),
      plot.subtitle    = ggplot2::element_text(size = font.size,
                                               colour = "grey30", hjust = 0),
      plot.caption     = ggplot2::element_text(size = half + 1,
                                               colour = "grey40", hjust = 0),
      plot.caption.position = "plot",
      plot.title.position   = "plot",
      legend.position  = legend.position,
      legend.key       = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(face = "bold"),
      axis.text        = ggplot2::element_text(colour = "black")
    )

  grid_lines <- ggplot2::element_line(colour = "grey90", linewidth = 0.3)
  thm <- thm + switch(
    grid,
    none = ggplot2::theme(),
    x    = ggplot2::theme(panel.grid.major.x = grid_lines),
    y    = ggplot2::theme(panel.grid.major.y = grid_lines),
    both = ggplot2::theme(panel.grid.major = grid_lines)
  )
  thm
}

#' Nature-style theme (compact, 8pt sans serif)
#' @inheritParams mp_theme_publication
#' @return A `ggplot2` theme object.
#' @export
mp_theme_nature <- function(font.size = 8, base.family = "") {
  mp_theme_publication(font.size = font.size, base.family = base.family) +
    ggplot2::theme(
      axis.line  = ggplot2::element_line(linewidth = 0.3),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      legend.key.size = ggplot2::unit(0.8, "lines")
    )
}

#' Cell-style theme (boxed panels, slightly larger type)
#' @inheritParams mp_theme_publication
#' @return A `ggplot2` theme object.
#' @export
mp_theme_cell <- function(font.size = 10, base.family = "") {
  ggplot2::theme_bw(base_size = font.size, base_family = base.family) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(colour = "black",
                                               linewidth = 0.6),
      strip.background = ggplot2::element_rect(fill = "grey95",
                                               colour = "black"),
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.caption     = ggplot2::element_text(colour = "grey40", hjust = 0),
      legend.key       = ggplot2::element_blank()
    )
}

#' Dark theme for slides and posters
#' @inheritParams mp_theme_publication
#' @return A `ggplot2` theme object.
#' @export
mp_theme_dark <- function(font.size = 12, base.family = "") {
  bg <- "#1B1B1B"
  fg <- "grey90"
  ggplot2::theme_minimal(base_size = font.size, base_family = base.family) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = bg, colour = NA),
      panel.grid       = ggplot2::element_line(colour = "grey25"),
      text             = ggplot2::element_text(colour = fg),
      axis.text        = ggplot2::element_text(colour = fg),
      axis.line        = ggplot2::element_line(colour = fg),
      legend.key       = ggplot2::element_blank(),
      plot.caption     = ggplot2::element_text(colour = "grey60", hjust = 0)
    )
}
