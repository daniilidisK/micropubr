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

#' Convert a finished figure to dark mode
#'
#' Recolours a finished figure for a dark background while keeping its layout.
#' Adding [mp_theme_dark()] to a finished plot replaces its whole theme, which
#' undoes layout that the plotting function set (rotated or wrapped axis
#' labels, bold facet strips, the absence of grid lines), and it cannot reach
#' dark text and lines that the layers draw themselves: feature labels,
#' significance brackets, tree branches. `mp_darken()` changes only colours:
#' the backgrounds and the theme's text and line colours, plus any near-black
#' colour a layer sets that is not mapped to data.
#'
#' @param plot A `ggplot` or `patchwork` figure.
#' @param background Background colour.
#' @param foreground Colour for text, axes and formerly dark marks.
#' @return The recoloured figure; `plot` itself is not modified.
#' @export
#' @seealso [mp_theme_dark()] to build a dark figure from the start.
#' @examples
#' df <- data.frame(feature = paste0("f", 1:30), log2FC = rnorm(30, sd = 2),
#'                  padj = runif(30)^3)
#' mp_darken(mp_diff(df, plot.type = "volcano"))
mp_darken <- function(plot, background = "#1B1B1B", foreground = "grey90") {
  if (inherits(plot, "patchwork")) {
    for (i in seq_along(plot)) {
      plot[[i]] <- mp_darken(plot[[i]], background, foreground)
    }
    return(plot & .mp_dark_overlay(NULL, background, foreground))
  }
  if (!.mp_is_ggplot(plot)) {
    rlang::abort("`plot` must be a ggplot or patchwork figure.",
                 class = "micropubr_bad_input")
  }
  plot$layers <- lapply(plot$layers, .mp_darken_layer, plot = plot,
                        foreground = foreground)
  plot + .mp_dark_overlay(plot, background, foreground)
}

#' Theme colours for dark mode, touching only elements the plot draws
#'
#' Adding an element to a theme where it is blank would switch it on (a
#' figure without grid lines would gain them), so line and strip elements are
#' recoloured only when the plot's own theme draws them. With `plot = NULL`
#' only backgrounds and text are set, which is safe on any figure.
#' @noRd
.mp_dark_overlay <- function(plot, background, foreground) {
  th <- ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = background, colour = NA),
    panel.background = ggplot2::element_rect(fill = background, colour = NA),
    legend.background = ggplot2::element_rect(fill = background, colour = NA),
    legend.key = ggplot2::element_rect(fill = background, colour = NA),
    text = ggplot2::element_text(colour = foreground),
    axis.text = ggplot2::element_text(colour = foreground),
    axis.title = ggplot2::element_text(colour = foreground),
    strip.text = ggplot2::element_text(colour = foreground),
    legend.text = ggplot2::element_text(colour = foreground),
    legend.title = ggplot2::element_text(colour = foreground),
    plot.title = ggplot2::element_text(colour = foreground),
    plot.subtitle = ggplot2::element_text(colour = "grey75"),
    plot.caption = ggplot2::element_text(colour = "grey65"))
  # ggplot2 >= 4.0 draws default text and lines in the theme's "ink".
  if (exists("element_geom", asNamespace("ggplot2"), inherits = FALSE)) {
    element_geom <- get("element_geom", envir = asNamespace("ggplot2"))
    th <- th + ggplot2::theme(geom = element_geom(ink = foreground,
                                                  paper = background))
  }
  if (is.null(plot)) return(th)

  full <- tryCatch(ggplot2::theme_get() + plot$theme,
                   error = function(e) ggplot2::theme_get())
  drawn <- function(...) {
    any(vapply(c(...), function(el) {
      e <- tryCatch(ggplot2::calc_element(el, full), error = function(e) NULL)
      !is.null(e) && !inherits(e, "element_blank")
    }, logical(1)))
  }
  if (drawn("axis.line.x.bottom", "axis.line.y.left")) {
    th <- th + ggplot2::theme(axis.line = ggplot2::element_line(colour = foreground))
  }
  if (drawn("axis.ticks.x.bottom", "axis.ticks.y.left")) {
    th <- th + ggplot2::theme(axis.ticks = ggplot2::element_line(colour = foreground))
  }
  if (drawn("panel.grid.major.x", "panel.grid.major.y")) {
    th <- th + ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey28"))
  }
  if (drawn("panel.grid.minor.x", "panel.grid.minor.y")) {
    th <- th + ggplot2::theme(panel.grid.minor = ggplot2::element_line(colour = "grey22"))
  }
  if (drawn("panel.border")) {
    th <- th + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "grey40", fill = NA))
  }
  if (drawn("strip.background.x", "strip.background.y")) {
    th <- th + ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey20", colour = NA))
  }
  th
}

#' A copy of `layer` whose fixed near-black colour is made light
#'
#' Colours mapped to data are left alone (unless a fixed colour overrides the
#' mapping), and the layer is cloned rather than
#' modified: ggplot layers are shared, so changing one in place would also
#' recolour the light figure it came from.
#' @noRd
.mp_darken_layer <- function(layer, plot, foreground) {
  # A fixed colour wins over a mapped one, so check it first; only without
  # one does a mapping (data-driven colour) make the layer off-limits.
  current <- layer$aes_params$colour %||% layer$aes_params$color
  if (is.null(current)) {
    mapped <- names(layer$mapping)
    if (isTRUE(layer$inherit.aes)) mapped <- c(mapped, names(plot$mapping))
    if (any(c("colour", "color") %in% mapped)) return(layer)
    current <- layer$geom$default_aes$colour
    # A default computed from the theme follows the overlay's ink instead.
    if (!is.character(current)) return(layer)
  }
  if (!.mp_is_dark(current)) return(layer)
  params <- layer$aes_params
  params$color <- NULL
  params$colour <- foreground
  ggplot2::ggproto(NULL, layer, aes_params = params)
}

#' Is a single colour near-black (relative luminance below 0.35)?
#' @noRd
.mp_is_dark <- function(col) {
  if (length(col) != 1 || is.na(col)) return(FALSE)
  rgb <- tryCatch(grDevices::col2rgb(col), error = function(e) NULL)
  if (is.null(rgb)) return(FALSE)
  sum(c(0.2126, 0.7152, 0.0722) * rgb[, 1]) / 255 < 0.35
}
