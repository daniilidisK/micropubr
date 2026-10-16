# Figure export ------------------------------------------------------------

# Standard journal column widths (millimetres).
.mp_column_widths <- c(single = 89, onehalf = 120, double = 183)

#' Export a figure at publication resolution and journal column widths
#'
#' Wraps [ggplot2::ggsave()] with sensible defaults for journal submission:
#' device inferred from the file extension, millimetre units, 300 dpi, and
#' optional standard column-width presets so figures are sized correctly for
#' single- or double-column layouts.
#'
#' @param plot A `ggplot` or `patchwork` object, or an interactive `plotly`
#'   widget (as returned with `interactive = TRUE`), which is saved as a
#'   self-contained `.html` page instead.
#' @param filename Output path. For static plots the extension selects the
#'   device (`.pdf`, `.png`, `.tiff`, `.svg`, `.eps`, `.jpg`); for an
#'   interactive widget it must be `.html`. PDF and EPS files are written with
#'   the cairo devices, which embed fonts and support Unicode. Where cairo
#'   cannot be loaded (on macOS it needs XQuartz), PDFs are written with the
#'   Quartz device on macOS and [grDevices::pdf()] elsewhere, and EPS files
#'   with [grDevices::postscript()].
#' @param width,height Dimensions in `units`. If `column` is set, `width` is
#'   taken from the preset; if `height` is `NULL` it defaults to
#'   `width * 0.75`.
#' @param units One of `"mm"` (default), `"cm"` or `"in"`.
#' @param dpi Raster resolution. Common choices: 300 or 600.
#' @param column Column-width preset: `"none"` (default), `"single"` (89 mm),
#'   `"onehalf"` (120 mm) or `"double"` (183 mm).
#' @param format Optional explicit device name overriding the extension.
#' @param ... Further arguments passed to [ggplot2::ggsave()].
#'
#' @return The `filename`, invisibly.
#' @export
#' @examples
#' \dontrun{
#' p <- mp_ordination(ps, group.by = "Group")
#' mp_export(p, "figure1.pdf", column = "single")
#' mp_export(p, "figure1.tiff", column = "double", dpi = 600)
#' }
mp_export <- function(plot,
                      filename,
                      width = NULL,
                      height = NULL,
                      units = c("mm", "cm", "in"),
                      dpi = 300,
                      column = c("none", "single", "onehalf", "double"),
                      format = NULL,
                      ...) {
  # Interactive plotly/htmlwidget objects are saved as a self-contained HTML
  # page (ggsave() cannot rasterise them); everything else goes through ggsave.
  if (methods::is(plot, "plotly") || methods::is(plot, "htmlwidget")) {
    return(.mp_export_widget(plot, filename))
  }
  if (!.mp_is_ggplot(plot) && !methods::is(plot, "patchwork")) {
    rlang::abort("`plot` must be a ggplot, patchwork or plotly object.",
                 class = "micropubr_bad_input")
  }
  units <- match.arg(units)
  column <- match.arg(column)

  if (column != "none") {
    width <- .mp_column_widths[[column]]
    units <- "mm"
  }
  if (is.null(width)) width <- 120
  if (is.null(height)) height <- width * 0.75

  ext <- tolower(tools::file_ext(filename))
  # NULL lets ggsave infer the device from the extension (png/tiff/jpeg/svg,
  # and its base pdf/postscript devices when cairo is unavailable).
  device <- format %||% switch(
    ext,
    pdf = .mp_pdf_device(),
    eps = if (.mp_has_cairo()) grDevices::cairo_ps,
    NULL
  )

  args <- list(filename = filename, plot = plot, width = width,
               height = height, units = units, dpi = dpi, ...)
  if (!is.null(device)) args$device <- device

  do.call(ggplot2::ggsave, args)
  if (getOption("micropubr.verbose", FALSE)) {
    rlang::inform(sprintf("Saved %s (%.0f x %.0f %s, %d dpi).",
                          filename, width, height, units, dpi))
  }
  invisible(filename)
}

#' PDF device: cairo when it loads, else Quartz on macOS, else ggsave's pdf()
#' @noRd
.mp_pdf_device <- function() {
  if (.mp_has_cairo()) return(grDevices::cairo_pdf)
  if (isTRUE(capabilities("aqua"))) {
    return(function(filename, ...) {
      grDevices::quartz(type = "pdf", file = filename, ...)
    })
  }
  NULL
}

#' Can the cairo devices be opened?
#'
#' `capabilities("cairo")` only reports that R was built with cairo. On macOS
#' the cairo module links against XQuartz, and without it `cairo_pdf()` warns
#' and writes nothing, so open a device once per session to find out.
#' @noRd
.mp_has_cairo <- function() {
  if (!is.null(.micropubr_state$cairo)) return(.micropubr_state$cairo)
  ok <- isTRUE(capabilities("cairo"))
  if (ok) {
    old <- grDevices::dev.cur()
    f <- tempfile(fileext = ".pdf")
    ok <- tryCatch({
      suppressWarnings(grDevices::cairo_pdf(f))
      grDevices::dev.cur() != old
    }, error = function(e) FALSE)
    if (ok) {
      grDevices::dev.off()
      if (old > 1) grDevices::dev.set(old)
    }
    unlink(f)
  }
  .micropubr_state$cairo <- ok
  ok
}

#' Save an interactive plotly widget as a self-contained HTML page
#' @noRd
.mp_export_widget <- function(plot, filename) {
  .mp_require("htmlwidgets", "Exporting interactive plots")
  ext <- tolower(tools::file_ext(filename))
  if (!ext %in% c("html", "htm")) {
    rlang::abort(
      c("Interactive plots can only be exported to an HTML file.",
        i = sprintf("Give `filename` an .html extension (got '%s'), or pass a static ggplot.",
                    ext)),
      class = "micropubr_bad_input")
  }
  # A self-contained single file needs pandoc; without it, fall back to an
  # HTML file plus a sidecar `_files/` directory of dependencies.
  selfcontained <- requireNamespace("rmarkdown", quietly = TRUE) &&
    rmarkdown::pandoc_available()
  # saveWidget() requires an absolute path for the accompanying assets.
  htmlwidgets::saveWidget(plot, normalizePath(filename, mustWork = FALSE),
                          selfcontained = selfcontained)
  if (!selfcontained) {
    rlang::warn(c(
      "Saved a non-self-contained widget (pandoc not found).",
      i = "The HTML file needs its accompanying `_files/` directory to render."),
      class = "micropubr_widget_not_selfcontained")
  } else if (getOption("micropubr.verbose", FALSE)) {
    rlang::inform(sprintf("Saved interactive widget to %s.", filename))
  }
  invisible(filename)
}
