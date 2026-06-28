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
#' @param plot A `ggplot` or `patchwork` object.
#' @param filename Output path. The extension selects the device (`.pdf`,
#'   `.png`, `.tiff`, `.svg`, `.eps`, `.jpg`).
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
  if (!.mp_is_ggplot(plot) && !methods::is(plot, "patchwork")) {
    rlang::abort("`plot` must be a ggplot or patchwork object.",
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
  device <- format %||% switch(
    ext,
    pdf = grDevices::cairo_pdf,
    eps = grDevices::cairo_ps,
    NULL  # let ggsave infer png/tiff/jpeg/svg from the extension
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
