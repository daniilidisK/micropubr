# Global defaults and customizability layer ------------------------------
# These let users set package-wide styling once (font size, palette, theme,
# verbosity) instead of repeating the same arguments on every mp_*() call,
# in the spirit of SCpubr's opinionated-but-overridable defaults.

#' Get, set and reset micropubr defaults
#'
#' `micropubr` reads a small set of package-wide defaults that every `mp_*()`
#' plotting function falls back to when the corresponding argument is left at
#' its default (`NULL`). Setting them once with `mp_set_defaults()` gives every
#' subsequent plot a consistent look without repeating arguments.
#'
#' @param font.size Base font size in points used by the default theme.
#' @param palette A palette name (see [mp_palette()]) or an explicit vector of
#'   colours used whenever a function's `palette` argument is `NULL`.
#' @param theme A `ggplot2` theme object used whenever a function's `theme`
#'   argument is `NULL`. If `NULL`, [mp_theme_publication()] is used.
#' @param verbose Logical; print the underlying method and parameters.
#'
#' @return `mp_get_defaults()` returns the current defaults as a named list;
#'   `mp_set_defaults()` and `mp_reset_defaults()` return the (updated)
#'   defaults invisibly.
#' @export
#' @examples
#' old <- mp_get_defaults()
#' mp_set_defaults(font.size = 9, palette = "okabe")
#' mp_get_defaults()
#' mp_reset_defaults()
mp_set_defaults <- function(font.size = NULL, palette = NULL, theme = NULL,
                            verbose = NULL) {
  if (!is.null(font.size)) {
    if (!is.numeric(font.size) || length(font.size) != 1) {
      rlang::abort("`font.size` must be a single number.",
                   class = "micropubr_bad_default")
    }
    options(micropubr.font.size = font.size)
  }
  if (!is.null(palette)) options(micropubr.palette = palette)
  if (!is.null(theme)) {
    if (!inherits(theme, "theme")) {
      rlang::abort("`theme` must be a ggplot2 theme object.",
                   class = "micropubr_bad_default")
    }
    options(micropubr.theme = theme)
  }
  if (!is.null(verbose)) options(micropubr.verbose = isTRUE(verbose))
  invisible(mp_get_defaults())
}

#' @rdname mp_set_defaults
#' @export
mp_get_defaults <- function() {
  list(
    font.size = getOption("micropubr.font.size", 11),
    palette   = getOption("micropubr.palette", NULL),
    theme     = getOption("micropubr.theme", NULL),
    verbose   = getOption("micropubr.verbose", FALSE)
  )
}

#' @rdname mp_set_defaults
#' @export
mp_reset_defaults <- function() {
  options(micropubr.font.size = 11, micropubr.palette = NULL,
          micropubr.theme = NULL, micropubr.verbose = FALSE)
  invisible(mp_get_defaults())
}

#' Resolve the active theme: explicit > global default > publication theme
#' @noRd
.mp_theme <- function(theme = NULL, font.size = 11) {
  theme %||% getOption("micropubr.theme", NULL) %||%
    mp_theme_publication(font.size = font.size)
}

#' Drop NULL elements so `...` can be spliced into a geom safely
#'
#' Lets every plotting function forward `...` to its primary geom via
#' `do.call()` while letting any explicitly-named argument win (the caller
#' removes clashing names before calling this).
#' @noRd
.mp_drop_null <- function(x) {
  x[!vapply(x, is.null, logical(1))]
}

#' Build a geom layer, merging fixed parameters with user `...` overrides
#'
#' `params` are the function's chosen defaults; anything in `dots` with the
#' same name overrides them, and extra `...` names are passed through. This is
#' how `...` reaches the primary geom of each `mp_*()` plot.
#' @noRd
.mp_geom <- function(geom, params = list(), dots = list(), mapping = NULL) {
  args <- utils::modifyList(params, .mp_drop_null(dots))
  if (!is.null(mapping)) args <- c(list(mapping = mapping), args)
  do.call(geom, args)
}
