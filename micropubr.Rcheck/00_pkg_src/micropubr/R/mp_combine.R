# Multi-panel composition --------------------------------------------------

#' Combine micropubr plots into a labelled multi-panel figure
#'
#' A thin convenience wrapper around [patchwork::wrap_plots()] that applies a
#' consistent layout, collects duplicate legends and adds automatic panel tags
#' (A, B, C ...), so figures composed from several `mp_*()` calls look
#' coherent.
#'
#' @param ... `ggplot` objects, or a single named list of them.
#' @param ncol,nrow Grid dimensions (passed to [patchwork::wrap_plots()]).
#' @param guides How to treat legends: `"collect"` (default), `"keep"` or
#'   `"auto"`.
#' @param tag.levels Panel-tag style passed to
#'   [patchwork::plot_annotation()], e.g. `"A"`, `"a"`, `"1"`. Use `NULL` for
#'   no tags.
#' @param title,subtitle,caption Optional figure-level annotations.
#' @param widths,heights Optional relative panel sizes.
#'
#' @return A `patchwork` object.
#' @export
#' @examples
#' \dontrun{
#' a <- mp_alpha(ps, group.by = "Group")
#' o <- mp_ordination(ps, group.by = "Group")
#' mp_combine(a, o, ncol = 2, tag.levels = "A")
#' }
mp_combine <- function(...,
                       ncol = NULL,
                       nrow = NULL,
                       guides = c("collect", "keep", "auto"),
                       tag.levels = "A",
                       title = NULL,
                       subtitle = NULL,
                       caption = NULL,
                       widths = NULL,
                       heights = NULL) {
  .mp_require("patchwork", "Combining plots")
  guides <- match.arg(guides)

  dots <- list(...)
  if (length(dots) == 1 && is.list(dots[[1]]) &&
      !.mp_is_ggplot(dots[[1]])) {
    plots <- dots[[1]]
  } else {
    plots <- dots
  }
  is_plot <- vapply(plots, .mp_is_ggplot, logical(1))
  if (!all(is_plot)) {
    rlang::abort("All inputs to `mp_combine()` must be ggplot objects.",
                 class = "micropubr_bad_input")
  }

  combined <- patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow,
                                    widths = widths, heights = heights,
                                    guides = guides)
  if (!is.null(tag.levels) || !is.null(title) || !is.null(subtitle) ||
      !is.null(caption)) {
    combined <- combined + patchwork::plot_annotation(
      title = title, subtitle = subtitle, caption = caption,
      tag_levels = tag.levels)
  }
  combined
}
