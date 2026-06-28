# Interactive (plotly) conversion layer -----------------------------------
# Every ggplot-based mp_*() function accepts `interactive = TRUE`, which routes
# its result through .mp_interactive() to return a `plotly` htmlwidget instead
# of a static ggplot. The conversion is centralised here (mirroring how
# .mp_theme()/.mp_geom() centralise styling) so all functions behave
# consistently: tooltips are enabled, and the subtitle/caption that ggplotly()
# would otherwise drop are re-injected as plotly title/annotation.
#
# Circlize-based figures (mp_chord(), mp_circos(), mp_coverage()) are *not*
# ggplot objects and therefore have no `interactive` argument -- there is no
# ggplot to convert. This is documented on those functions.

#' Convert an mp_*() plot result into an interactive plotly widget
#'
#' Accepts a single `ggplot`/`ggraph`/`ggtree` object, a `patchwork` of them,
#' or a bare list of them (multi-panel output) and returns a `plotly`
#' htmlwidget. Lists and patchworks are assembled into a single
#' [plotly::subplot()].
#'
#' @param plot A ggplot, patchwork, or list of ggplots.
#' @param tooltip Passed to [plotly::ggplotly()]; which aesthetics to show on
#'   hover. Defaults to `"all"`.
#' @param nrows Number of rows when assembling a multi-panel subplot.
#' @return A `plotly` htmlwidget.
#' @noRd
.mp_interactive <- function(plot, tooltip = "all", nrows = NULL) {
  .mp_require("plotly", "Interactive plots (`interactive = TRUE`)")

  # A bare list of panels (e.g. mp_alpha(combine = FALSE)) -> one subplot.
  if (is.list(plot) && !.mp_is_ggplot(plot) &&
        !methods::is(plot, "patchwork")) {
    return(.mp_interactive_panels(plot, tooltip = tooltip, nrows = nrows))
  }
  # patchwork (mp_combine(), multi-measure mp_alpha(combine = TRUE)).
  if (methods::is(plot, "patchwork")) {
    return(.mp_interactive_panels(.mp_patchwork_plots(plot),
                                  tooltip = tooltip, nrows = nrows))
  }
  .mp_one_interactive(plot, tooltip = tooltip)
}

#' Convert one ggplot to plotly, preserving subtitle and caption
#' @noRd
.mp_one_interactive <- function(p, tooltip = "all") {
  if (!.mp_is_ggplot(p)) {
    rlang::abort(
      "`interactive = TRUE` is only supported for ggplot-based plots.",
      class = "micropubr_not_interactive")
  }
  labs <- p$labels
  # ggplotly() emits warnings for geoms/aesthetics it cannot translate
  # (e.g. stat_pvalue brackets); keep the console clean and convert anyway.
  # Some third-party geoms (notably ggtree) can fail outright -- fall back to
  # the static plot with a warning rather than erroring on the user.
  widget <- tryCatch(
    suppressWarnings(plotly::ggplotly(p, tooltip = tooltip)),
    error = function(e) {
      rlang::warn(c(
        "Could not convert this plot to an interactive widget; returning the static ggplot.",
        i = conditionMessage(e)),
        class = "micropubr_interactive_fallback")
      NULL
    })
  if (is.null(widget)) return(p)
  .mp_plotly_relabel(widget, title = labs$title, subtitle = labs$subtitle,
                     caption = labs$caption)
}

#' Assemble several ggplots into a single interactive subplot
#' @noRd
.mp_interactive_panels <- function(plots, tooltip = "all", nrows = NULL) {
  plots <- Filter(.mp_is_ggplot, plots)
  if (length(plots) == 0) {
    rlang::abort("No ggplot panels available to make interactive.",
                 class = "micropubr_not_interactive")
  }
  if (length(plots) == 1) return(.mp_one_interactive(plots[[1]], tooltip))
  widgets <- lapply(plots, function(p) {
    suppressWarnings(plotly::ggplotly(p, tooltip = tooltip))
  })
  if (is.null(nrows)) nrows <- if (length(widgets) > 3) 2 else 1
  plotly::subplot(widgets, nrows = nrows, margin = 0.06,
                  titleX = TRUE, titleY = TRUE, shareX = FALSE,
                  shareY = FALSE)
}

#' Re-attach a ggplot title/subtitle/caption that ggplotly() drops
#'
#' plotly has a single title slot, so the subtitle is folded in as a `<sup>`
#' second line and the caption is added as a bottom-right annotation.
#' @noRd
.mp_plotly_relabel <- function(widget, title = NULL, subtitle = NULL,
                               caption = NULL) {
  title <- .mp_chr1(title)
  subtitle <- .mp_chr1(subtitle)
  caption <- .mp_chr1(caption)

  if (!is.null(title) || !is.null(subtitle)) {
    text <- .mp_html(title %||% "")
    if (!is.null(subtitle)) {
      sub <- sprintf("<sup>%s</sup>", .mp_html(subtitle))
      text <- if (nzchar(text)) paste0(text, "<br>", sub) else sub
    }
    widget <- plotly::layout(widget, title = list(text = text, x = 0.02,
                                                  xanchor = "left"))
  }
  if (!is.null(caption)) {
    widget <- plotly::layout(
      widget,
      margin = list(b = 70),
      annotations = list(list(
        text = .mp_html(caption), showarrow = FALSE,
        xref = "paper", yref = "paper", x = 1, y = -0.16,
        xanchor = "right", yanchor = "top",
        font = list(size = 10, color = "grey40"))))
  }
  widget
}

#' Coerce a label to a single non-empty character string or NULL
#' @noRd
.mp_chr1 <- function(x) {
  if (is.null(x) || length(x) == 0) return(NULL)
  x <- paste(as.character(x), collapse = " ")
  if (!nzchar(trimws(x))) NULL else x
}

#' Minimal HTML escaping + newline -> <br> for plotly text fields
#' @noRd
.mp_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\n", "<br>", x, fixed = TRUE)
}

#' Extract the constituent ggplots from a patchwork, in display order
#'
#' A patchwork is itself the most-recently-added plot, with the earlier plots
#' stored in `$patches$plots`. Strip the patchwork machinery off the top-level
#' plot so each panel converts as a plain ggplot.
#' @noRd
.mp_patchwork_plots <- function(pw) {
  earlier <- pw$patches$plots %||% list()
  top <- pw
  top$patches <- NULL
  class(top) <- setdiff(class(top), c("patchwork"))
  panels <- c(earlier, list(top))
  Filter(.mp_is_ggplot, panels)
}
