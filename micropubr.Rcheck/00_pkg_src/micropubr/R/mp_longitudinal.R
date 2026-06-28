# Longitudinal / paired-sample analysis ----------------------------------

#' Plot longitudinal and paired-sample microbiome data
#'
#' Visualises time-course or before/after designs: per-subject diversity
#' trajectories (spaghetti), ordination trajectories connecting paired
#' timepoints, alluvial flows of taxonomic shifts, and small-multiple
#' "horizon" area bands of taxon dynamics.
#'
#' @param object A [phyloseq::phyloseq] object.
#' @param time.var Metadata column giving the timepoint (numeric or ordered).
#' @param subject.var Metadata column identifying the subject / individual.
#' @param plot.type One of `"spaghetti"`, `"trajectory"`, `"alluvial"` or
#'   `"horizon"`.
#' @param group.by Optional metadata column for colouring (e.g. treatment
#'   arm).
#' @param measure Alpha-diversity index for `"spaghetti"` plots.
#' @param tax.level Taxonomic rank for `"alluvial"` / `"horizon"` plots.
#' @param top.n Number of taxa to display in `"alluvial"` / `"horizon"`.
#' @param distance Distance metric for `"trajectory"` ordination.
#' @param smooth Add a per-group LOESS smooth to `"spaghetti"` plots.
#' @param palette Palette name or colour vector.
#' @param theme A `ggplot2` theme; defaults to [mp_theme_publication()].
#' @param plot.title,plot.subtitle,plot.caption Optional labels.
#' @param font.size Base font size.
#' @param verbose Print the steps used.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' mp_longitudinal(ps, time.var = "Week", subject.var = "SubjectID",
#'                 group.by = "Arm", measure = "Shannon")
#' }
mp_longitudinal <- function(object,
                            time.var,
                            subject.var,
                            plot.type = c("spaghetti", "trajectory",
                                          "alluvial", "horizon"),
                            group.by = NULL,
                            measure = "Shannon",
                            tax.level = "Phylum",
                            top.n = 8,
                            distance = "bray",
                            smooth = TRUE,
                            palette = NULL,
                            theme = NULL,
                            plot.title = NULL,
                            plot.subtitle = NULL,
                            plot.caption = NULL,
                            font.size = getOption("micropubr.font.size", 11),
                            verbose = getOption("micropubr.verbose", FALSE)) {
  .mp_check_phyloseq(object)
  plot.type <- match.arg(plot.type)
  time.var <- .mp_check_var(object, time.var, "time.var", required = TRUE)
  subject.var <- .mp_check_var(object, subject.var, "subject.var",
                               required = TRUE)
  group.by <- .mp_check_var(object, group.by, "group.by")
  theme <- theme %||% mp_theme_publication(font.size = font.size)

  if (verbose) {
    rlang::inform(sprintf("mp_longitudinal(): %s over '%s' by '%s'.",
                          plot.type, time.var, subject.var))
  }

  p <- switch(
    plot.type,
    spaghetti  = .mp_long_spaghetti(object, time.var, subject.var, group.by,
                                    measure, smooth, palette),
    trajectory = .mp_long_trajectory(object, time.var, subject.var, group.by,
                                     distance, palette),
    alluvial   = .mp_long_alluvial(object, time.var, group.by, tax.level,
                                   top.n, palette),
    horizon    = .mp_long_horizon(object, time.var, tax.level, top.n, palette)
  )

  p + theme +
    ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                  caption = plot.caption)
}

#' Coerce a time variable to numeric where possible
#' @return list(values = numeric/integer codes, labels, numeric = logical)
#' @noRd
.mp_time_axis <- function(x) {
  if (is.numeric(x)) return(list(values = x, numeric = TRUE))
  num <- suppressWarnings(as.numeric(as.character(x)))
  if (!any(is.na(num))) return(list(values = num, numeric = TRUE))
  f <- if (is.ordered(x)) x else factor(x)
  list(values = f, numeric = FALSE)
}

#' Spaghetti plot of alpha diversity over time
#' @noRd
.mp_long_spaghetti <- function(object, time.var, subject.var, group.by,
                               measure, smooth, palette) {
  alpha <- .mp_estimate_alpha(object, measure)
  meta <- .mp_sample_data(object)
  df <- cbind(meta[rownames(alpha), , drop = FALSE], alpha)
  ta <- .mp_time_axis(df[[time.var]])
  df$.time <- ta$values

  aes_main <- ggplot2::aes(x = .data$.time, y = .data[[measure]],
                           group = .data[[subject.var]])
  if (!is.null(group.by)) {
    df[[group.by]] <- as.factor(df[[group.by]])
    aes_main$colour <- ggplot2::aes(colour = .data[[group.by]])$colour
  }

  p <- ggplot2::ggplot(df, aes_main) +
    ggplot2::geom_line(alpha = 0.4, linewidth = 0.5) +
    ggplot2::geom_point(size = 1.6, alpha = 0.7) +
    ggplot2::labs(x = time.var, y = measure)

  if (smooth) {
    smooth_aes <- if (!is.null(group.by)) {
      ggplot2::aes(x = .data$.time, y = .data[[measure]],
                   group = .data[[group.by]], colour = .data[[group.by]])
    } else {
      ggplot2::aes(x = .data$.time, y = .data[[measure]], group = 1)
    }
    p <- p + ggplot2::geom_smooth(
      mapping = smooth_aes, method = "loess", formula = y ~ x,
      se = TRUE, linewidth = 1, alpha = 0.15, inherit.aes = FALSE)
  }

  if (!is.null(group.by)) {
    cols <- .mp_resolve_palette(palette, nlevels(df[[group.by]]),
                                type = "qualitative")
    p <- p + ggplot2::scale_colour_manual(values = cols, name = group.by)
  }
  p
}

#' Ordination trajectory connecting paired timepoints per subject
#' @noRd
.mp_long_trajectory <- function(object, time.var, subject.var, group.by,
                                distance, palette) {
  ord <- .mp_ord_pcoa(object, distance, c(1, 2))
  meta <- .mp_sample_data(object)
  df <- cbind(ord$scores, meta[rownames(ord$scores), , drop = FALSE])
  xvar <- names(ord$scores)[1]
  yvar <- names(ord$scores)[2]
  ta <- .mp_time_axis(df[[time.var]])
  df$.time <- if (ta$numeric) ta$values else as.integer(ta$values)
  df <- df[order(df[[subject.var]], df$.time), , drop = FALSE]

  aes_main <- ggplot2::aes(x = .data[[xvar]], y = .data[[yvar]],
                           group = .data[[subject.var]])
  if (!is.null(group.by)) {
    df[[group.by]] <- as.factor(df[[group.by]])
    aes_main$colour <- ggplot2::aes(colour = .data[[group.by]])$colour
  }

  p <- ggplot2::ggplot(df, aes_main) +
    ggplot2::geom_path(arrow = ggplot2::arrow(
      length = ggplot2::unit(0.015, "npc"), type = "closed"),
      alpha = 0.6, linewidth = 0.5) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::labs(x = ord$xlab, y = ord$ylab,
                  caption = sprintf("PCoA trajectories (%s distance)", distance))

  if (!is.null(group.by)) {
    cols <- .mp_resolve_palette(palette, nlevels(df[[group.by]]),
                                type = "qualitative")
    p <- p + ggplot2::scale_colour_manual(values = cols, name = group.by)
  }
  p
}

#' Mean taxon composition per timepoint, long format (helper)
#' @noRd
.mp_long_composition <- function(object, time.var, group.by, tax.level,
                                 top.n) {
  df <- .mp_melt_taxa(object, tax.level = tax.level, top.n = top.n,
                      relative = TRUE)
  keep <- unique(c(".taxon", time.var, group.by))
  keep <- intersect(keep, names(df))
  f <- stats::as.formula(
    paste("Abundance ~", paste(sprintf("`%s`", keep), collapse = " + ")))
  stats::aggregate(f, data = df, FUN = mean, na.rm = TRUE)
}

#' Alluvial diagram of taxonomic shifts between timepoints
#' @noRd
.mp_long_alluvial <- function(object, time.var, group.by, tax.level, top.n,
                              palette) {
  .mp_require("ggalluvial", "Alluvial diagrams")
  agg <- .mp_long_composition(object, time.var, group.by, tax.level, top.n)
  agg[[time.var]] <- factor(agg[[time.var]])
  cols <- .mp_taxa_colors(palette, agg$.taxon, "Other")

  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(x = .data[[time.var]], y = .data$Abundance,
                 alluvium = .data$.taxon, stratum = .data$.taxon,
                 fill = .data$.taxon)) +
    ggalluvial::geom_flow(alpha = 0.6) +
    ggalluvial::geom_stratum(alpha = 0.95, colour = "grey30",
                             linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = cols, name = tax.level) +
    ggplot2::labs(x = time.var, y = "Mean relative abundance")

  if (!is.null(group.by)) {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", sprintf("`%s`", group.by))))
  }
  p
}

#' Small-multiple "horizon"-style area bands of taxon dynamics
#' @noRd
.mp_long_horizon <- function(object, time.var, tax.level, top.n, palette) {
  agg <- .mp_long_composition(object, time.var, group.by = NULL, tax.level,
                              top.n)
  ta <- .mp_time_axis(agg[[time.var]])
  agg$.time <- ta$values
  cols <- .mp_taxa_colors(palette, agg$.taxon, "Other")

  ggplot2::ggplot(agg, ggplot2::aes(x = .data$.time, y = .data$Abundance,
                                    fill = .data$.taxon)) +
    ggplot2::geom_area(alpha = 0.85, colour = NA) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$.taxon),
                        scales = "free_y", switch = "y") +
    ggplot2::labs(x = time.var, y = "Relative abundance") +
    ggplot2::theme(
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
      strip.placement = "outside",
      axis.text.y = ggplot2::element_text(size = 6),
      panel.spacing.y = ggplot2::unit(2, "pt"))
}
