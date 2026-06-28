# Taxonomic composition ---------------------------------------------------

#' Plot taxonomic composition
#'
#' Visualises community composition at a chosen taxonomic rank as stacked
#' bars, smoothed area charts, heatmaps or bubble plots. The top `top.n` taxa
#' are shown individually and the remainder are collapsed into an "Other"
#' category with a muted colour.
#'
#' @param object A [phyloseq::phyloseq] object.
#' @param tax.level Taxonomic rank to agglomerate to (e.g. `"Phylum"`,
#'   `"Genus"`). Must be one of [phyloseq::rank_names()].
#' @param plot.type One of `"bar"`, `"area"`, `"heatmap"`, `"bubble"` or
#'   `"treemap"` (the latter requires the `treemapify` package).
#' @param group.by Optional metadata column to facet samples by.
#' @param top.n Number of most-abundant taxa to display individually.
#' @param relative Show relative abundance (`TRUE`) or raw counts.
#' @param order.by How to order samples on the x-axis: a metadata column
#'   name, `"abundance"` (by the most dominant taxon, the default), or
#'   `"cluster"` (hierarchical clustering on Bray-Curtis distances). Ignored
#'   when `cluster.samples = TRUE`.
#' @param cluster.samples Order samples by hierarchical clustering; shorthand
#'   for `order.by = "cluster"`.
#' @param other.label Label for the collapsed low-abundance category.
#' @param palette Palette name or colour vector; defaults to a taxon palette.
#' @param theme A `ggplot2` theme; defaults to the active theme (see
#'   [mp_set_defaults()], otherwise [mp_theme_publication()]).
#' @param ... Additional arguments passed to the primary geom (e.g.
#'   [ggplot2::geom_col()] for bars).
#' @param plot.title,plot.subtitle,plot.caption Optional labels.
#' @param font.size Base font size.
#' @param return.data If `TRUE`, return the aggregated taxon-by-sample tidy
#'   data frame instead of a plot.
#' @param interactive If `TRUE`, return an interactive `plotly` widget instead
#'   of a static `ggplot` (requires `plotly`). See [mp_alpha()].
#' @param verbose Print the methods used.
#'
#' @return A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
#'   data frame when `return.data = TRUE`.
#' @export
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' mp_taxa(GlobalPatterns, tax.level = "Phylum", group.by = "SampleType")
#' }
mp_taxa <- function(object,
                    tax.level = "Phylum",
                    plot.type = c("bar", "area", "heatmap", "bubble",
                                  "treemap"),
                    group.by = NULL,
                    top.n = 10,
                    relative = TRUE,
                    order.by = "abundance",
                    cluster.samples = FALSE,
                    other.label = "Other",
                    palette = NULL,
                    theme = NULL,
                    ...,
                    plot.title = NULL,
                    plot.subtitle = NULL,
                    plot.caption = NULL,
                    font.size = getOption("micropubr.font.size", 11),
                    interactive = FALSE,
                    return.data = FALSE,
                    verbose = getOption("micropubr.verbose", FALSE)) {
  .mp_check_phyloseq(object)
  plot.type <- match.arg(plot.type)
  group.by <- .mp_check_var(object, group.by, "group.by")
  theme <- .mp_theme(theme, font.size)
  if (cluster.samples) order.by <- "cluster"

  if (verbose) {
    rlang::inform(sprintf(
      "mp_taxa(): %s plot at rank '%s', top %d taxa%s.",
      plot.type, tax.level, top.n,
      if (!is.null(group.by)) sprintf(", faceted by '%s'", group.by) else ""))
  }

  df <- .mp_melt_taxa(object, tax.level = tax.level, top.n = top.n,
                      relative = relative, other.label = other.label)

  keep_cols <- unique(c(".taxon", "Sample", group.by,
                        if (order.by %in% phyloseq::sample_variables(object))
                          order.by))
  agg <- .mp_taxa_aggregate(df, keep_cols)

  # Order samples on the x-axis.
  sample_order <- .mp_order_samples(object, agg, order.by, top.n,
                                    relative, other.label)
  agg$Sample <- factor(agg$Sample, levels = sample_order)

  if (isTRUE(return.data)) return(agg)

  ylab <- if (relative) "Relative abundance" else "Abundance"
  caption <- plot.caption %||%
    sprintf("Top %d %s shown; remainder grouped as \"%s\".",
            top.n, tax.level, other.label)

  dots <- list(...)
  p <- switch(
    plot.type,
    bar     = .mp_taxa_bar(agg, palette, ylab, other.label, area = FALSE,
                           dots = dots),
    area    = .mp_taxa_bar(agg, palette, ylab, other.label, area = TRUE,
                           dots = dots),
    heatmap = .mp_taxa_heatmap(agg, palette, ylab),
    bubble  = .mp_taxa_bubble(agg, palette, ylab),
    treemap = .mp_taxa_treemap(agg, palette, tax.level)
  )

  if (!is.null(group.by) && plot.type %in% c("bar", "area", "bubble",
                                             "heatmap")) {
    p <- p + ggplot2::facet_grid(
      stats::as.formula(paste("~", sprintf("`%s`", group.by))),
      scales = "free_x", space = "free_x"
    )
  }

  p <- p +
    ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                  caption = caption) +
    theme +
    .mp_rotate_x(levels(agg$Sample), font.size + 1, angle = 90)
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

#' Aggregate melted taxa to one row per taxon/sample/group combination
#' @noRd
.mp_taxa_aggregate <- function(df, keep_cols) {
  keep_cols <- intersect(keep_cols, names(df))
  f <- stats::as.formula(
    paste("Abundance ~", paste(sprintf("`%s`", keep_cols), collapse = " + "))
  )
  stats::aggregate(f, data = df, FUN = sum, na.rm = TRUE)
}

#' Determine sample ordering for the x-axis
#' @noRd
.mp_order_samples <- function(object, agg, order.by, top.n, relative,
                              other.label) {
  samples <- unique(as.character(agg$Sample))
  if (identical(order.by, "cluster")) {
    d <- .mp_compute_distance(object, "bray")
    hc <- stats::hclust(d, method = "average")
    return(hc$labels[hc$order])
  }
  if (identical(order.by, "abundance") || is.null(order.by)) {
    # Order by abundance of the single most dominant taxon overall.
    tot <- tapply(agg$Abundance, agg$.taxon, sum, na.rm = TRUE)
    top_taxon <- names(sort(tot, decreasing = TRUE))[1]
    sub <- agg[as.character(agg$.taxon) == top_taxon, ]
    val <- stats::setNames(sub$Abundance, as.character(sub$Sample))
    return(names(sort(val[samples], decreasing = TRUE)))
  }
  if (order.by %in% names(agg)) {
    meta_vals <- tapply(as.character(agg[[order.by]]),
                        as.character(agg$Sample),
                        function(x) x[1])
    return(names(sort(meta_vals[samples])))
  }
  samples
}

#' Stacked bar or area chart
#' @noRd
.mp_taxa_bar <- function(agg, palette, ylab, other.label, area = FALSE,
                         dots = list()) {
  ntax <- nlevels(agg$.taxon)
  cols <- .mp_taxa_colors(palette, agg$.taxon, other.label)

  base_aes <- ggplot2::aes(x = .data$Sample, y = .data$Abundance,
                           fill = .data$.taxon)
  geom <- if (area) {
    .mp_geom(ggplot2::geom_area,
             params = list(mapping = ggplot2::aes(group = .data$.taxon),
                           position = "fill", colour = NA),
             dots = dots)
  } else {
    .mp_geom(ggplot2::geom_col, params = list(width = 0.95), dots = dots)
  }
  ggplot2::ggplot(agg, base_aes) +
    geom +
    ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = ylab)
}

#' Abundance heatmap
#' @noRd
.mp_taxa_heatmap <- function(agg, palette, ylab) {
  cols <- .mp_resolve_palette(palette, 256, type = "sequential")
  ggplot2::ggplot(agg, ggplot2::aes(x = .data$Sample, y = .data$.taxon,
                                    fill = .data$Abundance)) +
    ggplot2::geom_tile(colour = "grey95", linewidth = 0.1) +
    ggplot2::scale_fill_gradientn(colours = cols, name = ylab) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(x = NULL, y = NULL)
}

#' Bubble plot (size = abundance)
#' @noRd
.mp_taxa_bubble <- function(agg, palette, ylab) {
  cols <- .mp_resolve_palette(palette, 256, type = "sequential")
  agg2 <- agg[agg$Abundance > 0, , drop = FALSE]
  ggplot2::ggplot(agg2, ggplot2::aes(x = .data$Sample, y = .data$.taxon,
                                     size = .data$Abundance,
                                     colour = .data$Abundance)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradientn(colours = cols, name = ylab) +
    ggplot2::scale_size_area(max_size = 8, name = ylab) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::guides(colour = ggplot2::guide_legend(),
                    size = ggplot2::guide_legend()) +
    ggplot2::labs(x = NULL, y = NULL)
}

#' Treemap (requires treemapify)
#' @noRd
.mp_taxa_treemap <- function(agg, palette, tax.level) {
  .mp_require("treemapify", "Treemap composition plots")
  tot <- stats::aggregate(Abundance ~ .taxon, data = agg, FUN = sum)
  cols <- .mp_taxa_colors(palette, tot$.taxon, "Other")
  ggplot2::ggplot(tot, ggplot2::aes(area = .data$Abundance,
                                    fill = .data$.taxon,
                                    label = .data$.taxon)) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(colour = "white", place = "centre",
                                  grow = FALSE, reflow = TRUE) +
    ggplot2::scale_fill_manual(values = cols, name = tax.level)
}

#' Assign taxon colours with a muted grey for the "Other" bucket
#' @noRd
.mp_taxa_colors <- function(palette, taxon_factor, other.label) {
  lvls <- levels(taxon_factor)
  non_other <- setdiff(lvls, other.label)
  cols <- .mp_resolve_palette(palette, length(non_other), type = "qualitative")
  out <- stats::setNames(cols, non_other)
  if (other.label %in% lvls) out[other.label] <- "grey70"
  out[lvls]
}
