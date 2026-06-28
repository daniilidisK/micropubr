# Phylogenetic trees ------------------------------------------------------

#' Plot an annotated phylogenetic tree
#'
#' Renders the phylogenetic tree stored in a [phyloseq::phyloseq] object using
#' `ggtree`, with optional tip colouring by taxonomy and abundance rings drawn
#' with `ggtreeExtra`. To keep large trees legible, agglomerate with
#' `tax.level` and/or restrict to the most abundant tips with `top.n`.
#'
#' @param object A [phyloseq::phyloseq] object containing a `phy_tree`.
#' @param layout Tree layout: `"rectangular"`, `"circular"`, `"fan"` or
#'   `"radial"`.
#' @param annotate Annotation layers to draw: any of `"abundance"` and
#'   `"taxonomy"`.
#' @param color.by Taxonomic rank used to colour tip points (e.g. `"Phylum"`).
#' @param group.by Metadata column; when set, the abundance ring becomes a
#'   per-group heatmap instead of a single bar.
#' @param tax.level Optional rank to agglomerate to before plotting.
#' @param top.n Optional number of most-abundant tips to retain.
#' @param tip.label Show tip labels.
#' @param tip.label.size Tip-label font size.
#' @param highlight.clades Named list mapping a taxonomic rank to one or more
#'   values whose clades should be highlighted, e.g.
#'   `list(Phylum = "Firmicutes")`.
#' @param open.angle Opening angle (degrees) for circular/fan layouts.
#' @param palette Palette name or colour vector.
#' @param theme A `ggtree` theme; defaults to [ggtree::theme_tree()].
#' @param plot.title,plot.subtitle,plot.caption Optional labels.
#' @param font.size Base font size.
#' @param verbose Print the steps used.
#'
#' @return A `ggtree`/`ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' gp <- phyloseq::tax_glom(GlobalPatterns, "Phylum")
#' mp_tree(gp, layout = "circular", color.by = "Phylum",
#'         group.by = "SampleType")
#' }
mp_tree <- function(object,
                    layout = c("rectangular", "circular", "fan", "radial"),
                    annotate = c("abundance", "taxonomy"),
                    color.by = NULL,
                    group.by = NULL,
                    tax.level = NULL,
                    top.n = NULL,
                    tip.label = FALSE,
                    tip.label.size = 2,
                    highlight.clades = NULL,
                    open.angle = 0,
                    palette = NULL,
                    theme = NULL,
                    plot.title = NULL,
                    plot.subtitle = NULL,
                    plot.caption = NULL,
                    font.size = getOption("micropubr.font.size", 11),
                    verbose = getOption("micropubr.verbose", FALSE)) {
  .mp_check_phyloseq(object, require_tree = TRUE)
  .mp_require("ggtree", "Phylogenetic tree plots", bioc = TRUE)
  layout <- match.arg(layout)
  annotate <- match.arg(annotate, several.ok = TRUE)
  group.by <- .mp_check_var(object, group.by, "group.by")

  if (!is.null(tax.level)) object <- .mp_agglomerate(object, tax.level)
  if (!is.null(top.n)) {
    object <- phyloseq::prune_taxa(.mp_top_taxa(object, top.n), object)
  }
  if (verbose) {
    rlang::inform(sprintf("mp_tree(): %s layout, %d tips, annotate = %s.",
                          layout, phyloseq::ntaxa(object),
                          paste(annotate, collapse = "+")))
  }

  tree <- phyloseq::phy_tree(object)
  tip_df <- .mp_tree_tip_data(object)

  gg_layout <- if (layout == "radial") "circular" else layout
  p <- ggtree::ggtree(tree, layout = gg_layout, open.angle = open.angle)
  if (layout == "radial") p <- p + ggtree::layout_dendrogram()
  # Attach tip metadata for aesthetic mappings.
  p <- p %<+% tip_df

  # Highlight clades by taxonomic rank value.
  if (!is.null(highlight.clades)) {
    p <- .mp_tree_highlight(p, object, tip_df, highlight.clades, palette)
  }

  # Tip colouring by taxonomy.
  if ("taxonomy" %in% annotate && !is.null(color.by)) {
    n <- length(unique(stats::na.omit(tip_df[[color.by]])))
    cols <- .mp_resolve_palette(palette, n, type = "qualitative")
    p <- p + ggtree::geom_tippoint(
      ggplot2::aes(colour = .data[[color.by]]), size = 1.8, na.rm = TRUE) +
      ggplot2::scale_colour_manual(values = cols, name = color.by,
                                   na.value = "grey80")
  }

  if (tip.label) {
    p <- p + ggtree::geom_tiplab(size = tip.label.size, offset = 0.002)
  }

  # Abundance ring (bar, or per-group heatmap).
  if ("abundance" %in% annotate) {
    p <- .mp_tree_abundance_ring(p, object, group.by, palette)
  }

  p <- p + ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                         caption = plot.caption)
  if (!is.null(theme)) p <- p + theme
  p
}

#' Build a tip-level metadata data.frame (label + taxonomy + total abundance)
#' @noRd
.mp_tree_tip_data <- function(object) {
  tax <- phyloseq::tax_table(object, errorIfNULL = FALSE)
  rel <- .mp_to_relative(object)
  ab <- phyloseq::taxa_sums(rel)
  df <- data.frame(label = phyloseq::taxa_names(object),
                   abundance = as.numeric(ab[phyloseq::taxa_names(object)]),
                   stringsAsFactors = FALSE)
  if (!is.null(tax)) {
    taxdf <- as.data.frame(as(tax, "matrix"), stringsAsFactors = FALSE)
    taxdf$label <- rownames(taxdf)
    df <- merge(df, taxdf, by = "label", all.x = TRUE, sort = FALSE)
  }
  df
}

#' Add an abundance ring via ggtreeExtra::geom_fruit
#' @noRd
.mp_tree_abundance_ring <- function(p, object, group.by, palette) {
  .mp_require("ggtreeExtra", "Abundance rings on trees", bioc = TRUE)
  # NOTE: ggtreeExtra calls rlang::as_name() on the mapping, so the aesthetics
  # must be *bare symbols* (not `.data[[...]]`) and `geom` a bare name.
  if (is.null(group.by)) {
    rel <- .mp_to_relative(object)
    ab <- data.frame(label = phyloseq::taxa_names(object),
                     value = as.numeric(phyloseq::taxa_sums(rel)),
                     stringsAsFactors = FALSE)
    return(p + ggtreeExtra::geom_fruit(
      data = ab, geom = geom_col,
      mapping = ggplot2::aes(y = label, x = value),
      orientation = "y", fill = "#3283FE", offset = 0.05, pwidth = 0.4))
  }
  # Per-group mean relative abundance heatmap ring; the group column is given
  # the fixed name `group_value` so it can be referenced as a bare symbol.
  ab_long <- .mp_tree_group_abundance(object, group.by)
  cols <- .mp_resolve_palette(palette, 256, type = "sequential")
  p + ggtreeExtra::geom_fruit(
    data = ab_long, geom = geom_tile,
    mapping = ggplot2::aes(y = label, x = group_value, fill = value),
    offset = 0.05, pwidth = 0.4) +
    ggplot2::scale_fill_gradientn(colours = cols, name = "Rel. abund.")
}

#' Long-format mean relative abundance per taxon per group
#' @noRd
.mp_tree_group_abundance <- function(object, group.by) {
  rel <- .mp_to_relative(object)
  mat <- .mp_otu_matrix(rel, samples_as_rows = TRUE)  # samples x taxa
  meta <- .mp_sample_data(object)
  grp <- meta[rownames(mat), group.by]
  means <- t(apply(mat, 2, function(col) tapply(col, grp, mean, na.rm = TRUE)))
  data.frame(
    label = rep(rownames(means), times = ncol(means)),
    group_value = rep(colnames(means), each = nrow(means)),
    value = as.numeric(means),
    stringsAsFactors = FALSE
  )
}

#' Highlight clades sharing a taxonomic rank value
#' @noRd
.mp_tree_highlight <- function(p, object, tip_df, highlight.clades, palette) {
  .mp_require("tidytree", "Clade highlighting", bioc = TRUE)
  rank <- names(highlight.clades)[1]
  values <- highlight.clades[[1]]
  cols <- .mp_resolve_palette(palette, length(values), type = "qualitative")
  for (i in seq_along(values)) {
    tips <- tip_df$label[tip_df[[rank]] %in% values[i]]
    if (length(tips) < 2) next
    node <- tryCatch(tidytree::MRCA(p, tips), error = function(e) NULL)
    if (!is.null(node)) {
      p <- p + ggtree::geom_hilight(node = node, fill = cols[i], alpha = 0.25)
    }
  }
  p
}

# `%<+%` is ggtree's data-attachment operator; re-export for internal use.
#' @noRd
`%<+%` <- function(pg, data) {
  op <- getExportedValue("ggtree", "%<+%")
  op(pg, data)
}
