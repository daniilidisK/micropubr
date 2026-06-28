# Co-occurrence / correlation networks -----------------------------------

#' Plot a taxon co-occurrence network
#'
#' Builds a correlation network between taxa and renders it with `ggraph`.
#' Correlations are computed on centred log-ratio (CLR) coordinates by
#' default; edges are filtered by correlation strength and adjusted p-value.
#' SparCC / SPIEC-EASI inference is available when the `SpiecEasi` package is
#' installed.
#'
#' @param object A [phyloseq::phyloseq] object.
#' @param method Association method: `"spearman"`, `"pearson"`, `"sparcc"` or
#'   `"spiec-easi"`. The latter two require the `SpiecEasi` package.
#' @param tax.level Optional rank to agglomerate to before building the
#'   network.
#' @param top.n Keep only the `top.n` most-abundant taxa (for tractability).
#' @param cor.threshold Minimum absolute correlation to draw an edge.
#' @param p.threshold Maximum adjusted p-value for an edge (correlation
#'   methods only).
#' @param p.adjust.method Multiple-testing correction (default `"BH"`).
#' @param layout `ggraph` layout: `"fr"`, `"stress"`, `"kk"` or `"circle"`.
#' @param node.color.by Taxonomic rank used to colour nodes.
#' @param node.size.by Size nodes by `"abundance"` or `"degree"`.
#' @param edge.color.by Colour edges by `"sign"` (default) or `"none"`.
#' @param label.nodes Label nodes with their taxon name.
#' @param label.size Node-label font size.
#' @param palette Palette name or colour vector.
#' @param theme A `ggplot2` theme; defaults to a clean void theme.
#' @param plot.title,plot.subtitle,plot.caption Optional labels.
#' @param font.size Base font size.
#' @param verbose Print the steps used.
#'
#' @return A `ggplot`/`ggraph` object.
#' @export
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' gp <- phyloseq::tax_glom(GlobalPatterns, "Genus")
#' mp_network(gp, top.n = 60, node.color.by = "Phylum")
#' }
mp_network <- function(object,
                       method = c("spearman", "pearson", "sparcc",
                                  "spiec-easi"),
                       tax.level = NULL,
                       top.n = 50,
                       cor.threshold = 0.6,
                       p.threshold = 0.05,
                       p.adjust.method = "BH",
                       layout = c("fr", "stress", "kk", "circle"),
                       node.color.by = NULL,
                       node.size.by = c("abundance", "degree"),
                       edge.color.by = c("sign", "none"),
                       label.nodes = FALSE,
                       label.size = 2.5,
                       palette = NULL,
                       theme = NULL,
                       plot.title = NULL,
                       plot.subtitle = NULL,
                       plot.caption = NULL,
                       font.size = getOption("micropubr.font.size", 11),
                       verbose = getOption("micropubr.verbose", FALSE)) {
  .mp_check_phyloseq(object)
  .mp_require(c("ggraph", "igraph"), "Network plots")
  method <- match.arg(method)
  layout <- match.arg(layout)
  node.size.by <- match.arg(node.size.by)
  edge.color.by <- match.arg(edge.color.by)

  if (!is.null(tax.level)) object <- .mp_agglomerate(object, tax.level)
  if (!is.null(top.n) && top.n < phyloseq::ntaxa(object)) {
    object <- phyloseq::prune_taxa(.mp_top_taxa(object, top.n), object)
  }

  edges <- .mp_network_edges(object, method, cor.threshold, p.threshold,
                             p.adjust.method)
  if (nrow(edges) == 0) {
    rlang::abort(
      c("No edges passed the thresholds.",
        i = sprintf("Try lowering `cor.threshold` (now %.2f) or raising `p.threshold` (now %g).",
                    cor.threshold, p.threshold)),
      class = "micropubr_empty_network")
  }

  nodes <- .mp_network_nodes(object, edges)
  graph <- igraph::graph_from_data_frame(edges, vertices = nodes,
                                         directed = FALSE)
  igraph::V(graph)$degree <- igraph::degree(graph)

  if (verbose) {
    rlang::inform(sprintf(
      "mp_network(): %s, %d nodes, %d edges (|r| >= %.2f, padj < %g).",
      method, igraph::vcount(graph), igraph::ecount(graph),
      cor.threshold, p.threshold))
  }

  theme <- theme %||% (ggplot2::theme_void(base_size = font.size) +
                         ggplot2::theme(
                           plot.caption = ggplot2::element_text(
                             colour = "grey40", hjust = 0)))
  .mp_network_plot(graph, layout, node.color.by, node.size.by,
                   edge.color.by, label.nodes, label.size, palette, theme,
                   plot.title, plot.subtitle, plot.caption)
}

#' Compute an edge list from a correlation or SparCC network
#' @noRd
.mp_network_edges <- function(object, method, cor.threshold, p.threshold,
                              p.adjust.method) {
  if (method %in% c("sparcc", "spiec-easi")) {
    return(.mp_network_spiec(object, method, cor.threshold))
  }
  mat <- .mp_clr(object)                     # samples x taxa
  C <- stats::cor(mat, method = method)
  n <- nrow(mat)
  ut <- upper.tri(C)
  r <- C[ut]
  # t-approximation for correlation p-values, then multiple-testing correction.
  tval <- r * sqrt((n - 2) / (1 - r^2))
  pval <- 2 * stats::pt(-abs(tval), df = n - 2)
  padj <- stats::p.adjust(pval, method = p.adjust.method)

  idx <- which(ut, arr.ind = TRUE)
  taxa <- colnames(C)
  keep <- abs(r) >= cor.threshold & padj < p.threshold & is.finite(r)
  data.frame(
    from = taxa[idx[keep, 1]],
    to   = taxa[idx[keep, 2]],
    r    = r[keep],
    abs_r = abs(r[keep]),
    sign = factor(ifelse(r[keep] > 0, "Positive", "Negative"),
                  levels = c("Positive", "Negative")),
    stringsAsFactors = FALSE
  )
}

#' SparCC / SPIEC-EASI edges via the SpiecEasi package
#' @noRd
.mp_network_spiec <- function(object, method, cor.threshold) {
  .mp_require("SpiecEasi", "SparCC / SPIEC-EASI network inference")
  mat <- round(.mp_otu_matrix(object, samples_as_rows = TRUE))
  taxa <- colnames(mat)
  if (method == "sparcc") {
    sp <- SpiecEasi::sparcc(mat)
    C <- sp$Cor
  } else {
    se <- SpiecEasi::spiec.easi(mat, method = "mb")
    C <- as.matrix(SpiecEasi::symBeta(SpiecEasi::getOptBeta(se)))
  }
  dimnames(C) <- list(taxa, taxa)
  ut <- upper.tri(C)
  r <- C[ut]
  idx <- which(ut, arr.ind = TRUE)
  keep <- abs(r) >= cor.threshold & r != 0
  data.frame(
    from = taxa[idx[keep, 1]], to = taxa[idx[keep, 2]],
    r = r[keep], abs_r = abs(r[keep]),
    sign = factor(ifelse(r[keep] > 0, "Positive", "Negative"),
                  levels = c("Positive", "Negative")),
    stringsAsFactors = FALSE
  )
}

#' Build a node table (connected taxa) with taxonomy and abundance attributes
#' @noRd
.mp_network_nodes <- function(object, edges) {
  ids <- union(edges$from, edges$to)
  rel <- .mp_to_relative(object)
  ab <- phyloseq::taxa_sums(rel) / phyloseq::nsamples(object)
  nodes <- data.frame(name = ids, abundance = as.numeric(ab[ids]),
                      stringsAsFactors = FALSE)
  tax <- phyloseq::tax_table(object, errorIfNULL = FALSE)
  if (!is.null(tax)) {
    taxdf <- as.data.frame(as(tax, "matrix"), stringsAsFactors = FALSE)
    taxdf$name <- rownames(taxdf)
    nodes <- merge(nodes, taxdf, by = "name", all.x = TRUE, sort = FALSE)
  }
  nodes[match(ids, nodes$name), , drop = FALSE]
}

#' Render the network with ggraph
#' @noRd
.mp_network_plot <- function(graph, layout, node.color.by, node.size.by,
                             edge.color.by, label.nodes, label.size, palette,
                             theme, plot.title, plot.subtitle, plot.caption) {
  p <- ggraph::ggraph(graph, layout = layout)

  # Edges.
  if (edge.color.by == "sign") {
    p <- p + ggraph::geom_edge_link(
      ggplot2::aes(edge_colour = .data$sign, edge_width = .data$abs_r),
      alpha = 0.5) +
      ggraph::scale_edge_colour_manual(
        values = c(Positive = "#2166AC", Negative = "#B2182B"),
        name = NULL) +
      ggraph::scale_edge_width(range = c(0.2, 1.5), guide = "none")
  } else {
    p <- p + ggraph::geom_edge_link(
      ggplot2::aes(edge_width = .data$abs_r), colour = "grey70",
      alpha = 0.5) +
      ggraph::scale_edge_width(range = c(0.2, 1.5), guide = "none")
  }

  # Nodes.
  size_var <- if (node.size.by == "degree") "degree" else "abundance"
  node_aes <- ggplot2::aes(size = .data[[size_var]])
  if (!is.null(node.color.by)) {
    node_aes$colour <- ggplot2::aes(colour = .data[[node.color.by]])$colour
  }
  p <- p + ggraph::geom_node_point(node_aes, alpha = 0.9)

  if (!is.null(node.color.by)) {
    n <- length(unique(stats::na.omit(
      igraph::vertex_attr(graph, node.color.by))))
    n <- max(n, 3)
    cols <- .mp_resolve_palette(palette, n, type = "qualitative")
    p <- p + ggplot2::scale_colour_manual(values = cols, name = node.color.by,
                                          na.value = "grey80")
  }
  p <- p + ggplot2::scale_size_continuous(
    range = c(1.5, 7),
    name = if (size_var == "degree") "Degree" else "Mean abund.")

  if (label.nodes) {
    p <- p + ggraph::geom_node_text(
      ggplot2::aes(label = .data$name), size = label.size, repel = TRUE,
      max.overlaps = 30)
  }

  p + theme +
    ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                  caption = plot.caption)
}
