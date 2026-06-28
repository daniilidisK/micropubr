# Functional pathway profiles ---------------------------------------------

#' Plot functional pathway profiles
#'
#' Visualises functional-prediction tables such as PICRUSt2 pathway abundances
#' or HUMAnN pathway/gene-family tables. Accepts either a
#' [phyloseq::phyloseq] object (whose "taxa" are pathways) or a plain
#' features-by-samples matrix / data frame plus optional sample metadata.
#'
#' @param object A [phyloseq::phyloseq] object, or a matrix / data frame with
#'   pathways in rows and samples in columns.
#' @param metadata Optional sample metadata data frame (rows = samples) used
#'   when `object` is a matrix.
#' @param plot.type One of `"heatmap"`, `"bar"` or `"bubble"`.
#' @param group.by Sample metadata column used to facet / summarise.
#' @param top.n Number of most-abundant pathways to display.
#' @param category Optional column (when `object` carries feature metadata, or
#'   a named vector) grouping pathways into higher-level categories for the
#'   bar plot (e.g. KEGG level 1/2).
#' @param relative Convert to per-sample relative abundance.
#' @param scale Row scaling for the heatmap: `"none"` or `"row"` (z-score).
#' @param palette Palette name or colour vector.
#' @param theme A `ggplot2` theme; defaults to [mp_theme_publication()].
#' @param plot.title,plot.subtitle,plot.caption Optional labels.
#' @param font.size Base font size.
#' @param interactive If `TRUE`, return an interactive `plotly` widget instead
#'   of a static `ggplot` (requires `plotly`). See [mp_alpha()].
#' @param verbose Print the steps used.
#'
#' @return A `ggplot` object, or a `plotly` widget when `interactive = TRUE`.
#' @export
#' @examples
#' set.seed(1)
#' m <- matrix(rpois(20 * 12, 50), nrow = 20,
#'             dimnames = list(paste0("path", 1:20), paste0("S", 1:12)))
#' meta <- data.frame(Group = rep(c("A", "B"), each = 6),
#'                    row.names = colnames(m))
#' mp_function(m, metadata = meta, plot.type = "heatmap", group.by = "Group")
mp_function <- function(object,
                        metadata = NULL,
                        plot.type = c("heatmap", "bar", "bubble"),
                        group.by = NULL,
                        top.n = 25,
                        category = NULL,
                        relative = TRUE,
                        scale = c("none", "row"),
                        palette = NULL,
                        theme = NULL,
                        plot.title = NULL,
                        plot.subtitle = NULL,
                        plot.caption = NULL,
                        font.size = getOption("micropubr.font.size", 11),
                        interactive = FALSE,
                        verbose = getOption("micropubr.verbose", FALSE)) {
  plot.type <- match.arg(plot.type)
  scale <- match.arg(scale)
  theme <- theme %||% mp_theme_publication(font.size = font.size)

  parsed <- .mp_function_input(object, metadata)
  mat <- parsed$mat          # pathways x samples
  meta <- parsed$meta

  if (relative) {
    cs <- colSums(mat)
    cs[cs == 0] <- 1
    mat <- sweep(mat, 2, cs, "/")
  }
  keep <- utils::head(names(sort(rowMeans(mat), decreasing = TRUE)), top.n)
  mat <- mat[keep, , drop = FALSE]

  if (!is.null(group.by) && !group.by %in% names(meta)) {
    rlang::abort(sprintf("`group.by = \"%s\"` not found in metadata.",
                         group.by), class = "micropubr_bad_var")
  }
  if (verbose) {
    rlang::inform(sprintf("mp_function(): %s of top %d pathways.",
                          plot.type, length(keep)))
  }

  long <- .mp_function_long(mat, meta, group.by)

  p <- switch(
    plot.type,
    heatmap = .mp_function_heatmap(mat, meta, group.by, scale, palette,
                                   font.size),
    bar     = .mp_function_bar(long, group.by, category, palette),
    bubble  = .mp_function_bubble(long, group.by, palette)
  )

  p <- p + theme +
    ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                  caption = plot.caption %||%
                    sprintf("Top %d pathways by mean abundance.", length(keep)))
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

#' Normalise mp_function input into a matrix + metadata
#' @noRd
.mp_function_input <- function(object, metadata) {
  if (methods::is(object, "phyloseq") || methods::is(object, "otu_table")) {
    mat <- t(.mp_otu_matrix(object, samples_as_rows = TRUE))  # features x samples
    meta <- .mp_sample_data(object)
    return(list(mat = mat, meta = meta))
  }
  mat <- as.matrix(object)
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("feature", seq_len(nrow(mat)))
  }
  meta <- if (!is.null(metadata)) {
    as.data.frame(metadata, stringsAsFactors = FALSE)
  } else {
    data.frame(row.names = colnames(mat))
  }
  meta$.sample <- rownames(meta)
  list(mat = mat, meta = meta)
}

#' Melt a feature matrix to long format, merging sample metadata
#' @noRd
.mp_function_long <- function(mat, meta, group.by) {
  long <- data.frame(
    pathway = rep(rownames(mat), times = ncol(mat)),
    .sample = rep(colnames(mat), each = nrow(mat)),
    value = as.numeric(mat),
    stringsAsFactors = FALSE
  )
  if (!is.null(group.by) && group.by %in% names(meta)) {
    long[[group.by]] <- meta[long$.sample, group.by]
  }
  long
}

#' @noRd
.mp_function_heatmap <- function(mat, meta, group.by, scale, palette,
                                 font.size) {
  if (scale == "row") {
    mat <- t(scale(t(mat)))
    mat[is.nan(mat)] <- 0
    fill_name <- "z-score"
    cols <- .mp_resolve_palette(palette, 256, type = "diverging")
  } else {
    fill_name <- "Abundance"
    cols <- .mp_resolve_palette(palette, 256, type = "sequential")
  }
  # Order samples by group for readability.
  if (!is.null(group.by)) {
    ord <- order(meta[colnames(mat), group.by])
    mat <- mat[, ord, drop = FALSE]
  }
  long <- data.frame(
    pathway = factor(rep(rownames(mat), times = ncol(mat)),
                     levels = rev(rownames(mat))),
    sample = factor(rep(colnames(mat), each = nrow(mat)),
                    levels = colnames(mat)),
    value = as.numeric(mat), stringsAsFactors = FALSE)
  if (!is.null(group.by)) long[[group.by]] <- meta[as.character(long$sample),
                                                   group.by]

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$sample, y = .data$pathway,
                                          fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = cols, name = fill_name) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90, hjust = 1, vjust = 0.5, size = font.size - 3))
  if (!is.null(group.by)) {
    p <- p + ggplot2::facet_grid(
      stats::as.formula(paste("~", sprintf("`%s`", group.by))),
      scales = "free_x", space = "free_x")
  }
  p
}

#' @noRd
.mp_function_bar <- function(long, group.by, category, palette) {
  x_var <- if (!is.null(group.by)) group.by else ".sample"
  agg <- stats::aggregate(
    stats::as.formula(sprintf("value ~ pathway + `%s`", x_var)),
    data = long, FUN = mean, na.rm = TRUE)
  n <- length(unique(agg$pathway))
  cols <- .mp_resolve_palette(palette, n, type = "qualitative")
  ggplot2::ggplot(agg, ggplot2::aes(x = .data[[x_var]], y = .data$value,
                                    fill = .data$pathway)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = cols, name = "Pathway") +
    ggplot2::labs(x = NULL, y = "Mean relative abundance") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1))
}

#' @noRd
.mp_function_bubble <- function(long, group.by, palette) {
  x_var <- if (!is.null(group.by)) group.by else ".sample"
  agg <- stats::aggregate(
    stats::as.formula(sprintf("value ~ pathway + `%s`", x_var)),
    data = long, FUN = mean, na.rm = TRUE)
  cols <- .mp_resolve_palette(palette, 256, type = "sequential")
  ggplot2::ggplot(agg, ggplot2::aes(x = .data[[x_var]], y = .data$pathway,
                                    size = .data$value,
                                    colour = .data$value)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradientn(colours = cols, name = "Abundance") +
    ggplot2::scale_size_area(max_size = 8, name = "Abundance") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
