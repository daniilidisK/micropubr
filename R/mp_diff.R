# Differential abundance --------------------------------------------------

#' Plot differential-abundance results
#'
#' Visualises the output of a differential-abundance analysis as a volcano,
#' MA, forest or ranked-dot plot, or compares two analyses with a concordance
#' plot. Results from `DESeq2`, `ALDEx2` and `ANCOM-BC`/`ANCOM-BC2` are parsed
#' automatically; any other tool can be supplied as a data frame with explicit
#' column mappings.
#'
#' @param result A differential-abundance result: a `DESeqResults` object, an
#'   `ALDEx2` data frame, an `ANCOM-BC(2)` result list, or a generic data
#'   frame (see `feature`/`log2fc`/`pvalue`/`padj`).
#' @param plot.type One of `"volcano"`, `"ma"`, `"forest"`, `"dot"` or
#'   `"concordance"`.
#' @param fc.threshold Absolute log2 fold-change cut-off for significance.
#' @param p.threshold Adjusted p-value cut-off.
#' @param y.value For volcano plots, plot `-log10()` of `"padj"` (default) or
#'   `"pvalue"` on the y-axis.
#' @param label.top Number of most-significant features to label.
#' @param color.up,color.down,color.ns Colours for enriched, depleted and
#'   non-significant features.
#' @param result2 A second result (same parsing rules) for
#'   `plot.type = "concordance"`.
#' @param method.names Length-2 labels for the two methods in a concordance
#'   plot.
#' @param feature,log2fc,pvalue,padj,mean,lfcSE Optional column-name overrides
#'   when `result` is a generic data frame.
#' @param theme A `ggplot2` theme; defaults to the active theme (see
#'   [mp_set_defaults()], otherwise [mp_theme_publication()]).
#' @param ... Additional arguments passed to the point geom.
#' @param plot.title,plot.subtitle,plot.caption Optional labels.
#' @param font.size Base font size.
#' @param return.data If `TRUE`, return the parsed differential-abundance
#'   data frame (with a `direction` column) instead of a plot.
#' @param interactive If `TRUE`, return an interactive `plotly` widget instead
#'   of a static `ggplot` (requires `plotly`). See [mp_alpha()].
#' @param verbose Print the detected result type and counts.
#'
#' @return A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
#'   data frame when `return.data = TRUE`.
#' @export
#' @examples
#' df <- data.frame(
#'   feature = paste0("ASV", 1:50),
#'   log2FC = rnorm(50, sd = 2),
#'   pvalue = runif(50)^2,
#'   padj = p.adjust(runif(50)^2, "BH")
#' )
#' mp_diff(df, plot.type = "volcano")
mp_diff <- function(result,
                    plot.type = c("volcano", "ma", "forest", "dot",
                                  "concordance"),
                    fc.threshold = 1,
                    p.threshold = 0.05,
                    y.value = c("padj", "pvalue"),
                    label.top = 10,
                    color.up = "#B2182B",
                    color.down = "#2166AC",
                    color.ns = "grey75",
                    result2 = NULL,
                    method.names = c("Method 1", "Method 2"),
                    feature = NULL, log2fc = NULL, pvalue = NULL,
                    padj = NULL, mean = NULL, lfcSE = NULL,
                    theme = NULL,
                    ...,
                    plot.title = NULL,
                    plot.subtitle = NULL,
                    plot.caption = NULL,
                    font.size = getOption("micropubr.font.size", 11),
                    interactive = FALSE,
                    return.data = FALSE,
                    verbose = getOption("micropubr.verbose", FALSE)) {
  plot.type <- match.arg(plot.type)
  y.value <- match.arg(y.value)
  theme <- .mp_theme(theme, font.size)

  da <- .mp_parse_da(result, feature, log2fc, pvalue, padj, mean, lfcSE,
                     verbose = verbose)
  da$direction <- .mp_da_direction(da, fc.threshold, p.threshold)

  if (isTRUE(return.data)) return(da)

  n_up <- sum(da$direction == "Enriched", na.rm = TRUE)
  n_dn <- sum(da$direction == "Depleted", na.rm = TRUE)
  if (verbose) {
    rlang::inform(sprintf("mp_diff(): %d enriched, %d depleted (padj < %g, |log2FC| > %g).",
                          n_up, n_dn, p.threshold, fc.threshold))
  }
  subtitle <- plot.subtitle %||%
    sprintf("%d enriched \u00b7 %d depleted (padj < %g, |log2FC| > %g)",
            n_up, n_dn, p.threshold, fc.threshold)

  cols <- c(Enriched = color.up, Depleted = color.down, NS = color.ns)

  dots <- list(...)
  p <- switch(
    plot.type,
    volcano = .mp_da_volcano(da, cols, fc.threshold, p.threshold,
                             y.value, label.top, font.size, dots = dots),
    ma      = .mp_da_ma(da, cols, label.top, font.size, dots = dots),
    forest  = .mp_da_forest(da, cols, font.size),
    dot     = .mp_da_dot(da, cols, font.size),
    concordance = .mp_da_concordance(
      da, result2, method.names, fc.threshold, p.threshold,
      feature, log2fc, pvalue, padj, mean, lfcSE, font.size, verbose)
  )

  p <- p + theme +
    ggplot2::labs(title = plot.title, subtitle = subtitle,
                  caption = plot.caption)
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

#' Parse a differential-abundance result into a standard data frame
#'
#' @return A data.frame with columns `feature`, `log2FC`, `pvalue`, `padj`
#'   and, where available, `mean` and `lfcSE`.
#' @noRd
.mp_parse_da <- function(result, feature = NULL, log2fc = NULL,
                         pvalue = NULL, padj = NULL, mean = NULL,
                         lfcSE = NULL, verbose = FALSE) {
  cls <- class(result)[1]

  # DESeq2::results() output ------------------------------------------------
  if (inherits(result, "DESeqResults") || cls == "DESeqResults") {
    df <- as.data.frame(result)
    out <- data.frame(
      feature = rownames(df),
      log2FC  = df$log2FoldChange,
      pvalue  = df$pvalue,
      padj    = df$padj,
      mean    = df$baseMean,
      lfcSE   = df$lfcSE,
      stringsAsFactors = FALSE
    )
    if (verbose) rlang::inform("mp_diff(): detected DESeq2 result.")
    return(out[stats::complete.cases(out$log2FC), ])
  }

  # ANCOM-BC / ANCOM-BC2 (list with $res) ----------------------------------
  if (is.list(result) && !is.data.frame(result) && "res" %in% names(result)) {
    res <- result$res
    feat <- res[["taxon"]] %||% rownames(res)
    out <- data.frame(
      feature = feat,
      log2FC  = .mp_first_col(res, "^lfc"),
      pvalue  = .mp_first_col(res, "^p_?val|^p_"),
      padj    = .mp_first_col(res, "^q_?val|^q_"),
      lfcSE   = .mp_first_col(res, "^se"),
      stringsAsFactors = FALSE
    )
    if (verbose) rlang::inform("mp_diff(): detected ANCOM-BC result.")
    return(out[stats::complete.cases(out$log2FC), ])
  }

  # Generic data frame (incl. ALDEx2) --------------------------------------
  if (is.data.frame(result)) {
    df <- as.data.frame(result)
    feat <- if (!is.null(feature)) df[[feature]] else
      .mp_pick_col(df, c("feature", "taxon", "OTU"), rownames(df))
    l2fc <- if (!is.null(log2fc)) df[[log2fc]] else
      .mp_pick_col(df, c("log2FC", "log2FoldChange", "effect", "diff.btw",
                         "logFC", "coef"))
    pv <- if (!is.null(pvalue)) df[[pvalue]] else
      .mp_pick_col(df, c("pvalue", "p", "pval", "p.value", "we.ep",
                         "wi.ep", "P.Value"))
    pa <- if (!is.null(padj)) df[[padj]] else
      .mp_pick_col(df, c("padj", "p.adj", "FDR", "qvalue", "q", "we.eBH",
                         "wi.eBH", "adj.P.Val"))
    mn <- if (!is.null(mean)) df[[mean]] else
      .mp_pick_col(df, c("mean", "baseMean", "rab.all", "AveExpr"),
                   default = NA_real_)
    se <- if (!is.null(lfcSE)) df[[lfcSE]] else
      .mp_pick_col(df, c("lfcSE", "se", "SE", "stderr"), default = NA_real_)
    if (is.null(l2fc)) {
      rlang::abort(
        c("Could not find a log fold-change column in the data frame.",
          i = "Pass `log2fc = \"<colname>\"` to specify it explicitly."),
        class = "micropubr_da_parse")
    }
    out <- data.frame(feature = as.character(feat), log2FC = l2fc,
                      pvalue = pv, padj = pa, mean = mn, lfcSE = se,
                      stringsAsFactors = FALSE)
    if (verbose) rlang::inform("mp_diff(): parsed generic data frame.")
    return(out[stats::complete.cases(out$log2FC), ])
  }

  rlang::abort(sprintf("Don't know how to parse a <%s> result.", cls),
               class = "micropubr_da_parse")
}

#' First column matching any regex pattern; NA vector otherwise
#' @noRd
.mp_first_col <- function(df, pattern) {
  hit <- grep(pattern, names(df), ignore.case = TRUE, value = TRUE)
  if (length(hit)) df[[hit[1]]] else rep(NA_real_, nrow(df))
}

#' Pick the first existing column from candidates
#' @noRd
.mp_pick_col <- function(df, candidates, default = NULL) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) df[[hit[1]]] else default
}

#' Classify features as Enriched / Depleted / NS
#' @noRd
.mp_da_direction <- function(da, fc.threshold, p.threshold) {
  padj <- da$padj
  padj[is.na(padj)] <- 1
  sig <- padj < p.threshold & abs(da$log2FC) > fc.threshold
  out <- rep("NS", nrow(da))
  out[sig & da$log2FC > 0] <- "Enriched"
  out[sig & da$log2FC < 0] <- "Depleted"
  factor(out, levels = c("Enriched", "Depleted", "NS"))
}

#' Pick the `label.top` most significant features for labelling
#' @noRd
.mp_da_top <- function(da, label.top) {
  sig <- da[da$direction != "NS", , drop = FALSE]
  if (nrow(sig) == 0 || label.top <= 0) return(sig[0, ])
  key <- sig$padj
  key[is.na(key)] <- 1
  sig[order(key)[seq_len(min(label.top, nrow(sig)))], , drop = FALSE]
}

#' @noRd
.mp_da_volcano <- function(da, cols, fc.threshold, p.threshold, y.value,
                           label.top, font.size, dots = list()) {
  da$y <- -log10(pmax(da[[y.value]], .Machine$double.xmin))
  ylab <- if (y.value == "padj") expression(-log[10]~adjusted~italic(p)) else
    expression(-log[10]~italic(p))
  p <- ggplot2::ggplot(da, ggplot2::aes(x = .data$log2FC, y = .data$y,
                                        colour = .data$direction)) +
    ggplot2::geom_vline(xintercept = c(-fc.threshold, fc.threshold),
                        linetype = "dashed", colour = "grey70",
                        linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(p.threshold),
                        linetype = "dashed", colour = "grey70",
                        linewidth = 0.3) +
    .mp_geom(ggplot2::geom_point,
             params = list(alpha = 0.8, size = 1.8), dots = dots) +
    ggplot2::scale_colour_manual(values = cols, name = NULL,
                                 drop = FALSE) +
    ggplot2::labs(x = expression(log[2]~fold-change), y = ylab)
  .mp_add_da_labels(p, da, label.top, font.size,
                    ggplot2::aes(label = .data$feature))
}

#' @noRd
.mp_da_ma <- function(da, cols, label.top, font.size, dots = list()) {
  if (all(is.na(da$mean))) {
    rlang::abort("MA plot requires a mean-abundance column.",
                 class = "micropubr_da_no_mean")
  }
  da$x <- log10(da$mean + 1)
  p <- ggplot2::ggplot(da, ggplot2::aes(x = .data$x, y = .data$log2FC,
                                        colour = .data$direction)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    .mp_geom(ggplot2::geom_point,
             params = list(alpha = 0.8, size = 1.8), dots = dots) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
    ggplot2::labs(x = expression(log[10]~mean~abundance),
                  y = expression(log[2]~fold-change))
  .mp_add_da_labels(p, da, label.top, font.size,
                    ggplot2::aes(label = .data$feature))
}

#' @noRd
.mp_da_forest <- function(da, cols, font.size) {
  sig <- da[da$direction != "NS", , drop = FALSE]
  if (nrow(sig) == 0) {
    rlang::abort("No significant features to display in the forest plot.",
                 class = "micropubr_da_empty")
  }
  if (all(is.na(sig$lfcSE))) sig$lfcSE <- 0
  sig$feature <- stats::reorder(sig$feature, sig$log2FC)
  ggplot2::ggplot(sig, ggplot2::aes(x = .data$log2FC, y = .data$feature,
                                    colour = .data$direction)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$log2FC - 1.96 * .data$lfcSE,
                   xmax = .data$log2FC + 1.96 * .data$lfcSE),
      orientation = "y", width = 0.25, linewidth = 0.4) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
    ggplot2::labs(x = expression(log[2]~fold-change~(95*"%"~CI)), y = NULL)
}

#' @noRd
.mp_da_dot <- function(da, cols, font.size) {
  sig <- da[da$direction != "NS", , drop = FALSE]
  if (nrow(sig) == 0) {
    rlang::abort("No significant features to display in the dot plot.",
                 class = "micropubr_da_empty")
  }
  sig$feature <- stats::reorder(sig$feature, sig$log2FC)
  sig$absp <- -log10(pmax(sig$padj, .Machine$double.xmin))
  ggplot2::ggplot(sig, ggplot2::aes(x = .data$log2FC, y = .data$feature,
                                    colour = .data$direction,
                                    size = .data$absp)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
    ggplot2::scale_size_area(max_size = 6,
                             name = expression(-log[10]~padj)) +
    ggplot2::labs(x = expression(log[2]~fold-change), y = NULL)
}

#' @noRd
.mp_da_concordance <- function(da1, result2, method.names, fc.threshold,
                               p.threshold, feature, log2fc, pvalue, padj,
                               mean, lfcSE, font.size, verbose) {
  if (is.null(result2)) {
    rlang::abort("`plot.type = \"concordance\"` requires `result2`.",
                 class = "micropubr_da_no_result2")
  }
  da2 <- .mp_parse_da(result2, feature, log2fc, pvalue, padj, mean, lfcSE,
                      verbose = verbose)
  da2$direction <- .mp_da_direction(da2, fc.threshold, p.threshold)

  m <- merge(
    data.frame(feature = da1$feature, fc1 = da1$log2FC, d1 = da1$direction),
    data.frame(feature = da2$feature, fc2 = da2$log2FC, d2 = da2$direction),
    by = "feature"
  )
  m$agreement <- ifelse(
    m$d1 != "NS" & m$d2 != "NS", "Both significant",
    ifelse(m$d1 != "NS" | m$d2 != "NS", "One significant", "Neither")
  )
  m$agreement <- factor(m$agreement,
                        levels = c("Both significant", "One significant",
                                   "Neither"))
  n_both <- sum(m$agreement == "Both significant")
  rho <- suppressWarnings(stats::cor(m$fc1, m$fc2, method = "spearman",
                                     use = "complete.obs"))
  ggplot2::ggplot(m, ggplot2::aes(x = .data$fc1, y = .data$fc2,
                                  colour = .data$agreement)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey70") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_point(alpha = 0.8, size = 1.8) +
    ggplot2::scale_colour_manual(
      values = c("Both significant" = "#B2182B",
                 "One significant" = "#E69F00", "Neither" = "grey75"),
      name = NULL, drop = FALSE) +
    ggplot2::labs(
      x = sprintf("log2FC \u2014 %s", method.names[1]),
      y = sprintf("log2FC \u2014 %s", method.names[2]),
      caption = sprintf("Spearman r = %.2f \u00b7 %d features significant in both",
                        rho, n_both))
}

#' Add repelled labels for the top features
#' @noRd
.mp_add_da_labels <- function(p, da, label.top, font.size, mapping) {
  top <- .mp_da_top(da, label.top)
  if (nrow(top) == 0) return(p)
  .mp_require("ggrepel", "Feature labelling")
  p + ggrepel::geom_text_repel(
    data = top, mapping = mapping, size = font.size / 3.5,
    max.overlaps = 30, show.legend = FALSE, colour = "black",
    min.segment.length = 0
  )
}
