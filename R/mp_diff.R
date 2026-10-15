# Differential abundance --------------------------------------------------

#' Plot differential-abundance results
#'
#' Visualises differential-abundance results as a volcano, MA, forest or
#' ranked-dot plot, as a heatmap of fold changes across contrasts, or -- when
#' several methods were run -- as a consensus matrix showing where the methods
#' agree. Results from [mp_run_da()], `DESeq2`, `ALDEx2` and
#' `ANCOM-BC`/`ANCOM-BC2` are parsed automatically; any other tool can be
#' supplied as a data frame with explicit column mappings, and a named list of
#' results from different tools is stacked for comparison.
#'
#' Results with several contrasts (multi-group designs) or several methods are
#' drawn as facets. ANCOM-BC(2) fold changes, which are natural logarithms,
#' are converted to log2. Taxa that ANCOM-BC2 classifies as structural zeros
#' (absent from a whole group) are not tested; they are counted in the
#' subtitle and shown by `plot.type = "zeros"`. When results carry ANCOM-BC2's
#' pseudo-count sensitivity analysis, taxa that fail it are not counted as
#' significant (as in ANCOM-BC2's own `diff_robust` output); they are counted in
#' the subtitle and marked in heatmaps. With `robust.only = FALSE` they are
#' kept and drawn as open symbols.
#'
#' @param result A differential-abundance result: the output of
#'   [mp_run_da()], a `DESeqResults` object, an `ALDEx2` data frame, an
#'   `ANCOM-BC(2)` result list, a generic data frame (see
#'   `feature`/`log2fc`/`pvalue`/`padj`), or a named list of any of these
#'   (one per method, e.g. `list(DESeq2 = res1, ALDEx2 = res2)`).
#' @param plot.type One of `"volcano"`, `"ma"`, `"forest"`, `"dot"`,
#'   `"heatmap"` (fold change of significant taxa across contrasts),
#'   `"consensus"` (agreement across methods), `"zeros"` (ANCOM-BC2
#'   structural zeros) or `"concordance"` (fold changes of two methods
#'   against each other).
#' @param fc.threshold Absolute log2 fold-change cut-off for significance.
#' @param p.threshold Adjusted p-value cut-off.
#' @param y.value For volcano plots, plot `-log10()` of `"padj"` (default) or
#'   `"pvalue"` on the y-axis.
#' @param label.top Number of most-significant features to label per panel.
#' @param color.up,color.down,color.ns Colours for enriched, depleted and
#'   non-significant features.
#' @param result2 A second result (same parsing rules) for
#'   `plot.type = "concordance"`. Not needed when `result` already holds
#'   exactly two methods.
#' @param method.names Length-2 labels for the two methods in a concordance
#'   plot.
#' @param contrast Optional character vector of contrasts (as in the
#'   `contrast` column, e.g. `"Treatment vs Control"`, or an ANCOM-BC2 term
#'   such as `"groupTreatment"`) to keep.
#' @param robust.only Treat significant features that fail ANCOM-BC2's
#'   pseudo-count sensitivity analysis as not significant (the default, which
#'   matches ANCOM-BC2's `diff_robust`). Results without sensitivity
#'   information are unaffected.
#' @param top.n Maximum number of features in forest, dot, heatmap and
#'   consensus plots (the most significant are kept; per panel for forest and
#'   dot plots).
#' @param feature,log2fc,pvalue,padj,mean,lfcSE Optional column-name overrides
#'   when `result` is a generic data frame.
#' @param theme A `ggplot2` theme; defaults to the active theme (see
#'   [mp_set_defaults()], otherwise [mp_theme_publication()]).
#' @param ... Additional arguments passed to the point geom.
#' @param plot.title,plot.subtitle,plot.caption Optional labels.
#' @param font.size Base font size.
#' @param return.data If `TRUE`, return the parsed differential-abundance
#'   data frame (with a `direction` column) instead of a plot; for
#'   `plot.type = "consensus"`, the per-feature consensus table.
#' @param interactive If `TRUE`, return an interactive `plotly` widget instead
#'   of a static `ggplot` (requires `plotly`). See [mp_alpha()].
#' @param verbose Print the detected result type and counts.
#'
#' @return A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
#'   data frame when `return.data = TRUE`.
#' @export
#' @seealso [mp_run_da()] to run the analysis.
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
                                  "heatmap", "consensus", "zeros",
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
                    contrast = NULL,
                    robust.only = TRUE,
                    top.n = 30,
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
  info <- attr(da, "mp_da")
  da <- .mp_da_select(da, contrast)
  # Significance before and after the sensitivity filter: the difference is
  # reported, not silently dropped.
  da$.sig_raw <- .mp_da_direction(da, fc.threshold, p.threshold) != "NS"
  da$direction <- .mp_da_direction(da, fc.threshold, p.threshold,
                                   robust.only)

  if (isTRUE(return.data)) {
    if (plot.type == "consensus") return(.mp_da_consensus_table(da))
    da$.sig_raw <- NULL
    return(da)
  }
  engines <- intersect(unique(da$method), c("ancombc2", "deseq2", "aldex2"))
  da$method <- .mp_method_label(da$method, info$random)
  if (plot.type == "consensus") tab <- .mp_da_consensus_table(da)

  tested <- da[!da$structural_zero, , drop = FALSE]
  n_up <- sum(tested$direction == "Enriched")
  n_dn <- sum(tested$direction == "Depleted")
  n_sz <- length(unique(da$feature_id[da$structural_zero]))
  n_sens <- sum(tested$.sig_raw & tested$robust %in% FALSE)
  if (verbose) {
    rlang::inform(sprintf("mp_diff(): %d enriched, %d depleted (padj < %g, |log2FC| > %g).",
                          n_up, n_dn, p.threshold, fc.threshold))
  }
  n_m <- length(unique(da$method))
  n_c <- length(unique(da$contrast))
  over <- if (n_m > 1 || n_c > 1) {
    sprintf("; summed over %s", .mp_and(c(
      if (n_m > 1) sprintf("%d methods", n_m),
      if (n_c > 1) sprintf("%d contrasts", n_c))))
  } else ""
  subtitle <- plot.subtitle %||% paste(c(
    sprintf("%d enriched \u00b7 %d depleted (padj < %g, |log2FC| > %g%s)",
            n_up, n_dn, p.threshold, fc.threshold, over),
    if (n_sens) sprintf(if (robust.only) "%d excluded as sensitive to pseudo-count" else
      "%d sensitive to pseudo-count", n_sens),
    if (n_sz) sprintf("%d structural zero%s", n_sz, if (n_sz > 1) "s" else "")),
    collapse = " \u00b7 ")

  cols <- c(Enriched = color.up, Depleted = color.down, NS = color.ns)
  dots <- list(...)
  caption <- plot.caption
  p <- switch(
    plot.type,
    volcano = .mp_da_volcano(tested, cols, fc.threshold, p.threshold,
                             y.value, label.top, font.size, dots = dots),
    ma      = .mp_da_ma(tested, cols, label.top, font.size, dots = dots),
    forest  = .mp_da_forest(tested, cols, font.size, top.n),
    dot     = .mp_da_dot(tested, cols, font.size, top.n),
    heatmap = .mp_da_heatmap(da, cols, top.n),
    consensus = .mp_da_consensus(da, tab, cols, top.n),
    zeros   = .mp_da_zeros(da, info, cols),
    concordance = .mp_da_concordance(
      da, result2, method.names, fc.threshold, p.threshold, robust.only,
      feature, log2fc, pvalue, padj, mean, lfcSE, font.size, verbose)
  )
  if (plot.type %in% c("consensus", "zeros", "concordance") &&
      is.null(plot.subtitle)) {
    subtitle <- NULL
  }
  if (is.null(caption)) caption <- attr(p, "mp_caption")

  p <- p + theme +
    ggplot2::labs(title = plot.title, subtitle = subtitle,
                  caption = caption)
  if (plot.type %in% c("heatmap", "consensus", "zeros")) {
    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0),
      panel.grid = ggplot2::element_blank())
  }

  sig <- tested[tested$direction != "NS", , drop = FALSE]
  top <- utils::head(unique(sig$feature[order(sig$padj)]), 3)
  p <- .mp_tag(
    p, "mp_diff",
    params = list(plot.type = plot.type, fc.threshold = fc.threshold,
                  p.threshold = p.threshold, robust.only = robust.only,
                  da = info),
    stats = list(n_up = n_up, n_down = n_dn, n_structural_zero = n_sz,
                 methods = unique(da$method),
                 contrasts = unique(da$contrast),
                 top = as.character(top)),
    notes = info$notes,
    packages = .mp_da_packages(engines, NULL),
    emit = FALSE)
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

# --- Parsing -----------------------------------------------------------------

#' Standard long DA table
#' @noRd
.mp_da_table <- function(feature, log2FC, pvalue = NA, padj = NA, mean = NA,
                         lfcSE = NA, feature_id = feature,
                         contrast = "Comparison", method = "result",
                         robust = NA, structural_zero = FALSE) {
  n <- length(feature)
  out <- data.frame(
    feature = as.character(feature),
    feature_id = as.character(feature_id),
    contrast = rep_len(as.character(contrast), n),
    method = rep_len(as.character(method), n),
    log2FC = as.numeric(log2FC),
    lfcSE = rep_len(as.numeric(lfcSE), n),
    pvalue = rep_len(as.numeric(pvalue), n),
    padj = rep_len(as.numeric(padj), n),
    mean = rep_len(as.numeric(mean), n),
    robust = rep_len(as.logical(robust), n),
    structural_zero = rep_len(as.logical(structural_zero), n),
    stringsAsFactors = FALSE)
  out$structural_zero[is.na(out$structural_zero)] <- FALSE
  out[!is.na(out$log2FC) | out$structural_zero, , drop = FALSE]
}

#' Parse a differential-abundance result into the standard long table
#'
#' @return A data.frame with columns `feature`, `feature_id`, `contrast`,
#'   `method`, `log2FC`, `lfcSE`, `pvalue`, `padj`, `mean`, `robust` and
#'   `structural_zero`, with the run's design in `attr(, "mp_da")` when known.
#' @noRd
.mp_parse_da <- function(result, feature = NULL, log2fc = NULL,
                         pvalue = NULL, padj = NULL, mean = NULL,
                         lfcSE = NULL, verbose = FALSE) {
  cls <- class(result)[1]

  # DESeq2::results() output ------------------------------------------------
  if (inherits(result, "DESeqResults") || cls == "DESeqResults") {
    df <- as.data.frame(result)
    if (verbose) rlang::inform("mp_diff(): detected DESeq2 result.")
    return(.mp_da_table(rownames(df), df$log2FoldChange, df$pvalue, df$padj,
                        df$baseMean, df$lfcSE,
                        contrast = .mp_deseq_contrast(result),
                        method = "deseq2"))
  }

  # ANCOM-BC / ANCOM-BC2 (list with $res) ----------------------------------
  if (is.list(result) && !is.data.frame(result) && "res" %in% names(result)) {
    if (verbose) rlang::inform("mp_diff(): detected ANCOM-BC result.")
    return(.mp_parse_ancombc(result))
  }

  # A named list of results, one per method ---------------------------------
  if (is.list(result) && !is.data.frame(result)) {
    nm <- names(result)
    if (is.null(nm) || any(!nzchar(nm))) {
      rlang::abort(
        c("A list of results must be named, one name per method.",
          i = "e.g. list(DESeq2 = res1, ALDEx2 = res2)."),
        class = "micropubr_da_parse")
    }
    parts <- lapply(seq_along(result), function(i) {
      d <- .mp_parse_da(result[[i]], feature, log2fc, pvalue, padj, mean,
                        lfcSE, verbose = verbose)
      d$method <- rep(nm[i], nrow(d))
      d
    })
    return(do.call(rbind, parts))
  }

  # Generic data frame (incl. mp_run_da() output and ALDEx2) ----------------
  if (is.data.frame(result)) {
    df <- as.data.frame(result)
    info <- attr(result, "mp_da")
    is_aldex <- all(c("diff.btw", "we.ep") %in% names(df))
    feat <- if (!is.null(feature)) df[[feature]] else
      .mp_pick_col(df, c("feature", "taxon", "OTU"), rownames(df))
    # ALDEx2's `effect` is a standardised effect size, not a fold change, so
    # its log2 difference (diff.btw) must win.
    l2fc <- if (!is.null(log2fc)) df[[log2fc]] else
      .mp_pick_col(df, c("log2FC", "log2FoldChange", "diff.btw", "logFC",
                         "coef", "effect"))
    pv <- if (!is.null(pvalue)) df[[pvalue]] else
      .mp_pick_col(df, c("pvalue", "p", "pval", "p.value", "we.ep",
                         "wi.ep", "P.Value"), default = NA_real_)
    pa <- if (!is.null(padj)) df[[padj]] else
      .mp_pick_col(df, c("padj", "p.adj", "FDR", "qvalue", "q", "we.eBH",
                         "wi.eBH", "adj.P.Val"), default = NA_real_)
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
    method <- .mp_pick_col(df, "method",
                           if (is_aldex) "aldex2" else
                             if (!is.null(info)) "mp_run_da" else "result")
    if (verbose) rlang::inform("mp_diff(): parsed a data frame.")
    out <- .mp_da_table(
      feat, l2fc, pv, pa, mn, se,
      feature_id = .mp_pick_col(df, "feature_id", feat),
      contrast = .mp_pick_col(df, "contrast", "Comparison"),
      method = method,
      robust = .mp_pick_col(df, "robust", NA),
      structural_zero = .mp_pick_col(df, "structural_zero", FALSE))
    attr(out, "mp_da") <- info
    return(out)
  }

  rlang::abort(sprintf("Don't know how to parse a <%s> result.", cls),
               class = "micropubr_da_parse")
}

#' "B vs A" from a DESeqResults description, if it is there
#' @noRd
.mp_deseq_contrast <- function(result) {
  desc <- tryCatch(S4Vectors::mcols(result)$description[2],
                   error = function(e) NULL)
  if (is.null(desc) || is.na(desc)) return("Comparison")
  sub("^.*?: *(\\S+ )?", "", desc)
}

#' ANCOM-BC2 (`res` is a data frame) and ANCOM-BC (`res` is a list of
#' data frames), natural-log fold changes converted to log2
#' @noRd
.mp_parse_ancombc <- function(result) {
  res <- result$res
  zeros <- .mp_parse_zero_ind(result$zero_ind)
  if (is.data.frame(res)) {
    method <- "ancombc2"
    terms <- sub("^lfc_", "", grep("^lfc_", names(res), value = TRUE))
    get <- function(prefix, term) {
      nm <- paste0(prefix, term)
      if (nm %in% names(res)) res[[nm]] else rep(NA, nrow(res))
    }
    part <- c(lfc = "lfc_", se = "se_", p = "p_", q = "q_")
    taxa <- res$taxon
  } else if (is.list(res) && is.data.frame(res$lfc)) {
    method <- "ancombc"
    terms <- setdiff(names(res$lfc), "taxon")
    get <- function(prefix, term) {
      d <- res[[prefix]]
      if (!is.null(d) && term %in% names(d)) d[[term]] else
        rep(NA, nrow(res$lfc))
    }
    part <- c(lfc = "lfc", se = "se", p = "p_val", q = "q_val")
    taxa <- res$lfc$taxon %||% rownames(res$lfc)
  } else {
    rlang::abort("Unrecognised ANCOM-BC result: `res` has no fold changes.",
                 class = "micropubr_da_parse")
  }
  # Never plot the intercept; with a known grouping variable keep only its
  # terms, so covariate effects do not masquerade as group contrasts.
  terms <- setdiff(terms, "(Intercept)")
  group_terms <- if (!is.null(zeros)) paste0(zeros$var, zeros$levels) else
    character(0)
  if (any(terms %in% group_terms)) terms <- terms[terms %in% group_terms]
  if (!length(terms)) {
    rlang::abort("The ANCOM-BC result has no non-intercept terms.",
                 class = "micropubr_da_parse")
  }

  # The reference is the one group level without a term of its own.
  ref <- if (!is.null(zeros)) {
    setdiff(zeros$levels, zeros$levels[group_terms %in% terms])
  } else character(0)
  parts <- lapply(terms, function(t) {
    lvl <- if (t %in% group_terms) zeros$levels[match(t, group_terms)] else NA
    label <- if (!is.na(lvl) && length(ref) == 1) sprintf("%s vs %s", lvl, ref) else t
    sz <- rep(FALSE, length(taxa))
    if (!is.na(lvl) && length(ref) == 1) {
      zm <- zeros$matrix[match(taxa, rownames(zeros$matrix)), , drop = FALSE]
      sz <- xor(zm[, lvl], zm[, ref])
      sz[is.na(sz)] <- FALSE
    }
    robust <- if (method == "ancombc2") get("passed_ss_", t) else NA
    lfc <- get(part[["lfc"]], t) / log(2)
    lfc[sz] <- NA
    out <- .mp_da_table(taxa, lfc, get(part[["p"]], t), get(part[["q"]], t),
                        lfcSE = get(part[["se"]], t) / log(2),
                        contrast = label, method = method, robust = robust,
                        structural_zero = sz)
    # ANCOM-BC2 drops structural zeros from `res`; add them back.
    if (method == "ancombc2" && !is.na(lvl) && length(ref) == 1) {
      zm <- zeros$matrix
      ids <- setdiff(rownames(zm)[xor(zm[, lvl], zm[, ref])], taxa)
      if (length(ids)) {
        out <- rbind(out, .mp_da_table(ids, NA, contrast = label,
                                       method = method,
                                       structural_zero = TRUE))
      }
    }
    out
  })
  out <- do.call(rbind, parts)
  if (!is.null(zeros)) {
    attr(out, "mp_da") <- list(zeros = zeros$matrix, method = method)
  }
  out
}

#' ANCOM-BC's zero_ind table -> list(var, levels, matrix taxa x levels)
#' @noRd
.mp_parse_zero_ind <- function(zero_ind) {
  if (is.null(zero_ind) || !is.data.frame(zero_ind) || !nrow(zero_ind)) {
    return(NULL)
  }
  pat <- "^structural_zero \\((.+) = (.+)\\)$"
  cols <- grep(pat, names(zero_ind), value = TRUE)
  if (!length(cols)) return(NULL)
  var <- unique(sub(pat, "\\1", cols))
  if (length(var) != 1) return(NULL)
  levels <- sub(pat, "\\2", cols)
  m <- as.matrix(zero_ind[, cols, drop = FALSE])
  storage.mode(m) <- "logical"
  dimnames(m) <- list(as.character(zero_ind$taxon %||% rownames(zero_ind)),
                      levels)
  list(var = var, levels = levels, matrix = m)
}

#' Pick the first existing column from candidates
#' @noRd
.mp_pick_col <- function(df, candidates, default = NULL) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) df[[hit[1]]] else default
}

#' Display names for method identifiers (other names pass through)
#' @noRd
.mp_method_label <- function(m, random = NULL) {
  lab <- c(wilcox = "Wilcoxon (CLR)",
           lm = if (is.null(random)) "Linear model (CLR)" else
             "Linear mixed model (CLR)",
           ancombc2 = "ANCOM-BC2", ancombc = "ANCOM-BC", deseq2 = "DESeq2",
           aldex2 = "ALDEx2")
  out <- unname(lab[m])
  out[is.na(out)] <- m[is.na(out)]
  out
}

#' Keep the requested contrasts
#' @noRd
.mp_da_select <- function(da, contrast) {
  if (is.null(contrast)) return(da)
  keep <- da$contrast %in% contrast
  if (!any(keep)) {
    rlang::abort(
      c(sprintf("None of the requested contrasts were found: %s.",
                paste(contrast, collapse = ", ")),
        i = sprintf("Available: %s.",
                    paste(unique(da$contrast), collapse = ", "))),
      class = "micropubr_da_no_contrast")
  }
  out <- da[keep, , drop = FALSE]
  attr(out, "mp_da") <- attr(da, "mp_da")
  out
}

#' Classify features as Enriched / Depleted / NS
#' @noRd
.mp_da_direction <- function(da, fc.threshold, p.threshold,
                             robust.only = FALSE) {
  padj <- da$padj
  padj[is.na(padj)] <- 1
  sig <- padj < p.threshold & abs(da$log2FC) > fc.threshold &
    !da$structural_zero
  sig[is.na(sig)] <- FALSE
  if (robust.only) sig <- sig & !(da$robust %in% FALSE)
  out <- rep("NS", nrow(da))
  out[sig & da$log2FC > 0] <- "Enriched"
  out[sig & da$log2FC < 0] <- "Depleted"
  factor(out, levels = c("Enriched", "Depleted", "NS"))
}

#' Panel label: contrast, method, or both
#' @noRd
.mp_da_panels <- function(da) {
  multi_c <- length(unique(da$contrast)) > 1
  multi_m <- length(unique(da$method)) > 1
  if (multi_c && multi_m) paste(da$method, da$contrast, sep = " \u00b7 ") else
    if (multi_m) da$method else if (multi_c) da$contrast else NULL
}

#' Facet a DA plot by contrast and/or method
#' @noRd
.mp_da_facet <- function(p, da, free_y = FALSE) {
  multi_c <- length(unique(da$contrast)) > 1
  multi_m <- length(unique(da$method)) > 1
  if (!multi_c && !multi_m) return(p)
  scales <- if (free_y) "free_y" else "fixed"
  if (multi_c && multi_m) {
    p + ggplot2::facet_grid(rows = ggplot2::vars(.data$method),
                            cols = ggplot2::vars(.data$contrast),
                            scales = scales)
  } else {
    v <- if (multi_m) "method" else "contrast"
    p + ggplot2::facet_wrap(ggplot2::vars(.data[[v]]), scales = scales)
  }
}

#' Pick the `label.top` most significant features (per panel) for labelling
#' @noRd
.mp_da_top <- function(da, label.top) {
  sig <- da[da$direction != "NS", , drop = FALSE]
  if (nrow(sig) == 0 || label.top <= 0) return(sig[0, ])
  panel <- .mp_da_panels(sig) %||% rep("all", nrow(sig))
  pieces <- lapply(split(sig, panel), function(s) {
    key <- s$padj
    key[is.na(key)] <- 1
    s[order(key, -abs(s$log2FC))[seq_len(min(label.top, nrow(s)))], ,
      drop = FALSE]
  })
  do.call(rbind, pieces)
}

#' Shape scale distinguishing results that fail the sensitivity analysis
#' @noRd
.mp_da_robust_shape <- function(da) {
  lab <- ifelse(da$direction != "NS" & da$robust %in% FALSE,
                "Sensitive to pseudo-count", "Robust")
  factor(lab, levels = c("Robust", "Sensitive to pseudo-count"))
}

# --- Plots -------------------------------------------------------------------

#' @noRd
.mp_da_volcano <- function(da, cols, fc.threshold, p.threshold, y.value,
                           label.top, font.size, dots = list()) {
  da$y <- -log10(pmax(da[[y.value]], .Machine$double.xmin))
  da <- da[!is.na(da$y), , drop = FALSE]
  ylab <- if (y.value == "padj") expression(-log[10]~adjusted~italic(p)) else
    expression(-log[10]~italic(p))
  mapping <- ggplot2::aes(x = .data$log2FC, y = .data$y,
                          colour = .data$direction)
  has_robust <- any(da$direction != "NS" & da$robust %in% FALSE)
  if (has_robust) {
    da$.robust <- .mp_da_robust_shape(da)
    mapping$shape <- ggplot2::aes(shape = .data$.robust)$shape
  }
  p <- ggplot2::ggplot(da, mapping) +
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
  if (has_robust) {
    p <- p + ggplot2::scale_shape_manual(
      values = c(Robust = 16, `Sensitive to pseudo-count` = 1), name = NULL)
  }
  p <- .mp_da_facet(p, da)
  .mp_add_da_labels(p, da, label.top, font.size,
                    ggplot2::aes(label = .data$feature))
}

#' @noRd
.mp_da_ma <- function(da, cols, label.top, font.size, dots = list()) {
  if (all(is.na(da$mean))) {
    rlang::abort("MA plot requires a mean-abundance column.",
                 class = "micropubr_da_no_mean")
  }
  # ALDEx2 reports mean CLR values (which can be negative) rather than
  # abundances; plot those on their own scale.
  da <- da[!is.na(da$mean), , drop = FALSE]
  clr_scale <- any(da$mean < 0, na.rm = TRUE)
  da$x <- if (clr_scale) da$mean else log10(da$mean + 1)
  xlab <- if (clr_scale) "Mean CLR abundance" else
    expression(log[10]~mean~abundance)
  p <- ggplot2::ggplot(da, ggplot2::aes(x = .data$x, y = .data$log2FC,
                                        colour = .data$direction)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    .mp_geom(ggplot2::geom_point,
             params = list(alpha = 0.8, size = 1.8), dots = dots) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
    ggplot2::labs(x = xlab, y = expression(log[2]~fold-change))
  p <- .mp_da_facet(p, da)
  .mp_add_da_labels(p, da, label.top, font.size,
                    ggplot2::aes(label = .data$feature))
}

#' The `top.n` most significant features per panel, for forest and dot plots
#'
#' A forest of 140 taxa is unreadable, so only the strongest are drawn; the
#' returned caption says how many of how many are shown.
#' @noRd
.mp_da_top_sig <- function(da, top.n, what) {
  sig <- da[da$direction != "NS", , drop = FALSE]
  if (nrow(sig) == 0) {
    rlang::abort(sprintf("No significant features to display in the %s plot.",
                         what),
                 class = "micropubr_da_empty")
  }
  shown <- .mp_da_top(sig, top.n)
  caption <- if (nrow(shown) < nrow(sig)) {
    sprintf("Showing the %d most significant of %d significant features%s.",
            top.n, nrow(sig),
            if (!is.null(.mp_da_panels(sig))) " per panel" else "")
  }
  list(data = shown, caption = caption)
}

#' @noRd
.mp_da_forest <- function(da, cols, font.size, top.n = Inf) {
  top <- .mp_da_top_sig(da, top.n, "forest")
  sig <- top$data
  sig$lfcSE[is.na(sig$lfcSE)] <- 0
  sig$feature <- stats::reorder(sig$feature, sig$log2FC)
  mapping <- ggplot2::aes(x = .data$log2FC, y = .data$feature,
                          colour = .data$direction)
  has_robust <- any(sig$robust %in% FALSE)
  if (has_robust) {
    sig$.robust <- .mp_da_robust_shape(sig)
    mapping$shape <- ggplot2::aes(shape = .data$.robust)$shape
  }
  p <- ggplot2::ggplot(sig, mapping) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$log2FC - 1.96 * .data$lfcSE,
                   xmax = .data$log2FC + 1.96 * .data$lfcSE),
      orientation = "y", width = 0.25, linewidth = 0.4) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_colour_manual(values = cols, name = NULL,
                                 breaks = c("Enriched", "Depleted")) +
    ggplot2::labs(x = expression(log[2]~fold-change~(95*"%"~CI)), y = NULL)
  if (has_robust) {
    p <- p + ggplot2::scale_shape_manual(
      values = c(Robust = 16, `Sensitive to pseudo-count` = 1), name = NULL)
  }
  p <- .mp_da_facet(p, sig, free_y = TRUE)
  attr(p, "mp_caption") <- top$caption
  p
}

#' @noRd
.mp_da_dot <- function(da, cols, font.size, top.n = Inf) {
  top <- .mp_da_top_sig(da, top.n, "dot")
  sig <- top$data
  sig$feature <- stats::reorder(sig$feature, sig$log2FC)
  sig$absp <- -log10(pmax(sig$padj, .Machine$double.xmin))
  p <- ggplot2::ggplot(sig, ggplot2::aes(x = .data$log2FC, y = .data$feature,
                                         colour = .data$direction,
                                         size = .data$absp)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_colour_manual(values = cols, name = NULL,
                                 breaks = c("Enriched", "Depleted")) +
    ggplot2::scale_size_area(max_size = 6,
                             name = expression(-log[10]~padj)) +
    ggplot2::labs(x = expression(log[2]~fold-change), y = NULL)
  p <- .mp_da_facet(p, sig, free_y = TRUE)
  attr(p, "mp_caption") <- top$caption
  p
}

#' Fold changes of significant features across contrasts (x) per method
#' @noRd
.mp_da_heatmap <- function(da, cols, top.n) {
  sig_ids <- unique(da$feature_id[da$direction != "NS"])
  if (!length(sig_ids)) {
    rlang::abort("No significant features to display in the heatmap.",
                 class = "micropubr_da_empty")
  }
  # Keep the `top.n` features with the smallest adjusted p in any contrast.
  best <- tapply(ifelse(is.na(da$padj), 1, da$padj), da$feature_id, min)
  keep <- utils::head(names(sort(best[sig_ids])), top.n)
  d <- da[da$feature_id %in% keep, , drop = FALSE]

  # Order rows by clustering the fold-change profile (contrasts x methods).
  d$.col <- paste(d$method, d$contrast, sep = "\r")
  m <- tapply(d$log2FC, list(d$feature_id, d$.col), mean)
  m[is.na(m)] <- 0
  ord <- if (nrow(m) > 2 && ncol(m) > 1) {
    rownames(m)[stats::hclust(stats::dist(m))$order]
  } else {
    rownames(m)[order(rowMeans(m))]
  }
  labs <- d$feature[match(ord, d$feature_id)]
  d$feature <- factor(d$feature, levels = unique(labs))
  d$contrast <- factor(d$contrast, levels = unique(da$contrast))

  d$.mark <- ifelse(d$structural_zero, "\u00d8",
                    ifelse(d$.sig_raw & d$robust %in% FALSE, "\u2218",
                           ifelse(d$direction == "NS", "", "*")))
  lim <- max(abs(d$log2FC), na.rm = TRUE)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$contrast, y = .data$feature,
                                       fill = .data$log2FC)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$.mark), size = 4.5,
                       colour = "grey10", vjust = 0.75) +
    ggplot2::scale_fill_gradient2(low = cols[["Depleted"]],
                                  mid = "white", high = cols[["Enriched"]],
                                  midpoint = 0, limits = c(-lim, lim),
                                  na.value = "grey90",
                                  name = expression(log[2]~FC)) +
    ggplot2::labs(x = NULL, y = NULL)
  if (length(unique(d$method)) > 1) {
    p <- p + ggplot2::facet_grid(cols = ggplot2::vars(.data$method))
  }
  marks <- c("* adjusted p below threshold",
             if (any(d$.mark == "\u2218")) "\u2218 sensitive to pseudo-count (ANCOM-BC2 sensitivity analysis)",
             if (any(d$.mark == "\u00d8")) "\u00d8 structural zero (absent from one group)")
  attr(p, "mp_caption") <- paste(marks, collapse = "; ")
  p
}

#' Per-feature agreement across methods
#' @noRd
.mp_da_consensus_table <- function(da) {
  methods <- unique(da$method)
  if (length(methods) < 2) {
    rlang::abort(
      c("`plot.type = \"consensus\"` needs results from at least two methods.",
        i = "Use mp_run_da(method = c(\"wilcox\", \"lm\", \"ancombc2\")) or pass a named list of results."),
      class = "micropubr_da_one_method")
  }
  d <- da[!da$structural_zero, , drop = FALSE]
  key <- paste(d$feature_id, d$contrast, sep = "\r")
  pieces <- lapply(split(d, key), function(s) {
    up <- sum(s$direction == "Enriched")
    dn <- sum(s$direction == "Depleted")
    data.frame(
      feature = s$feature[1], feature_id = s$feature_id[1],
      contrast = s$contrast[1], n_methods = length(unique(s$method)),
      n_significant = up + dn, n_enriched = up, n_depleted = dn,
      mean_log2FC = mean(s$log2FC, na.rm = TRUE),
      consensus = if (up && dn) "Conflicting" else if (up) "Enriched" else
        if (dn) "Depleted" else "NS",
      stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, pieces)
  tab$agreement <- tab$n_significant / length(methods)
  tab <- tab[order(-tab$n_significant, -abs(tab$mean_log2FC)), ,
             drop = FALSE]
  rownames(tab) <- NULL
  attr(tab, "methods") <- methods
  tab
}

#' Feature x method matrix of calls, with the agreement count alongside
#' @noRd
.mp_da_consensus <- function(da, tab, cols, top.n) {
  methods <- attr(tab, "methods")
  m <- length(methods)
  sig <- tab[tab$n_significant > 0, , drop = FALSE]
  if (!nrow(sig)) {
    rlang::abort("No feature is significant in any method.",
                 class = "micropubr_da_empty")
  }
  sig <- utils::head(sig, top.n)
  key <- paste(sig$feature_id, sig$contrast, sep = "\r")
  # A full feature x method grid, so a method that did not test a feature
  # (filtered, or an ANCOM-BC2 structural zero) shows as such, not as a gap.
  called <- da[!da$structural_zero, , drop = FALSE]
  called$.key <- paste(called$feature_id, called$contrast, sep = "\r")
  d <- expand.grid(.key = key, method = methods, stringsAsFactors = FALSE)
  hit <- match(paste(d$.key, d$method), paste(called$.key, called$method))
  d$direction <- as.character(called$direction[hit])
  d$direction[is.na(d$direction)] <- "Not tested"
  d$direction <- factor(d$direction,
                        levels = c("Enriched", "Depleted", "NS", "Not tested"))
  d$feature <- sig$feature[match(d$.key, key)]
  d$contrast <- sig$contrast[match(d$.key, key)]
  lev <- unique(sig$feature[order(sig$n_significant, abs(sig$mean_log2FC))])
  d$feature <- factor(d$feature, levels = lev)
  d$method <- factor(d$method, levels = c(methods, "Agreement"))
  sig$feature <- factor(sig$feature, levels = lev)
  sig$method <- factor("Agreement", levels = levels(d$method))
  sig$.lab <- sprintf("%d/%d", sig$n_significant, m)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$method, y = .data$feature)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$direction),
                       colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(data = sig, ggplot2::aes(label = .data$.lab),
                       size = 3.2, colour = "grey15") +
    ggplot2::scale_fill_manual(values = c(cols, `Not tested` = "grey97"),
                               name = NULL, drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL)
  if (length(unique(d$contrast)) > 1) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data$contrast),
                                 scales = "free_y", space = "free_y")
  }

  all_m <- sum(tab$n_significant == m)
  half <- sum(tab$n_significant >= m / 2 & tab$n_significant < m)
  one <- sum(tab$n_significant == 1)
  conflict <- sum(tab$consensus == "Conflicting")
  rho <- .mp_da_method_correlation(da, methods)
  attr(p, "mp_caption") <- paste(c(
    sprintf("%d feature%s significant in all %d methods, %d in at least half, %d in only one.",
            all_m, if (all_m == 1) "" else "s", m, half, one),
    if (conflict) sprintf("%d with opposite directions between methods.", conflict),
    if (!is.null(rho)) sprintf("Spearman r of log2FC between methods: %s.", rho)),
    collapse = " ")
  p
}

#' Range of pairwise Spearman correlations of fold changes between methods
#' @noRd
.mp_da_method_correlation <- function(da, methods) {
  d <- da[!da$structural_zero, , drop = FALSE]
  d$.key <- paste(d$feature_id, d$contrast, sep = "\r")
  w <- lapply(methods, function(m) {
    s <- d[d$method == m, , drop = FALSE]
    stats::setNames(s$log2FC, s$.key)
  })
  pairs <- utils::combn(length(methods), 2)
  r <- apply(pairs, 2, function(ij) {
    k <- intersect(names(w[[ij[1]]]), names(w[[ij[2]]]))
    if (length(k) < 3) return(NA_real_)
    suppressWarnings(stats::cor(w[[ij[1]]][k], w[[ij[2]]][k],
                                method = "spearman"))
  })
  r <- r[!is.na(r)]
  if (!length(r)) return(NULL)
  if (length(r) == 1) sprintf("%.2f", r) else
    sprintf("%.2f to %.2f", min(r), max(r))
}

#' Presence/absence of structurally zero taxa across groups
#' @noRd
.mp_da_zeros <- function(da, info, cols) {
  zm <- info$zeros
  if (is.null(zm)) {
    rlang::abort(
      c("No structural-zero information in this result.",
        i = "Run mp_run_da(method = \"ancombc2\") or pass an ANCOM-BC2 result with struc_zero = TRUE."),
      class = "micropubr_da_no_zeros")
  }
  zm <- zm[rowSums(zm) > 0 & rowSums(zm) < ncol(zm), , drop = FALSE]
  if (!nrow(zm)) {
    rlang::abort("No structural zeros were detected.",
                 class = "micropubr_da_empty")
  }
  labels <- info$labels %||% stats::setNames(da$feature, da$feature_id)
  feat <- labels[rownames(zm)]
  feat[is.na(feat)] <- rownames(zm)[is.na(feat)]
  d <- data.frame(
    feature = rep(unname(feat), ncol(zm)),
    group = rep(colnames(zm), each = nrow(zm)),
    status = ifelse(as.vector(zm), "Absent", "Present"),
    stringsAsFactors = FALSE)
  d$group <- factor(d$group, levels = colnames(zm))
  d$feature <- factor(d$feature, levels = unique(feat[order(rowSums(zm))]))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$group, y = .data$feature,
                                       fill = .data$status)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c(Present = cols[["Enriched"]], Absent = "grey92"),
      name = NULL) +
    ggplot2::labs(x = NULL, y = NULL)
  attr(p, "mp_caption") <- sprintf(
    "%d taxa with structural zeros (ANCOM-BC2): systematically absent from at least one group, so not tested for differential abundance.",
    nrow(zm))
  p
}

#' @noRd
.mp_da_concordance <- function(da1, result2, method.names, fc.threshold,
                               p.threshold, robust.only, feature, log2fc,
                               pvalue, padj, mean, lfcSE, font.size,
                               verbose) {
  if (is.null(result2)) {
    methods <- unique(da1$method)
    if (length(methods) != 2) {
      rlang::abort(
        "`plot.type = \"concordance\"` requires `result2`, or a result holding exactly two methods.",
        class = "micropubr_da_no_result2")
    }
    da2 <- da1[da1$method == methods[2], , drop = FALSE]
    da1 <- da1[da1$method == methods[1], , drop = FALSE]
    if (identical(method.names, c("Method 1", "Method 2"))) {
      method.names <- methods
    }
  } else {
    da2 <- .mp_parse_da(result2, feature, log2fc, pvalue, padj, mean, lfcSE,
                        verbose = verbose)
    da2$direction <- .mp_da_direction(da2, fc.threshold, p.threshold,
                                      robust.only)
  }
  da1 <- da1[!da1$structural_zero, , drop = FALSE]
  da2 <- da2[!da2$structural_zero, , drop = FALSE]
  key <- function(d) {
    if (length(unique(d$contrast)) > 1) paste(d$feature_id, d$contrast) else
      d$feature_id
  }
  m <- merge(
    data.frame(key = key(da1), fc1 = da1$log2FC, d1 = da1$direction),
    data.frame(key = key(da2), fc2 = da2$log2FC, d2 = da2$direction),
    by = "key"
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
  p <- ggplot2::ggplot(m, ggplot2::aes(x = .data$fc1, y = .data$fc2,
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
      y = sprintf("log2FC \u2014 %s", method.names[2]))
  attr(p, "mp_caption") <- sprintf(
    "Spearman r = %.2f \u00b7 %d features significant in both", rho, n_both)
  p
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
