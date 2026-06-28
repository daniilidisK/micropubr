# Microbiome statistical analyses ----------------------------------------
# Each public function runs a standard microbial-ecology test and returns a
# publication-ready ggplot with the result baked into the caption, following
# the same conventions as the other mp_*() functions. mp_run_da() is the
# exception: it *runs* a differential-abundance method and returns a tidy
# table for mp_diff().

# --- Beta-dispersion and group-difference tests -------------------------

#' Multivariate dispersion (PERMDISP) with optional group tests
#'
#' Runs [vegan::betadisper()] to test homogeneity of multivariate dispersions
#' (PERMDISP) -- the standard companion to PERMANOVA -- and plots each sample's
#' distance to its group centroid. Additional group-difference statistics
#' (ANOSIM, MRPP, PERMANOVA) can be appended to the caption.
#'
#' @param object A [phyloseq::phyloseq] object.
#' @param group.by Sample-metadata column defining the groups.
#' @param distance Distance metric (see [mp_ordination()]; default `"bray"`).
#' @param add.test Optional character vector of further tests to summarise in
#'   the caption: any of `"anosim"`, `"mrpp"`, `"permanova"`.
#' @param permutations Number of permutations.
#' @param plot.type One of `"box"` or `"violin"`.
#' @param point.size Size of jittered sample points.
#' @param palette,theme,plot.title,plot.subtitle,plot.caption,font.size See
#'   [mp_alpha()].
#' @param ... Passed to the primary geom.
#' @param interactive If `TRUE`, return an interactive `plotly` widget instead
#'   of a static `ggplot` (requires `plotly`). See [mp_alpha()].
#' @param return.data If `TRUE`, return the per-sample distance-to-centroid
#'   data frame (with the betadisper object attached) instead of a plot.
#' @param verbose Print the methods used.
#'
#' @return A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
#'   data frame when `return.data = TRUE`.
#' @export
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' mp_dispersion(GlobalPatterns, group.by = "SampleType",
#'               add.test = "permanova")
#' }
mp_dispersion <- function(object,
                          group.by = NULL,
                          distance = "bray",
                          add.test = NULL,
                          permutations = 999,
                          plot.type = c("box", "violin"),
                          point.size = 1.5,
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
  .mp_require("vegan", "Multivariate dispersion (PERMDISP)")
  .mp_check_phyloseq(object)
  plot.type <- match.arg(plot.type)
  group.by <- .mp_check_var(object, group.by, "group.by", required = TRUE)
  theme <- .mp_theme(theme, font.size)
  if (!is.null(add.test)) {
    add.test <- match.arg(add.test, c("anosim", "mrpp", "permanova"),
                          several.ok = TRUE)
  }

  d <- .mp_compute_distance(object, distance)
  meta <- .mp_sample_data(object)[labels(d), , drop = FALSE]
  g <- factor(meta[[group.by]])
  keep <- !is.na(g)
  if (any(!keep)) {
    d <- stats::as.dist(as.matrix(d)[keep, keep])
    g <- droplevels(g[keep])
  }
  if (nlevels(g) < 2) {
    rlang::abort("`group.by` must have at least two non-missing groups.",
                 class = "micropubr_one_group")
  }

  bd <- vegan::betadisper(d, g)
  pt <- vegan::permutest(bd, permutations = permutations)
  fval <- pt$tab[["F"]][1]
  pval <- pt$tab[["Pr(>F)"]][1]

  df <- data.frame(.sample = names(bd$distances),
                   distance = as.numeric(bd$distances),
                   group = g, stringsAsFactors = FALSE)
  names(df)[names(df) == "group"] <- group.by

  if (isTRUE(return.data)) {
    attr(df, "betadisper") <- bd
    return(df)
  }

  caption <- plot.caption %||% sprintf(
    "PERMDISP: F = %.2f, p = %s (%d permutations)",
    fval, format.pval(pval, digits = 2), permutations)
  caption <- paste(c(caption,
                     .mp_extra_beta_tests(d, meta[keep, , drop = FALSE],
                                          group.by, g, add.test,
                                          permutations)),
                   collapse = "\n")
  if (verbose) rlang::inform(caption)

  ng <- nlevels(g)
  cols <- .mp_resolve_palette(palette, ng, type = "qualitative")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[group.by]],
                                        y = .data$distance,
                                        fill = .data[[group.by]]))
  if (plot.type == "violin") {
    p <- p + .mp_geom(ggplot2::geom_violin,
                      params = list(trim = FALSE, alpha = 0.6, colour = NA),
                      dots = list(...)) +
      ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white")
  } else {
    p <- p + .mp_geom(ggplot2::geom_boxplot,
                      params = list(outlier.shape = NA, alpha = 0.7),
                      dots = list(...))
  }
  p <- p +
    ggplot2::geom_jitter(width = 0.18, height = 0, size = point.size,
                         alpha = 0.6, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = cols, name = group.by) +
    ggplot2::labs(x = NULL, y = "Distance to group centroid",
                  title = plot.title, subtitle = plot.subtitle,
                  caption = caption) +
    theme
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

#' Caption fragments for ANOSIM / MRPP / PERMANOVA
#' @noRd
.mp_extra_beta_tests <- function(d, meta, group.by, g, add.test,
                                 permutations) {
  out <- character(0)
  if ("anosim" %in% add.test) {
    an <- tryCatch(vegan::anosim(d, g, permutations = permutations),
                   error = function(e) NULL)
    if (!is.null(an)) {
      out <- c(out, sprintf("ANOSIM: R = %.3f, p = %s",
                            an$statistic, format.pval(an$signif, digits = 2)))
    }
  }
  if ("mrpp" %in% add.test) {
    mr <- tryCatch(vegan::mrpp(d, g, permutations = permutations),
                   error = function(e) NULL)
    if (!is.null(mr)) {
      out <- c(out, sprintf("MRPP: A = %.3f, p = %s",
                            mr$A, format.pval(mr$Pvalue, digits = 2)))
    }
  }
  if ("permanova" %in% add.test) {
    perm <- .mp_permanova(d, meta, group.by, permutations)
    if (!is.null(perm)) out <- c(out, perm$caption)
  }
  out
}

# --- Mantel test / correlogram ------------------------------------------

#' Mantel test between community and explanatory distances
#'
#' Correlates the community dissimilarity matrix with a second distance matrix
#' -- typically an environmental/spatial distance derived from numeric metadata
#' -- via [vegan::mantel()], and visualises the relationship as a
#' distance-distance scatter or a Mantel correlogram.
#'
#' @param object A [phyloseq::phyloseq] object.
#' @param distance Community distance metric (default `"bray"`).
#' @param vars Numeric metadata column(s) used to build the second
#'   (Euclidean, standardised) distance matrix. Ignored if `dist2` is given.
#' @param dist2 An explicit second distance (a `dist` object or square matrix).
#' @param method Correlation method: `"pearson"` or `"spearman"`.
#' @param permutations Number of permutations.
#' @param plot.type `"scatter"` or `"correlogram"`.
#' @param point.size Point size for the scatter plot.
#' @param palette,theme,plot.title,plot.subtitle,plot.caption,font.size See
#'   [mp_alpha()].
#' @param ... Passed to the primary geom.
#' @param interactive If `TRUE`, return an interactive `plotly` widget instead
#'   of a static `ggplot` (requires `plotly`). See [mp_alpha()].
#' @param return.data If `TRUE`, return the underlying data frame instead of a
#'   plot (pairwise distances for `"scatter"`, the correlogram table for
#'   `"correlogram"`).
#' @param verbose Print the methods used.
#'
#' @return A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
#'   data frame when `return.data = TRUE`.
#' @export
mp_mantel <- function(object,
                      distance = "bray",
                      vars = NULL,
                      dist2 = NULL,
                      method = c("pearson", "spearman"),
                      permutations = 999,
                      plot.type = c("scatter", "correlogram"),
                      point.size = 1.2,
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
  .mp_require("vegan", "Mantel test")
  .mp_check_phyloseq(object)
  method <- match.arg(method)
  plot.type <- match.arg(plot.type)
  theme <- .mp_theme(theme, font.size)

  d1 <- .mp_compute_distance(object, distance)
  labs <- labels(d1)
  meta <- .mp_sample_data(object)[labs, , drop = FALSE]

  if (is.null(dist2)) {
    if (is.null(vars)) {
      rlang::abort(
        c("Supply `vars` (numeric metadata) or an explicit `dist2`.",
          i = "e.g. mp_mantel(ps, vars = c(\"pH\", \"Temperature\"))"),
        class = "micropubr_mantel_input")
    }
    vars <- .mp_check_var(object, vars, "vars")
    num <- vars[vapply(vars, function(v) is.numeric(meta[[v]]), logical(1))]
    if (length(num) == 0) {
      rlang::abort("`vars` must include at least one numeric column.",
                   class = "micropubr_mantel_nonnumeric")
    }
    env <- scale(as.matrix(meta[, num, drop = FALSE]))
    d2 <- stats::dist(env, method = "euclidean")
    d2lab <- sprintf("Euclidean distance (%s)", paste(num, collapse = ", "))
  } else {
    d2 <- stats::as.dist(as.matrix(dist2))
    d2lab <- "Second distance"
  }

  mant <- vegan::mantel(d1, d2, method = method, permutations = permutations)
  rr <- mant$statistic
  pp <- mant$signif
  caption <- plot.caption %||% sprintf(
    "Mantel (%s): r = %.3f, p = %s (%d permutations)",
    method, rr, format.pval(pp, digits = 2), permutations)
  if (verbose) rlang::inform(caption)

  if (plot.type == "correlogram") {
    mc <- vegan::mantel.correlog(d1, D.geo = d2, nperm = permutations)
    res <- as.data.frame(mc$mantel.res)
    res <- res[stats::complete.cases(res[, "Mantel.cor"]), , drop = FALSE]
    res$.sig <- !is.na(res[["Pr(corrected)"]]) &
      res[["Pr(corrected)"]] < 0.05
    if (isTRUE(return.data)) return(res)
    p <- ggplot2::ggplot(res, ggplot2::aes(x = .data[["class.index"]],
                                           y = .data[["Mantel.cor"]])) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60",
                          linewidth = 0.3) +
      ggplot2::geom_line(colour = "grey40", linewidth = 0.4) +
      .mp_geom(ggplot2::geom_point,
               params = list(mapping = ggplot2::aes(fill = .data$.sig),
                             size = 3, shape = 21, colour = "grey20"),
               dots = list(...)) +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = "#B2182B", "FALSE" = "white"),
        name = "p < 0.05", labels = c("TRUE" = "yes", "FALSE" = "no")) +
      ggplot2::labs(x = "Distance class", y = "Mantel correlation",
                    title = plot.title, subtitle = plot.subtitle,
                    caption = caption) +
      theme
    if (isTRUE(interactive)) return(.mp_interactive(p))
    return(p)
  }

  df <- data.frame(community = as.numeric(d1), explanatory = as.numeric(d2))
  if (isTRUE(return.data)) return(df)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$explanatory,
                                        y = .data$community)) +
    .mp_geom(ggplot2::geom_point,
             params = list(alpha = 0.4, size = point.size,
                           colour = "#0072B2"),
             dots = list(...)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                         colour = "grey20", linewidth = 0.5) +
    ggplot2::labs(x = d2lab,
                  y = sprintf("%s dissimilarity", distance),
                  title = plot.title, subtitle = plot.subtitle,
                  caption = caption) +
    theme
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

# --- Variance partitioning ----------------------------------------------

#' Variance partitioning of community structure
#'
#' Partitions the variation in community composition (Hellinger-transformed)
#' among two or three sets of explanatory metadata variables using
#' [vegan::varpart()], and plots the adjusted-R-squared of each fraction.
#'
#' @param object A [phyloseq::phyloseq] object.
#' @param groups A **named** list mapping each explanatory set to a character
#'   vector of metadata columns, e.g.
#'   `list(Host = c("Age", "Sex"), Diet = "Treatment")`. Length 2 or 3.
#' @param transform Community standardisation passed to
#'   [vegan::decostand()] (default `"hellinger"`).
#' @param palette,theme,plot.title,plot.subtitle,plot.caption,font.size See
#'   [mp_alpha()].
#' @param ... Passed to the primary geom.
#' @param interactive If `TRUE`, return an interactive `plotly` widget instead
#'   of a static `ggplot` (requires `plotly`). See [mp_alpha()].
#' @param return.data If `TRUE`, return the fractions table instead of a plot.
#' @param verbose Print the methods used.
#'
#' @return A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
#'   data frame when `return.data = TRUE`.
#' @export
mp_varpart <- function(object,
                       groups,
                       transform = "hellinger",
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
  .mp_require("vegan", "Variance partitioning")
  .mp_check_phyloseq(object)
  if (!is.list(groups) || is.null(names(groups)) ||
      length(groups) < 2 || length(groups) > 3) {
    rlang::abort(
      c("`groups` must be a named list of length 2 or 3.",
        i = "e.g. list(Host = c(\"Age\",\"Sex\"), Diet = \"Treatment\")"),
      class = "micropubr_varpart_groups")
  }
  for (g in groups) .mp_check_var(object, g, "groups")
  theme <- .mp_theme(theme, font.size)

  comm <- .mp_otu_matrix(object, samples_as_rows = TRUE)
  Y <- vegan::decostand(comm, method = transform)
  meta <- .mp_sample_data(object)[rownames(comm), , drop = FALSE]
  forms <- lapply(groups, function(v)
    stats::reformulate(sprintf("`%s`", v)))

  # varpart() takes the explanatory tables positionally (X, ...); the formulas
  # must be unnamed or they are swallowed as named arguments and X is empty.
  vp <- do.call(vegan::varpart,
                c(list(Y), unname(forms), list(data = meta)))
  fr <- vp$part$indfract
  df <- data.frame(fraction = rownames(fr),
                   adj.r2 = fr$Adj.R.square,
                   stringsAsFactors = FALSE)
  df$adj.r2[is.na(df$adj.r2)] <- 0
  if (verbose) {
    rlang::inform(sprintf("mp_varpart(): partitioned among %s.",
                          paste(names(groups), collapse = ", ")))
  }
  if (isTRUE(return.data)) {
    attr(df, "varpart") <- vp
    return(df)
  }

  cols <- .mp_resolve_palette(palette, nrow(df), type = "qualitative")
  df$fraction <- factor(df$fraction, levels = df$fraction)
  caption <- plot.caption %||% sprintf(
    "Variance partitioning of community structure among %s.",
    paste(names(groups), collapse = ", "))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$fraction, y = .data$adj.r2,
                                        fill = .data$fraction)) +
    .mp_geom(ggplot2::geom_col, params = list(width = 0.8), dots = list(...)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60",
                        linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::labs(x = NULL, y = "Adjusted R\u00b2",
                  title = plot.title,
                  subtitle = plot.subtitle %||%
                    paste(sprintf("%s: %s", names(groups),
                                  vapply(groups, paste, character(1),
                                         collapse = "+")),
                          collapse = "   "),
                  caption = caption) +
    theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

# --- Prevalence / core microbiome ---------------------------------------

#' Prevalence-abundance and core-microbiome plots
#'
#' Summarises how widely each taxon is detected. `plot.type = "scatter"` plots
#' prevalence against mean relative abundance (a standard QC / core-selection
#' view); `plot.type = "core"` draws a core-microbiome heatmap of prevalence
#' across a gradient of detection thresholds.
#'
#' @param object A [phyloseq::phyloseq] object (counts).
#' @param tax.level Optional rank to agglomerate to first.
#' @param plot.type `"scatter"` or `"core"`.
#' @param detection Detection threshold (relative abundance) defining presence
#'   for the scatter plot's prevalence and the core threshold line.
#' @param prevalence Prevalence threshold (fraction of samples) for the core
#'   line / core membership.
#' @param color.by For `"scatter"`, a taxonomy rank to colour points by.
#' @param label.top Number of most-prevalent taxa to label (scatter).
#' @param palette,theme,plot.title,plot.subtitle,plot.caption,font.size See
#'   [mp_alpha()].
#' @param ... Passed to the primary geom.
#' @param interactive If `TRUE`, return an interactive `plotly` widget instead
#'   of a static `ggplot` (requires `plotly`). See [mp_alpha()].
#' @param return.data If `TRUE`, return the prevalence/abundance table.
#' @param verbose Print the methods used.
#'
#' @return A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
#'   data frame when `return.data = TRUE`.
#' @export
mp_prevalence <- function(object,
                          tax.level = NULL,
                          plot.type = c("scatter", "core"),
                          detection = 0.001,
                          prevalence = 0.5,
                          color.by = NULL,
                          label.top = 10,
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
  theme <- .mp_theme(theme, font.size)

  obj <- .mp_agglomerate(object, tax.level)
  rel <- .mp_to_relative(obj)
  mat <- .mp_otu_matrix(rel, samples_as_rows = TRUE)   # samples x taxa
  nsmp <- nrow(mat)

  prev <- colSums(mat > detection) / nsmp
  meanab <- colMeans(mat)
  tab <- data.frame(feature = colnames(mat), prevalence = prev,
                    mean_abundance = meanab, stringsAsFactors = FALSE)
  tab <- .mp_attach_taxonomy(tab, obj, tax.level, color.by)

  if (plot.type == "core") {
    core_df <- .mp_core_grid(mat, prevalence)
    if (isTRUE(return.data)) return(core_df)
    ph <- .mp_core_heatmap(core_df, palette, theme, plot.title,
                           plot.subtitle, plot.caption, prevalence,
                           list(...))
    if (isTRUE(interactive)) return(.mp_interactive(ph))
    return(ph)
  }

  if (isTRUE(return.data)) return(tab)

  tab$is_core <- tab$prevalence >= prevalence & tab$mean_abundance >= detection
  ncore <- sum(tab$is_core)
  caption <- plot.caption %||% sprintf(
    "%d core taxa (prevalence \u2265 %.0f%%, abundance \u2265 %s).",
    ncore, 100 * prevalence, format(detection, scientific = FALSE))

  colvar <- attr(tab, "color.by")
  mapping <- if (!is.null(colvar)) {
    ggplot2::aes(x = .data$mean_abundance, y = .data$prevalence,
                 colour = .data[[colvar]])
  } else {
    ggplot2::aes(x = .data$mean_abundance, y = .data$prevalence)
  }
  p <- ggplot2::ggplot(tab, mapping) +
    ggplot2::geom_hline(yintercept = prevalence, linetype = "dashed",
                        colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = detection, linetype = "dashed",
                        colour = "grey70", linewidth = 0.3) +
    .mp_geom(ggplot2::geom_point,
             params = list(alpha = 0.8, size = 2), dots = list(...)) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Mean relative abundance (log scale)",
                  y = "Prevalence (fraction of samples)",
                  title = plot.title, subtitle = plot.subtitle,
                  caption = caption) +
    theme
  if (!is.null(colvar)) {
    ng <- length(unique(stats::na.omit(tab[[colvar]])))
    cols <- .mp_resolve_palette(palette, ng, type = "qualitative")
    p <- p + ggplot2::scale_colour_manual(values = cols, name = colvar)
  }
  top <- tab[order(-tab$prevalence, -tab$mean_abundance), , drop = FALSE]
  top <- utils::head(top[top$is_core, , drop = FALSE], label.top)
  if (nrow(top) > 0 && label.top > 0) {
    .mp_require("ggrepel", "Feature labelling")
    p <- p + ggrepel::geom_text_repel(
      data = top, ggplot2::aes(label = .data$feature),
      size = font.size / 3.5, max.overlaps = 30, show.legend = FALSE,
      inherit.aes = TRUE)
  }
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

#' Attach taxonomy labels (feature names) and optional colour rank
#' @noRd
.mp_attach_taxonomy <- function(tab, object, tax.level, color.by) {
  tt <- tryCatch(as.data.frame(unclass(phyloseq::tax_table(object)),
                               stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (!is.null(tt)) {
    # Colour rank is read off the original OTU ids, before relabelling.
    if (!is.null(color.by) && color.by %in% names(tt)) {
      tab[[color.by]] <- tt[match(tab$feature, rownames(tt)), color.by]
      attr(tab, "color.by") <- color.by
    }
    if (!is.null(tax.level) && tax.level %in% names(tt)) {
      lab <- tt[match(tab$feature, rownames(tt)), tax.level]
      tab$feature <- ifelse(is.na(lab), tab$feature, as.character(lab))
    }
  }
  tab
}

#' Build a detection-threshold x prevalence grid for the core heatmap
#' @noRd
.mp_core_grid <- function(mat, prevalence) {
  dets <- signif(10^seq(-4, -1, length.out = 8), 2)
  taxa <- colnames(mat)
  nsmp <- nrow(mat)
  grid <- expand.grid(feature = taxa, detection = dets,
                      stringsAsFactors = FALSE)
  grid$prevalence <- mapply(function(f, d) sum(mat[, f] > d) / nsmp,
                            grid$feature, grid$detection)
  # Keep taxa that are core somewhere along the gradient.
  core_taxa <- unique(grid$feature[grid$prevalence >= prevalence])
  grid <- grid[grid$feature %in% core_taxa, , drop = FALSE]
  ord <- names(sort(tapply(grid$prevalence, grid$feature, max)))
  grid$feature <- factor(grid$feature, levels = ord)
  grid
}

#' Render the core-microbiome heatmap
#' @noRd
.mp_core_heatmap <- function(core_df, palette, theme, plot.title,
                             plot.subtitle, plot.caption, prevalence, dots) {
  if (nrow(core_df) == 0) {
    rlang::abort(
      c("No core taxa found at the requested prevalence.",
        i = "Lower `prevalence` or `detection`."),
      class = "micropubr_no_core")
  }
  cols <- .mp_resolve_palette(palette, 256, type = "sequential")
  ggplot2::ggplot(core_df, ggplot2::aes(x = factor(.data$detection),
                                        y = .data$feature,
                                        fill = .data$prevalence)) +
    .mp_geom(ggplot2::geom_tile,
             params = list(colour = "grey95", linewidth = 0.1), dots = dots) +
    ggplot2::scale_fill_gradientn(colours = cols, name = "Prevalence",
                                  limits = c(0, 1)) +
    ggplot2::labs(x = "Detection threshold (relative abundance)", y = NULL,
                  title = plot.title, subtitle = plot.subtitle,
                  caption = plot.caption %||% sprintf(
                    "Core microbiome (prevalence \u2265 %.0f%%).",
                    100 * prevalence)) +
    theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# --- Indicator species --------------------------------------------------

#' Indicator species analysis
#'
#' Identifies taxa significantly associated with groups (or group
#' combinations) via `indicspecies::multipatt()` (IndVal.g), and plots the
#' strongest indicators coloured by their associated group.
#'
#' @param object A [phyloseq::phyloseq] object.
#' @param group.by Sample-metadata column defining the groups.
#' @param tax.level Optional rank to agglomerate to first.
#' @param top.n Maximum number of indicator taxa to display.
#' @param p.threshold Significance cut-off on the permutation p-value.
#' @param permutations Number of permutations.
#' @param palette,theme,plot.title,plot.subtitle,plot.caption,font.size See
#'   [mp_alpha()].
#' @param ... Passed to the primary geom.
#' @param interactive If `TRUE`, return an interactive `plotly` widget instead
#'   of a static `ggplot` (requires `plotly`). See [mp_alpha()].
#' @param return.data If `TRUE`, return the indicator-statistics table.
#' @param verbose Print the methods used.
#'
#' @return A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
#'   data frame when `return.data = TRUE`.
#' @export
mp_indicator <- function(object,
                         group.by = NULL,
                         tax.level = NULL,
                         top.n = 20,
                         p.threshold = 0.05,
                         permutations = 999,
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
  .mp_require(c("indicspecies", "permute"), "Indicator species analysis")
  .mp_check_phyloseq(object)
  group.by <- .mp_check_var(object, group.by, "group.by", required = TRUE)
  theme <- .mp_theme(theme, font.size)

  obj <- .mp_agglomerate(object, tax.level)
  comm <- .mp_otu_matrix(obj, samples_as_rows = TRUE)
  meta <- .mp_sample_data(obj)[rownames(comm), , drop = FALSE]
  g <- meta[[group.by]]
  keep <- !is.na(g)
  comm <- comm[keep, , drop = FALSE]
  g <- factor(g[keep])

  mp <- indicspecies::multipatt(
    comm, g, func = "IndVal.g", duleg = TRUE,
    control = permute::how(nperm = permutations))
  sign <- mp$sign
  sign$feature <- rownames(sign)
  grpcols <- grep("^s\\.", names(sign), value = TRUE)
  groups_of <- sub("^s\\.", "", grpcols)
  sign$association <- apply(sign[, grpcols, drop = FALSE], 1, function(r)
    paste(groups_of[r == 1], collapse = "+"))

  sig <- sign[!is.na(sign$p.value) & sign$p.value < p.threshold, ,
              drop = FALSE]
  sig <- sig[order(-sig$stat), , drop = FALSE]
  sig <- utils::head(sig, top.n)
  sig <- .mp_relabel_features(sig, obj, tax.level)

  if (verbose) {
    rlang::inform(sprintf("mp_indicator(): %d significant indicators (p < %g).",
                          nrow(sig), p.threshold))
  }
  if (isTRUE(return.data)) return(sig)
  if (nrow(sig) == 0) {
    rlang::abort("No significant indicator species found.",
                 class = "micropubr_no_indicator")
  }

  sig$feature <- stats::reorder(sig$feature, sig$stat)
  ng <- length(unique(sig$association))
  cols <- .mp_resolve_palette(palette, ng, type = "qualitative")
  p <- ggplot2::ggplot(sig, ggplot2::aes(x = .data$stat, y = .data$feature,
                                         fill = .data$association)) +
    .mp_geom(ggplot2::geom_col, params = list(width = 0.8), dots = list(...)) +
    ggplot2::scale_fill_manual(values = cols, name = group.by) +
    ggplot2::labs(x = "Indicator value (IndVal.g)", y = NULL,
                  title = plot.title, subtitle = plot.subtitle,
                  caption = plot.caption %||% sprintf(
                    "Indicator species (p < %g, %d permutations).",
                    p.threshold, permutations)) +
    theme
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

#' Map OTU ids in a `feature` column to a taxonomy label
#' @noRd
.mp_relabel_features <- function(df, object, tax.level) {
  if (is.null(tax.level)) return(df)
  tt <- tryCatch(as.data.frame(unclass(phyloseq::tax_table(object)),
                               stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (!is.null(tt) && tax.level %in% names(tt)) {
    lab <- tt[df$feature, tax.level]
    df$feature <- ifelse(is.na(lab), df$feature, as.character(lab))
  }
  df
}

# --- Run differential abundance -----------------------------------------

#' Run a differential-abundance analysis from a phyloseq object
#'
#' Executes a two-group differential-abundance test and returns a tidy data
#' frame (`feature`, `log2FC`, `pvalue`, `padj`, `mean`) ready to pass to
#' [mp_diff()]. The default `"wilcox"` method (Wilcoxon test on centred
#' log-ratios) needs no extra packages; `"deseq2"`, `"aldex2"` and `"ancombc"`
#' wrap the respective Bioconductor packages.
#'
#' @param object A [phyloseq::phyloseq] object with counts.
#' @param group.by Sample-metadata column defining the groups.
#' @param method One of `"wilcox"`, `"deseq2"`, `"aldex2"`, `"ancombc"`.
#' @param groups Optional length-2 vector naming the two groups to contrast
#'   (reference first). Required when `group.by` has more than two levels.
#' @param tax.level Optional rank to agglomerate to first.
#' @param verbose Print the methods used.
#'
#' @return A data frame compatible with [mp_diff()].
#' @export
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' res <- mp_run_da(GlobalPatterns, group.by = "SampleType",
#'                  groups = c("Feces", "Soil"), method = "wilcox")
#' mp_diff(res, plot.type = "volcano")
#' }
mp_run_da <- function(object,
                      group.by = NULL,
                      method = c("wilcox", "deseq2", "aldex2", "ancombc"),
                      groups = NULL,
                      tax.level = NULL,
                      verbose = getOption("micropubr.verbose", FALSE)) {
  .mp_check_phyloseq(object)
  method <- match.arg(method)
  group.by <- .mp_check_var(object, group.by, "group.by", required = TRUE)
  obj <- .mp_agglomerate(object, tax.level)
  meta <- .mp_sample_data(obj)
  lv <- levels(factor(meta[[group.by]]))
  if (!is.null(groups)) {
    if (length(groups) != 2 || !all(groups %in% lv)) {
      rlang::abort("`groups` must name two levels of `group.by`.",
                   class = "micropubr_bad_groups")
    }
    lv <- groups
  } else if (length(lv) != 2) {
    rlang::abort(
      c(sprintf("`group.by` has %d levels; specify `groups` to pick two.",
                length(lv)),
        i = sprintf("Levels: %s.", paste(lv, collapse = ", "))),
      class = "micropubr_need_groups")
  }
  if (verbose) {
    rlang::inform(sprintf("mp_run_da(): %s, %s vs %s.",
                          method, lv[1], lv[2]))
  }
  out <- switch(
    method,
    wilcox  = .mp_da_run_wilcox(obj, group.by, lv, tax.level),
    deseq2  = .mp_da_run_deseq2(obj, group.by, lv),
    aldex2  = .mp_da_run_aldex2(obj, group.by, lv),
    ancombc = .mp_da_run_ancombc(obj, group.by, lv))
  out[stats::complete.cases(out$log2FC), , drop = FALSE]
}

#' Wilcoxon test on CLR values (dependency-free default)
#' @noRd
.mp_da_run_wilcox <- function(object, group.by, lv, tax.level) {
  clr <- .mp_clr(object)                                  # samples x taxa
  rel <- .mp_otu_matrix(.mp_to_relative(object), samples_as_rows = TRUE)
  meta <- .mp_sample_data(object)[rownames(clr), , drop = FALSE]
  g <- meta[[group.by]]
  keep <- g %in% lv & !is.na(g)
  clr <- clr[keep, , drop = FALSE]
  rel <- rel[keep, , drop = FALSE]
  g <- factor(g[keep], levels = lv)
  i1 <- g == lv[1]
  i2 <- g == lv[2]
  # log2 fold-change as the difference in mean CLR (natural log) over log(2).
  l2fc <- (colMeans(clr[i2, , drop = FALSE]) -
             colMeans(clr[i1, , drop = FALSE])) / log(2)
  pv <- apply(clr, 2, function(x) tryCatch(
    stats::wilcox.test(x[i1], x[i2])$p.value,
    error = function(e) NA_real_))
  out <- data.frame(feature = colnames(clr),
                    log2FC = as.numeric(l2fc),
                    pvalue = as.numeric(pv),
                    padj = stats::p.adjust(pv, method = "BH"),
                    mean = as.numeric(colMeans(rel)),
                    stringsAsFactors = FALSE)
  .mp_relabel_features(out, object, tax.level)
}

#' DESeq2 wrapper
#' @noRd
.mp_da_run_deseq2 <- function(object, group.by, lv) {
  .mp_require("DESeq2", "DESeq2 differential abundance", bioc = TRUE)
  sub <- phyloseq::prune_samples(
    phyloseq::sample_data(object)[[group.by]] %in% lv, object)
  phyloseq::sample_data(sub)[[group.by]] <-
    factor(phyloseq::sample_data(sub)[[group.by]], levels = lv)
  dds <- phyloseq::phyloseq_to_deseq2(
    sub, stats::as.formula(sprintf("~ `%s`", group.by)))
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c(group.by, lv[2], lv[1]))
  .mp_parse_da(res)
}

#' ALDEx2 wrapper
#' @noRd
.mp_da_run_aldex2 <- function(object, group.by, lv) {
  .mp_require("ALDEx2", "ALDEx2 differential abundance", bioc = TRUE)
  meta <- .mp_sample_data(object)
  g <- meta[[group.by]]
  keep <- g %in% lv & !is.na(g)
  counts <- .mp_otu_matrix(object, samples_as_rows = FALSE)[, keep,
                                                            drop = FALSE]
  conds <- as.character(g[keep])
  ax <- ALDEx2::aldex(round(counts), conds, test = "t", effect = TRUE,
                      verbose = FALSE)
  ax$feature <- rownames(ax)
  .mp_parse_da(ax)
}

#' ANCOM-BC wrapper
#' @noRd
.mp_da_run_ancombc <- function(object, group.by, lv) {
  .mp_require("ANCOMBC", "ANCOM-BC differential abundance", bioc = TRUE)
  sub <- phyloseq::prune_samples(
    phyloseq::sample_data(object)[[group.by]] %in% lv, object)
  res <- ANCOMBC::ancombc(data = sub, formula = group.by,
                          group = group.by)
  .mp_parse_da(res)
}
