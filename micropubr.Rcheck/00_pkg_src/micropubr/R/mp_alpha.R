# Alpha diversity ---------------------------------------------------------

#' Plot alpha (within-sample) diversity
#'
#' Computes one or more alpha-diversity indices from a [phyloseq::phyloseq]
#' object and draws publication-ready box, violin or dot plots, optionally
#' annotated with pairwise statistical tests. Rarefaction curves are also
#' supported.
#'
#' @param object A [phyloseq::phyloseq] object containing **count** data.
#' @param measure Character vector of indices. Any of `"Observed"`,
#'   `"Chao1"`, `"ACE"`, `"Shannon"`, `"Simpson"`, `"InvSimpson"`,
#'   `"Fisher"`, `"Pielou"` (evenness) and `"Faith"` (phylogenetic diversity;
#'   requires a tree and the `picante` package).
#' @param group.by Sample-metadata column to compare across (x-axis / fill).
#'   If `NULL`, all samples are pooled into a single group and no test is run.
#' @param plot.type One of `"box"`, `"violin"`, `"dot"` or `"rarefaction"`.
#' @param test Statistical test for group comparisons: `"wilcoxon"`,
#'   `"kruskal"`, `"t.test"` or `"anova"`. Pairwise post-hoc tests are Dunn
#'   (for `"kruskal"`) and Tukey HSD (for `"anova"`).
#' @param add.stats Add significance brackets / p-values when `group.by` has
#'   at least two levels.
#' @param comparisons Optional list of length-2 character vectors specifying
#'   which group pairs to test. If `NULL`, all pairwise comparisons are used.
#' @param p.adjust.method Multiple-testing correction (default `"BH"`).
#' @param facet.by Optional metadata column to facet each panel by.
#' @param palette Palette name, explicit colour vector, or `NULL` for the
#'   default colourblind-safe palette. See [mp_palette()].
#' @param point.size Size of jittered sample points.
#' @param jitter Overlay individual sample points.
#' @param theme A `ggplot2` theme; defaults to the active theme (see
#'   [mp_set_defaults()], otherwise [mp_theme_publication()]).
#' @param ... Additional arguments passed to the primary geom (boxplot,
#'   violin or, for `plot.type = "dot"`, the point summary).
#' @param plot.title,plot.subtitle,plot.caption Optional labels. When `NULL`,
#'   sensible defaults are generated (e.g. the omnibus test result as
#'   subtitle).
#' @param font.size Base font size in points.
#' @param combine If multiple `measure`s are requested, combine the panels
#'   into a single [patchwork][patchwork::patchwork-package] object (`TRUE`)
#'   or return a named list of `ggplot` objects (`FALSE`).
#' @param return.data If `TRUE`, return the tidy data frame of per-sample
#'   diversity estimates and metadata instead of a plot.
#' @param verbose Print the methods used.
#'
#' @return A `ggplot` object (single measure), a `patchwork` (multiple
#'   measures, `combine = TRUE`) or a named list of `ggplot`s; or a data frame
#'   when `return.data = TRUE`.
#' @export
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' mp_alpha(GlobalPatterns, measure = c("Shannon", "Observed"),
#'          group.by = "SampleType")
#' }
mp_alpha <- function(object,
                     measure = c("Shannon", "Observed"),
                     group.by = NULL,
                     plot.type = c("box", "violin", "dot", "rarefaction"),
                     test = c("wilcoxon", "kruskal", "t.test", "anova"),
                     add.stats = TRUE,
                     comparisons = NULL,
                     p.adjust.method = "BH",
                     facet.by = NULL,
                     palette = NULL,
                     point.size = 1.5,
                     jitter = TRUE,
                     theme = NULL,
                     ...,
                     plot.title = NULL,
                     plot.subtitle = NULL,
                     plot.caption = NULL,
                     font.size = getOption("micropubr.font.size", 11),
                     combine = TRUE,
                     return.data = FALSE,
                     verbose = getOption("micropubr.verbose", FALSE)) {
  .mp_check_phyloseq(object)
  plot.type <- match.arg(plot.type)
  test <- match.arg(test)
  group.by <- .mp_check_var(object, group.by, "group.by")
  facet.by <- .mp_check_var(object, facet.by, "facet.by")
  theme <- .mp_theme(theme, font.size)

  if (plot.type == "rarefaction") {
    return(.mp_rarefaction(object, group.by = group.by, palette = palette,
                           theme = theme, plot.title = plot.title,
                           plot.subtitle = plot.subtitle,
                           plot.caption = plot.caption, verbose = verbose))
  }

  alpha <- .mp_estimate_alpha(object, measure, verbose = verbose)
  meta <- .mp_sample_data(object)
  df <- cbind(meta[rownames(alpha), , drop = FALSE], alpha)

  if (isTRUE(return.data)) return(df)

  if (is.null(group.by)) {
    df[[".all"]] <- "All samples"
    group.by <- ".all"
    add.stats <- FALSE
  }
  df[[group.by]] <- as.factor(df[[group.by]])

  if (verbose) {
    rlang::inform(sprintf(
      "mp_alpha(): %s plot of %s; grouped by '%s'%s.",
      plot.type, paste(measure, collapse = ", "), group.by,
      if (add.stats) sprintf(", %s test", test) else ""))
  }

  plots <- stats::setNames(vector("list", length(measure)), measure)
  for (m in measure) {
    plots[[m]] <- .mp_alpha_panel(
      df, measure = m, group = group.by, facet.by = facet.by,
      plot.type = plot.type, test = test, add.stats = add.stats,
      comparisons = comparisons, p.adjust.method = p.adjust.method,
      palette = palette, point.size = point.size, jitter = jitter,
      theme = theme, plot.subtitle = plot.subtitle,
      plot.caption = plot.caption, dots = list(...)
    )
  }

  if (length(plots) == 1) {
    p <- plots[[1]]
    if (!is.null(plot.title)) p <- p + ggplot2::labs(title = plot.title)
    return(p)
  }
  if (combine) {
    .mp_require("patchwork", "Combining alpha-diversity panels")
    combined <- patchwork::wrap_plots(plots) +
      patchwork::plot_layout(guides = "collect")
    if (!is.null(plot.title)) {
      combined <- combined +
        patchwork::plot_annotation(title = plot.title)
    }
    return(combined)
  }
  plots
}

#' Estimate alpha diversity, including derived Pielou and Faith indices
#' @noRd
.mp_estimate_alpha <- function(object, measure, verbose = FALSE) {
  core <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson",
            "InvSimpson", "Fisher")
  want_core <- intersect(measure, core)
  # "Pielou" needs Shannon + Observed; ensure both are computed.
  need <- want_core
  if ("Pielou" %in% measure) need <- union(need, c("Shannon", "Observed"))

  est <- if (length(need)) {
    suppressWarnings(phyloseq::estimate_richness(object, measures = need))
  } else {
    data.frame(row.names = phyloseq::sample_names(object))
  }
  # phyloseq mangles sample names with make.names(); restore original order.
  rownames(est) <- phyloseq::sample_names(object)

  out <- est[, intersect(want_core, colnames(est)), drop = FALSE]

  if ("Pielou" %in% measure) {
    out[["Pielou"]] <- est[["Shannon"]] / log(est[["Observed"]])
  }
  if ("Faith" %in% measure) {
    out[["Faith"]] <- .mp_faith_pd(object, verbose = verbose)
  }
  out[, measure[measure %in% colnames(out)], drop = FALSE]
}

#' Faith's phylogenetic diversity via picante (optional)
#' @noRd
.mp_faith_pd <- function(object, verbose = FALSE) {
  .mp_check_phyloseq(object, require_tree = TRUE)
  .mp_require("picante", "Faith's phylogenetic diversity")
  comm <- .mp_otu_matrix(object, samples_as_rows = TRUE)
  tree <- phyloseq::phy_tree(object)
  pd <- picante::pd(comm, tree, include.root = FALSE)
  stats::setNames(pd$PD, rownames(comm))[phyloseq::sample_names(object)]
}

#' Build a single alpha-diversity panel
#' @noRd
.mp_alpha_panel <- function(df, measure, group, facet.by, plot.type,
                            test, add.stats, comparisons, p.adjust.method,
                            palette, point.size, jitter, theme,
                            plot.subtitle, plot.caption, dots = list()) {
  ng <- nlevels(df[[group]])
  cols <- .mp_resolve_palette(palette, ng, type = "qualitative")

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[group]], y = .data[[measure]],
                 fill = .data[[group]])
  )

  if (plot.type == "violin") {
    p <- p +
      .mp_geom(ggplot2::geom_violin,
               params = list(trim = FALSE, alpha = 0.55, colour = NA),
               dots = dots) +
      ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA,
                            fill = "white", alpha = 0.8)
  } else if (plot.type == "dot") {
    p <- p + .mp_geom(
      ggplot2::stat_summary,
      params = list(fun = mean, geom = "point", size = 3,
                    mapping = ggplot2::aes(colour = .data[[group]])),
      dots = dots
    ) + ggplot2::stat_summary(
      fun.data = ggplot2::mean_sdl, fun.args = list(mult = 1),
      geom = "errorbar", width = 0.2,
      ggplot2::aes(colour = .data[[group]])
    )
  } else {
    p <- p + .mp_geom(ggplot2::geom_boxplot,
                      params = list(outlier.shape = NA, alpha = 0.7),
                      dots = dots)
  }

  if (jitter && plot.type != "dot") {
    p <- p + ggplot2::geom_jitter(width = 0.18, height = 0,
                                  size = point.size, alpha = 0.6,
                                  ggplot2::aes(colour = .data[[group]]),
                                  show.legend = FALSE)
  }

  p <- p +
    ggplot2::scale_fill_manual(values = cols, name = group) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::labs(x = NULL, y = measure) +
    theme

  if (!is.null(facet.by)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet.by)),
                                 scales = "free_x")
  }

  # Hide redundant axis text/legend when there is a single pooled group.
  if (identical(group, ".all")) {
    p <- p + ggplot2::theme(legend.position = "none",
                            axis.text.x = ggplot2::element_blank(),
                            axis.ticks.x = ggplot2::element_blank())
  }

  subtitle <- plot.subtitle
  if (add.stats && ng >= 2 && is.null(facet.by)) {
    stat_df <- tryCatch(
      .mp_pairwise_stats(df, measure, group, test = test,
                         p.adjust.method = p.adjust.method,
                         comparisons = comparisons),
      error = function(e) NULL
    )
    if (!is.null(stat_df) && nrow(stat_df) > 0) {
      p <- p + ggpubr::stat_pvalue_manual(
        stat_df, label = "p.adj.signif", tip.length = 0.01,
        hide.ns = TRUE, inherit.aes = FALSE
      )
    }
    if (is.null(subtitle)) {
      subtitle <- .mp_global_test(df, measure, group, test)
    }
  }
  p + ggplot2::labs(subtitle = subtitle, caption = plot.caption)
}

#' Rarefaction-curve panel built from vegan::rarecurve output
#' @noRd
.mp_rarefaction <- function(object, group.by, palette, theme,
                            plot.title, plot.subtitle, plot.caption,
                            step = NULL, verbose = FALSE) {
  .mp_require("vegan", "Rarefaction curves")
  mat <- round(.mp_otu_matrix(object, samples_as_rows = TRUE))
  storage.mode(mat) <- "integer"
  if (is.null(step)) step <- max(1, floor(min(rowSums(mat)) / 20))

  curves <- vegan::rarecurve(mat, step = step, label = FALSE, tidy = TRUE)
  names(curves) <- c(".sample", "Sample_size", "Species")

  meta <- .mp_sample_data(object)
  if (!is.null(group.by)) {
    curves[[group.by]] <- meta[as.character(curves$.sample), group.by]
  }

  aes_line <- if (!is.null(group.by)) {
    ggplot2::aes(x = .data$Sample_size, y = .data$Species,
                 group = .data$.sample, colour = .data[[group.by]])
  } else {
    ggplot2::aes(x = .data$Sample_size, y = .data$Species,
                 group = .data$.sample)
  }

  p <- ggplot2::ggplot(curves, aes_line) +
    ggplot2::geom_line(alpha = 0.8, linewidth = 0.5) +
    ggplot2::labs(x = "Sequencing depth", y = "Observed richness",
                  title = plot.title, subtitle = plot.subtitle,
                  caption = plot.caption %||%
                    "Rarefaction curves (observed richness)") +
    theme

  if (!is.null(group.by)) {
    ng <- length(unique(stats::na.omit(curves[[group.by]])))
    cols <- .mp_resolve_palette(palette, ng, type = "qualitative")
    p <- p + ggplot2::scale_colour_manual(values = cols, name = group.by)
  }
  p
}
