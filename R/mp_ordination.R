# Beta diversity and ordination ------------------------------------------

#' Ordination plot with ellipses and PERMANOVA
#'
#' Computes an ordination from a [phyloseq::phyloseq] object and renders a
#' publication-ready scatter plot. Distance-based methods (PCoA, NMDS, dbRDA)
#' accept any [phyloseq::distance()] metric; PCA operates on centred
#' log-ratio (Aitchison) coordinates; RDA and CCA support constrained
#' ordination against metadata. Group ellipses and a PERMANOVA result are
#' added automatically.
#'
#' @param object A [phyloseq::phyloseq] object.
#' @param method One of `"PCoA"`, `"NMDS"`, `"PCA"`, `"RDA"`, `"CCA"`,
#'   `"dbRDA"`.
#' @param distance Distance metric for distance-based methods. Any value
#'   accepted by [phyloseq::distance()] (e.g. `"bray"`, `"jaccard"`,
#'   `"unifrac"`, `"wunifrac"`, `"canberra"`) plus `"aitchison"`.
#' @param group.by Metadata column mapped to point colour (and used for
#'   ellipses and PERMANOVA).
#' @param shape.by Optional metadata column mapped to point shape.
#' @param ellipse Draw group ellipses.
#' @param ellipse.type Ellipse type passed to [ggplot2::stat_ellipse()].
#' @param ellipse.level Confidence level for ellipses.
#' @param permanova Run [vegan::adonis2()] and append the result to the
#'   caption (distance-based methods only).
#' @param permutations Number of PERMANOVA permutations.
#' @param constraints For `"RDA"`, `"CCA"` and `"dbRDA"`, a character vector
#'   of metadata columns to constrain the ordination on.
#' @param add.vectors Overlay constraint/biplot arrows (constrained methods).
#' @param envfit Optional character vector of (numeric) metadata columns to
#'   fit onto the ordination as vectors via [vegan::envfit()]. Significant
#'   vectors are overlaid as arrows and summarised in the caption.
#' @param axes Length-2 integer vector of axes to plot.
#' @param label.samples Label individual sample points (via `ggrepel`).
#' @param point.size Point size.
#' @param point.alpha Point opacity.
#' @param palette Palette name, colour vector, or `NULL`. See [mp_palette()].
#' @param theme A `ggplot2` theme; defaults to the active theme (see
#'   [mp_set_defaults()], otherwise [mp_theme_publication()]).
#' @param ... Additional arguments passed to the point geom
#'   ([ggplot2::geom_point()]).
#' @param plot.title,plot.subtitle,plot.caption Optional labels; the caption
#'   is auto-generated when `NULL`.
#' @param font.size Base font size.
#' @param return.data If `TRUE`, return the underlying data frame of sample
#'   scores and metadata (with the ordination object and any PERMANOVA/envfit
#'   results attached as attributes) instead of a plot.
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
#' mp_ordination(GlobalPatterns, method = "PCoA", distance = "bray",
#'               group.by = "SampleType")
#' }
mp_ordination <- function(object,
                          method = c("PCoA", "NMDS", "PCA", "RDA",
                                     "CCA", "dbRDA"),
                          distance = "bray",
                          group.by = NULL,
                          shape.by = NULL,
                          ellipse = TRUE,
                          ellipse.type = c("norm", "t", "euclid"),
                          ellipse.level = 0.95,
                          permanova = TRUE,
                          permutations = 999,
                          constraints = NULL,
                          add.vectors = FALSE,
                          envfit = NULL,
                          axes = c(1, 2),
                          label.samples = FALSE,
                          point.size = 2.5,
                          point.alpha = 0.85,
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
  method <- match.arg(method)
  ellipse.type <- match.arg(ellipse.type)
  group.by <- .mp_check_var(object, group.by, "group.by")
  shape.by <- .mp_check_var(object, shape.by, "shape.by")
  if (!is.null(constraints)) {
    constraints <- .mp_check_var(object, constraints, "constraints")
  }
  if (!is.null(envfit)) envfit <- .mp_check_var(object, envfit, "envfit")
  theme <- .mp_theme(theme, font.size)
  stopifnot(length(axes) == 2)

  ord <- .mp_ordinate(object, method = method, distance = distance,
                      constraints = constraints, axes = axes,
                      verbose = verbose)

  meta <- .mp_sample_data(object)
  df <- cbind(ord$scores, meta[rownames(ord$scores), , drop = FALSE])
  xvar <- names(ord$scores)[1]
  yvar <- names(ord$scores)[2]

  # Optional PERMANOVA, computed up front so it can be attached to the data
  # frame (return.data) as well as the caption.
  perm <- NULL
  if (permanova && !is.null(group.by) && !is.null(ord$dist)) {
    perm <- .mp_permanova(ord$dist, meta, group.by, permutations)
  }
  # Optional envfit of metadata vectors onto the two plotted axes.
  ef <- NULL
  if (!is.null(envfit)) {
    ef <- .mp_envfit(ord$scores, meta, envfit, permutations)
  }

  if (isTRUE(return.data)) {
    attr(df, "ordination") <- ord[setdiff(names(ord), "scores")]
    if (!is.null(perm)) attr(df, "permanova") <- perm$model
    if (!is.null(ef)) attr(df, "envfit") <- ef$table
    return(df)
  }

  mapping <- ggplot2::aes(x = .data[[xvar]], y = .data[[yvar]])
  if (!is.null(group.by)) {
    df[[group.by]] <- as.factor(df[[group.by]])
    mapping$colour <- ggplot2::aes(colour = .data[[group.by]])$colour
  }
  if (!is.null(shape.by)) {
    df[[shape.by]] <- as.factor(df[[shape.by]])
    mapping$shape <- ggplot2::aes(shape = .data[[shape.by]])$shape
  }

  p <- ggplot2::ggplot(df, mapping) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey85", linewidth = 0.3) +
    .mp_geom(ggplot2::geom_point,
             params = list(size = point.size, alpha = point.alpha),
             dots = list(...)) +
    ggplot2::labs(x = ord$xlab, y = ord$ylab) +
    theme

  if (!is.null(group.by)) {
    ng <- nlevels(df[[group.by]])
    cols <- .mp_resolve_palette(palette, ng, type = "qualitative")
    p <- p + ggplot2::scale_colour_manual(values = cols, name = group.by)
    if (ellipse && ng >= 2) {
      # stat_ellipse needs >= 3 points per group; restrict to those groups
      # to avoid noisy "Too few points" warnings.
      big <- names(which(table(df[[group.by]]) >= 3))
      ell_df <- df[df[[group.by]] %in% big, , drop = FALSE]
      if (nrow(ell_df) > 0) {
        p <- p + ggplot2::stat_ellipse(
          data = ell_df, type = ellipse.type,
          level = ellipse.level, linewidth = 0.5
        )
      }
    }
  }

  if (add.vectors && !is.null(ord$biplot)) {
    p <- .mp_add_biplot(p, ord$biplot, xvar, yvar)
  }

  if (!is.null(ef) && !is.null(ef$arrows) && nrow(ef$arrows) > 0) {
    p <- .mp_add_biplot(p, ef$arrows, xvar, yvar)
  }

  if (label.samples) {
    .mp_require("ggrepel", "Sample labelling")
    df[[".lab"]] <- rownames(df)
    p <- p + ggrepel::geom_text_repel(
      data = df, ggplot2::aes(label = .data[[".lab"]]),
      size = font.size / 4, max.overlaps = 20, inherit.aes = TRUE,
      show.legend = FALSE
    )
  }

  # Auto caption: methods line + optional PERMANOVA + optional envfit.
  caption <- plot.caption %||% ord$caption
  if (!is.null(perm)) {
    caption <- paste(caption, perm$caption, sep = "\n")
    if (verbose) rlang::inform(perm$caption)
  }
  if (!is.null(ef) && !is.null(ef$caption)) {
    caption <- paste(caption, ef$caption, sep = "\n")
    if (verbose) rlang::inform(ef$caption)
  }

  p <- p + ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                         caption = caption)
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

#' Dispatch ordination computation by method
#'
#' @return A list with `scores` (data.frame, 2 cols named by axis labels),
#'   `xlab`, `ylab`, `caption`, optional `dist` (for PERMANOVA) and optional
#'   `biplot` (constraint arrows).
#' @noRd
.mp_ordinate <- function(object, method, distance, constraints, axes,
                         verbose = FALSE) {
  if (verbose) {
    rlang::inform(sprintf("mp_ordination(): %s%s.", method,
                          if (method %in% c("PCoA", "NMDS", "dbRDA"))
                            sprintf(" on %s distances", distance) else ""))
  }
  switch(
    method,
    PCoA  = .mp_ord_pcoa(object, distance, axes),
    NMDS  = .mp_ord_nmds(object, distance),
    PCA   = .mp_ord_pca(object, axes),
    RDA   = .mp_ord_constrained(object, "rda", constraints, axes),
    CCA   = .mp_ord_constrained(object, "cca", constraints, axes),
    dbRDA = .mp_ord_constrained(object, "dbrda", constraints, axes,
                                distance = distance)
  )
}

#' Compute a sample distance, including Aitchison
#' @noRd
.mp_compute_distance <- function(object, distance) {
  if (identical(distance, "aitchison")) {
    return(stats::dist(.mp_clr(object), method = "euclidean"))
  }
  phyloseq::distance(object, method = distance)
}

#' @noRd
.mp_ord_pcoa <- function(object, distance, axes) {
  d <- .mp_compute_distance(object, distance)
  k <- max(axes)
  cmd <- stats::cmdscale(d, k = k, eig = TRUE)
  pts <- cmd$points[, axes, drop = FALSE]
  eig <- cmd$eig
  pct <- 100 * eig / sum(eig[eig > 0])
  lab <- sprintf("PCoA%d (%.1f%%)", axes, pct[axes])
  scores <- as.data.frame(pts)
  names(scores) <- lab
  rownames(scores) <- labels(d)
  list(scores = scores, xlab = lab[1], ylab = lab[2], dist = d,
       caption = sprintf("PCoA of %s distances (n = %d).",
                         distance, nrow(scores)))
}

#' @noRd
.mp_ord_nmds <- function(object, distance) {
  .mp_require("vegan", "NMDS ordination")
  d <- .mp_compute_distance(object, distance)
  nmds <- vegan::metaMDS(d, k = 2, trace = 0, autotransform = FALSE)
  pts <- vegan::scores(nmds, display = "sites")
  scores <- as.data.frame(pts)
  names(scores) <- c("NMDS1", "NMDS2")
  rownames(scores) <- labels(d)
  list(scores = scores, xlab = "NMDS1", ylab = "NMDS2", dist = d,
       caption = sprintf("NMDS of %s distances (stress = %.3f, n = %d).",
                         distance, nmds$stress, nrow(scores)))
}

#' @noRd
.mp_ord_pca <- function(object, axes) {
  clr <- .mp_clr(object)
  pr <- stats::prcomp(clr, center = TRUE, scale. = FALSE)
  pct <- 100 * pr$sdev^2 / sum(pr$sdev^2)
  pts <- pr$x[, axes, drop = FALSE]
  lab <- sprintf("PC%d (%.1f%%)", axes, pct[axes])
  scores <- as.data.frame(pts)
  names(scores) <- lab
  rownames(scores) <- rownames(clr)
  list(scores = scores, xlab = lab[1], ylab = lab[2], dist = NULL,
       caption = sprintf("PCA of CLR-transformed (Aitchison) data (n = %d).",
                         nrow(scores)))
}

#' Constrained ordination: rda / cca / dbrda
#' @noRd
.mp_ord_constrained <- function(object, engine, constraints, axes,
                                distance = "bray") {
  .mp_require("vegan", "Constrained ordination")
  if (is.null(constraints)) {
    rlang::abort(
      c(sprintf("`%s` requires `constraints`.", toupper(engine)),
        i = "Supply one or more metadata columns, e.g. constraints = \"Group\"."),
      class = "micropubr_missing_constraints"
    )
  }
  meta <- .mp_sample_data(object)
  rhs <- paste(sprintf("`%s`", constraints), collapse = " + ")

  if (engine == "dbrda") {
    d <- .mp_compute_distance(object, distance)
    meta <- meta[labels(d), , drop = FALSE]
    f <- stats::as.formula(paste("d ~", rhs))
    mod <- vegan::dbrda(f, data = meta)
    prefix <- "dbRDA"
    dist_out <- d
  } else {
    comm <- .mp_otu_matrix(object, samples_as_rows = TRUE)
    meta <- meta[rownames(comm), , drop = FALSE]
    f <- stats::as.formula(paste("comm ~", rhs))
    mod <- if (engine == "rda") vegan::rda(f, data = meta) else
      vegan::cca(f, data = meta)
    prefix <- toupper(engine)
    dist_out <- NULL
  }

  sit <- vegan::scores(mod, display = "sites", choices = axes)
  eig <- vegan::eigenvals(mod)
  pct <- 100 * eig / sum(eig)
  lab <- sprintf("%s%d (%.1f%%)", prefix, axes, pct[axes])
  scores <- as.data.frame(sit)
  names(scores) <- lab

  bip <- tryCatch(
    as.data.frame(vegan::scores(mod, display = "bp", choices = axes)),
    error = function(e) NULL
  )
  if (!is.null(bip)) {
    names(bip) <- lab
    bip$.label <- rownames(bip)
  }

  list(scores = scores, xlab = lab[1], ylab = lab[2], dist = dist_out,
       biplot = bip,
       caption = sprintf("%s constrained by %s (n = %d).",
                         prefix, paste(constraints, collapse = " + "),
                         nrow(scores)))
}

#' Overlay constraint/biplot arrows on an ordination plot
#' @noRd
.mp_add_biplot <- function(p, bip, xvar, yvar) {
  .mp_require("ggrepel", "Biplot arrow labels")
  mult <- 0.8 * min(
    diff(range(p$data[[xvar]])) / diff(range(c(0, bip[[xvar]]))),
    diff(range(p$data[[yvar]])) / diff(range(c(0, bip[[yvar]])))
  )
  bip$.xend <- bip[[xvar]] * mult
  bip$.yend <- bip[[yvar]] * mult
  p +
    ggplot2::geom_segment(
      data = bip,
      ggplot2::aes(x = 0, y = 0, xend = .data$.xend, yend = .data$.yend),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc")),
      colour = "grey25", linewidth = 0.4, inherit.aes = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = bip,
      ggplot2::aes(x = .data$.xend, y = .data$.yend, label = .data$.label),
      colour = "grey15", size = 3, inherit.aes = FALSE
    )
}

#' Fit metadata vectors onto an ordination via vegan::envfit
#'
#' Only numeric metadata columns are fitted (as vectors). Returns a list with
#' `arrows` (a biplot-ready data.frame of significant vectors), `table` (the
#' full r2/p table) and a `caption` string, or `NULL` if nothing is fitted.
#' @noRd
.mp_envfit <- function(scores, meta, vars, permutations = 999,
                       p.cutoff = 0.05) {
  .mp_require("vegan", "Environmental vector fitting (envfit)")
  meta <- meta[rownames(scores), , drop = FALSE]
  num <- vars[vapply(vars, function(v) is.numeric(meta[[v]]), logical(1))]
  if (length(num) == 0) {
    rlang::warn("`envfit` skipped: no numeric metadata columns supplied.",
                class = "micropubr_envfit_nonnumeric")
    return(NULL)
  }
  env <- meta[, num, drop = FALSE]
  ef <- tryCatch(
    vegan::envfit(scores, env, permutations = permutations, na.rm = TRUE),
    error = function(e) NULL
  )
  if (is.null(ef) || is.null(ef$vectors)) return(NULL)

  vec <- as.data.frame(vegan::scores(ef, display = "vectors"))
  names(vec) <- names(scores)
  vec$r2 <- ef$vectors$r
  vec$p <- ef$vectors$pvals
  vec$.label <- rownames(vec)

  sig <- vec[!is.na(vec$p) & vec$p < p.cutoff, , drop = FALSE]
  # Scale arrow length by strength of association (sqrt r2), as is conventional.
  if (nrow(sig) > 0) {
    sig[[names(scores)[1]]] <- sig[[names(scores)[1]]] * sqrt(sig$r2)
    sig[[names(scores)[2]]] <- sig[[names(scores)[2]]] * sqrt(sig$r2)
  }
  caption <- if (nrow(sig) > 0) {
    sprintf("envfit (p < %.2g): %s", p.cutoff,
            paste(sprintf("%s R\u00b2=%.2f", sig$.label, sig$r2),
                  collapse = ", "))
  } else {
    sprintf("envfit: no vectors significant at p < %.2g.", p.cutoff)
  }
  list(arrows = sig[, c(names(scores), ".label"), drop = FALSE],
       table = vec[, c(".label", "r2", "p")], caption = caption)
}
