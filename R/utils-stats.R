# Internal statistical helpers -------------------------------------------
# Thin wrappers over rstatix that return tidy results ready for ggpubr
# bracket placement, plus PERMANOVA formatting for ordination captions.

#' Pairwise group comparisons returning a stat_pvalue_manual-ready table
#'
#' @param df Long data.frame.
#' @param value Numeric response column name.
#' @param group Grouping column name (coerced to factor).
#' @param test One of "wilcoxon", "t.test", "kruskal", "anova". For the
#'   omnibus tests ("kruskal"/"anova") the pairwise post-hoc is Dunn / Tukey.
#' @param p.adjust.method Passed to rstatix (default "BH").
#' @param comparisons Optional list of length-2 character vectors.
#' @return A data.frame with `group1`, `group2`, `p`, `p.adj`,
#'   `p.adj.signif` and `y.position`, or `NULL` when fewer than two groups
#'   are present.
#' @noRd
.mp_pairwise_stats <- function(df, value, group,
                               test = "wilcoxon",
                               p.adjust.method = "BH",
                               comparisons = NULL,
                               step.increase = 0.12) {
  df <- df[!is.na(df[[value]]) & !is.na(df[[group]]), , drop = FALSE]
  df[[group]] <- as.factor(df[[group]])
  if (nlevels(df[[group]]) < 2) {
    return(NULL)
  }
  # A constant response makes every test degenerate (and crashes wilcox.test).
  if (stats::sd(df[[value]], na.rm = TRUE) == 0) {
    return(NULL)
  }
  f <- stats::as.formula(sprintf("`%s` ~ `%s`", value, group))

  res <- switch(
    test,
    "wilcoxon" = rstatix::wilcox_test(df, f, comparisons = comparisons,
                                      p.adjust.method = p.adjust.method),
    "t.test"   = rstatix::t_test(df, f, comparisons = comparisons,
                                 p.adjust.method = p.adjust.method),
    "kruskal"  = rstatix::dunn_test(df, f, p.adjust.method = p.adjust.method),
    "anova"    = rstatix::tukey_hsd(df, f),
    rlang::abort(sprintf("Unknown `test` value: \"%s\".", test),
                 class = "micropubr_bad_test")
  )
  if (!is.null(comparisons) && test %in% c("kruskal", "anova")) {
    res <- .mp_filter_comparisons(res, comparisons)
  }
  # rstatix rounds these p-values to three significant digits; recompute them
  # exactly so the multiple-testing adjustment (and the stars) start from the
  # same numbers stats:: would give. Dunn's test is not rounded.
  res <- .mp_exact_pairwise_p(res, df, value, group, test, p.adjust.method)
  if (!"p.adj.signif" %in% names(res)) {
    res <- rstatix::add_significance(res, p.col = "p.adj",
                                     output.col = "p.adj.signif")
  }
  res <- rstatix::add_xy_position(res, x = group, fun = "max",
                                  step.increase = step.increase)
  res
}

#' Replace rstatix's rounded pairwise p-values with exact ones
#'
#' Wilcoxon and t-test p-values are recomputed with stats:: and re-adjusted;
#' Tukey HSD adjusted p-values are taken from stats::TukeyHSD(). The
#' significance stars are then rebuilt from the exact values.
#' @noRd
.mp_exact_pairwise_p <- function(res, df, value, group, test,
                                 p.adjust.method) {
  if (!nrow(res) || test == "kruskal") {
    if (!"p.adj" %in% names(res)) res$p.adj <- res$p
    return(res)
  }
  if (test == "anova") {
    tk <- stats::TukeyHSD(stats::aov(
      stats::as.formula(sprintf("`%s` ~ `%s`", value, group)), data = df))[[1]]
    key <- paste(res$group2, res$group1, sep = "-")
    hit <- match(key, rownames(tk))
    res$p.adj[!is.na(hit)] <- tk[hit[!is.na(hit)], "p adj"]
  } else {
    res$p <- unname(mapply(function(a, b) {
      x <- df[[value]][df[[group]] == a]
      y <- df[[value]][df[[group]] == b]
      tryCatch(
        if (test == "t.test") stats::t.test(x, y)$p.value else
          suppressWarnings(stats::wilcox.test(x, y)$p.value),
        error = function(e) NA_real_)
    }, as.character(res$group1), as.character(res$group2)))
    res$p.adj <- stats::p.adjust(res$p, method = p.adjust.method)
  }
  res$p.adj.signif <- NULL
  res
}

#' Keep only user-requested pairwise comparisons from a post-hoc table
#' @noRd
.mp_filter_comparisons <- function(res, comparisons) {
  wanted <- vapply(comparisons, function(p) paste(sort(p), collapse = "\r"),
                   character(1))
  have <- mapply(function(a, b) paste(sort(c(a, b)), collapse = "\r"),
                 res$group1, res$group2)
  res[have %in% wanted, , drop = FALSE]
}

#' Convert numeric p-values to conventional significance stars
#' @noRd
.mp_signif_stars <- function(p) {
  cut(p,
      breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 5e-2, Inf),
      labels = c("****", "***", "**", "*", "ns"))
}

#' Omnibus test across all groups, returned as a caption-ready string
#' @noRd
.mp_global_test <- function(df, value, group, test) {
  df <- df[!is.na(df[[value]]) & !is.na(df[[group]]), , drop = FALSE]
  df[[group]] <- as.factor(df[[group]])
  ng <- nlevels(df[[group]])
  if (ng < 2 || stats::sd(df[[value]], na.rm = TRUE) == 0) return(NULL)
  f <- stats::as.formula(sprintf("`%s` ~ `%s`", value, group))
  # stats:: rather than rstatix, which rounds p to three significant digits
  # before it is formatted. Degenerate inputs can still trip the tests; fail
  # soft to NULL.
  tryCatch({
    if (ng == 2) {
      label <- switch(test, "t.test" = "t-test", "Wilcoxon")
      p <- if (test == "t.test") stats::t.test(f, data = df)$p.value else
        suppressWarnings(stats::wilcox.test(f, data = df)$p.value)
      return(sprintf("%s, p = %s", label, format.pval(p, digits = 2)))
    }
    if (test == "anova") {
      p <- summary(stats::aov(f, data = df))[[1]][["Pr(>F)"]][1]
      sprintf("ANOVA, p = %s", format.pval(p, digits = 2))
    } else {
      p <- stats::kruskal.test(f, data = df)$p.value
      sprintf("Kruskal-Wallis, p = %s", format.pval(p, digits = 2))
    }
  }, error = function(e) NULL)
}

#' Run PERMANOVA (adonis2) and return both the model and a caption string
#'
#' @param dist A `dist` object.
#' @param meta Sample metadata data.frame (rows aligned to `dist` labels).
#' @param group Grouping column name.
#' @param permutations Number of permutations.
#' @return List with `$model` (adonis2 result) and `$caption` (string), or
#'   `NULL` if it cannot be computed.
#' @noRd
.mp_permanova <- function(dist, meta, group, permutations = 999) {
  labs <- labels(dist)
  meta <- meta[labs, , drop = FALSE]
  ok <- !is.na(meta[[group]])
  if (sum(ok) < length(labs)) {
    dist <- as.dist(as.matrix(dist)[ok, ok])
    meta <- meta[ok, , drop = FALSE]
  }
  if (length(unique(meta[[group]])) < 2) return(NULL)
  f <- stats::as.formula(sprintf("dist ~ `%s`", group))
  model <- tryCatch(
    vegan::adonis2(f, data = meta, permutations = permutations),
    error = function(e) NULL
  )
  if (is.null(model)) return(NULL)
  r2 <- model$R2[1]
  fstat <- model$F[1]
  p <- model$`Pr(>F)`[1]
  caption <- sprintf(
    "PERMANOVA: R\u00b2 = %.3f, F = %.2f, p = %s (%d permutations)",
    r2, fstat, format.pval(p, digits = 2), permutations
  )
  list(model = model, caption = caption)
}
