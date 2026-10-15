# Differential-abundance engine --------------------------------------------
# mp_run_da() runs one or several differential-abundance methods and returns
# one long, tidy table that mp_diff() plots. The design (groups, reference,
# covariates, random effect, prevalence filter) is resolved once and shared by
# every engine, so results from different methods line up row for row and can
# be compared directly -- which is the point: benchmarks show DA methods
# disagree, and agreement between them is the most useful signal.

#' Run a differential-abundance analysis
#'
#' Tests every taxon for a difference in abundance between the levels of
#' `group.by`, with one or several methods, and returns a long, tidy table
#' ready for [mp_diff()]. Running several methods in one call and looking at
#' where they agree (`mp_diff(plot.type = "consensus")`) is the recommended
#' practice: independent benchmarks show that differential-abundance methods
#' frequently disagree, and that adjusting for confounders matters.
#'
#' @section Methods:
#' \describe{
#'   \item{`"wilcox"`}{Wilcoxon rank-sum test on centred log-ratio (CLR)
#'     transformed counts (pseudocount 1). Dependency-free and fast, but it
#'     cannot adjust for covariates. `lfcSE` is the Welch standard error of
#'     the mean CLR difference.}
#'   \item{`"lm"`}{Linear model on CLR-transformed counts with the `adjust`
#'     covariates as additional terms. With `random`, a linear mixed model
#'     with a random intercept per level of `random` (requires `lmerTest`;
#'     Satterthwaite degrees of freedom).}
#'   \item{`"ancombc2"`}{ANCOM-BC2 (package `ANCOMBC`): a log-linear model
#'     with correction for sample-specific sampling fractions, detection of
#'     structural zeros (taxa absent from a whole group) and a pseudo-count
#'     sensitivity analysis, reported in the `robust` column. Supports
#'     `adjust` and `random`; with `pairwise = TRUE` the contrasts control the
#'     mixed directional FDR.}
#'   \item{`"deseq2"`}{DESeq2 negative-binomial Wald test with `"poscounts"`
#'     size factors, which tolerate the many zeros of microbiome data.
#'     Supports `adjust`.}
#'   \item{`"aldex2"`}{ALDEx2 Welch's t-test on Monte-Carlo CLR instances,
#'     run separately for each contrast. Cannot adjust for covariates.}
#' }
#' `"ancombc"` is accepted as an alias of `"ancombc2"`; the original ANCOM-BC
#' has been superseded by its authors.
#'
#' Fold changes are always on the log2 scale and positive when a feature is
#' more abundant in the first-named group of the contrast (`"B vs A"` means
#' higher in B). ANCOM-BC2's natural-log estimates are converted.
#'
#' @param object A [phyloseq::phyloseq] or `(Tree)SummarizedExperiment`
#'   object with counts.
#' @param group.by Sample-metadata column defining the groups.
#' @param method One or more of `"wilcox"`, `"lm"`, `"ancombc2"`, `"deseq2"`
#'   and `"aldex2"`. Several methods are run on the same design and stacked.
#' @param groups Optional subset of `group.by` levels to compare (at least
#'   two). The first is the reference unless `reference` is given, so
#'   `groups = c("Control", "Treatment")` is the classic two-group comparison.
#' @param reference Reference level. Defaults to the first of `groups`, or to
#'   the first level of `group.by`.
#' @param adjust Optional metadata columns to adjust for (covariates or
#'   confounders such as age, sex or medication). Used by `"lm"`,
#'   `"ancombc2"` and `"deseq2"`; the other methods warn and run unadjusted.
#' @param random Optional metadata column identifying repeated measures (for
#'   example a subject id), modelled as a random intercept by `"lm"` and
#'   `"ancombc2"`.
#' @param pairwise If `TRUE`, test every pair of levels; by default each level
#'   is compared with `reference`.
#' @param tax.level Optional rank to agglomerate to first.
#' @param prevalence Minimum fraction of samples a taxon must be detected in
#'   to be tested (default `0.1`, as in ANCOM-BC2).
#' @param p.adjust.method Multiple-testing correction, applied within each
#'   method and contrast (default `"BH"`); see [stats::p.adjust()].
#' @param engine.args Extra arguments for the underlying engine, e.g.
#'   `list(n_cl = 4)` for [ANCOMBC::ancombc2()]. With several methods, nest
#'   them by method: `list(ancombc2 = list(n_cl = 4))`.
#' @param verbose Print the methods used.
#'
#' @return A data frame with one row per feature, contrast and method:
#'   `feature` (display label), `feature_id` (taxon id), `contrast`,
#'   `method`, `log2FC`, `lfcSE`, `pvalue`, `padj`, `mean` (mean relative
#'   abundance), `robust` (passed ANCOM-BC2's pseudo-count sensitivity
#'   analysis; `NA` for other methods) and `structural_zero` (absent from one
#'   group of the contrast, so not tested). The design is recorded in
#'   `attr(x, "mp_da")` and used by [mp_diff()] and [mp_methods()].
#' @export
#' @seealso [mp_diff()] to plot the result.
#' @examples
#' \donttest{
#' data(GlobalPatterns, package = "phyloseq")
#' gp <- phyloseq::subset_samples(GlobalPatterns,
#'                                SampleType %in% c("Feces", "Skin", "Tongue"))
#' # Two groups, one method
#' res <- mp_run_da(gp, group.by = "SampleType", groups = c("Skin", "Feces"),
#'                  tax.level = "Genus")
#' mp_diff(res, plot.type = "volcano")
#'
#' # All groups against a reference, two methods, then their agreement
#' res2 <- mp_run_da(gp, group.by = "SampleType", reference = "Skin",
#'                   method = c("wilcox", "lm"), tax.level = "Phylum")
#' mp_diff(res2, plot.type = "consensus")
#' }
mp_run_da <- function(object,
                      group.by = NULL,
                      method = "wilcox",
                      groups = NULL,
                      reference = NULL,
                      adjust = NULL,
                      random = NULL,
                      pairwise = FALSE,
                      tax.level = NULL,
                      prevalence = 0.1,
                      p.adjust.method = "BH",
                      engine.args = list(),
                      verbose = getOption("micropubr.verbose", FALSE)) {
  object <- .mp_check_phyloseq(object)
  method <- .mp_da_methods(method)
  group.by <- .mp_check_var(object, group.by, "group.by", required = TRUE)
  adjust <- .mp_check_var(object, adjust, "adjust")
  random <- .mp_check_var(object, random, "random")
  if (length(random) > 1) {
    rlang::abort("`random` must name a single metadata column.",
                 class = "micropubr_bad_var")
  }
  if (group.by %in% c(adjust, random)) {
    rlang::abort("`group.by` cannot also appear in `adjust` or `random`.",
                 class = "micropubr_bad_var")
  }
  p.adjust.method <- rlang::arg_match(p.adjust.method, stats::p.adjust.methods)
  # Resolve per-method engine arguments now, so a malformed list fails before
  # any (possibly slow) method has run.
  engine_args <- lapply(stats::setNames(method, method), .mp_engine_args,
                        engine.args = engine.args, methods = method)

  obj <- .mp_agglomerate(object, tax.level)
  design <- .mp_da_design(obj, group.by, groups, reference, adjust, random,
                          pairwise)
  counts <- .mp_da_counts(design$object, prevalence)
  design$mf <- design$mf[colnames(counts), , drop = FALSE]

  if (verbose) {
    rlang::inform(sprintf(
      "mp_run_da(): %s; %s; %d taxa tested.",
      paste(method, collapse = ", "),
      paste(design$contrasts$label, collapse = ", "), nrow(counts)))
  }

  notes <- c(design$notes, attr(counts, "notes"))
  zeros <- NULL
  parts <- lapply(method, function(m) {
    out <- .mp_da_engine(m, counts, design, p.adjust.method,
                         engine_args[[m]])
    notes <<- c(notes, attr(out, "notes"))
    zeros <<- attr(out, "zeros") %||% zeros
    out$method <- rep(m, nrow(out))
    out
  })
  out <- do.call(rbind, parts)

  labels <- .mp_feature_labels(design$object, rownames(counts), tax.level)
  out <- .mp_da_finalise(out, counts, labels)
  notes <- unique(notes)
  attr(out, "mp_da") <- list(
    group.by = group.by,
    levels = design$levels,
    reference = design$reference,
    contrasts = design$contrasts$label,
    pairwise = pairwise,
    adjust = adjust,
    random = random,
    method = method,
    n = design$n,
    n_taxa = nrow(counts),
    n_filtered = attr(counts, "n_filtered"),
    prevalence = prevalence,
    p.adjust.method = p.adjust.method,
    tax.level = tax.level,
    labels = labels,
    zeros = zeros,
    notes = notes,
    versions = .mp_versions(c("micropubr", .mp_da_packages(method, random)))
  )
  .mp_guard_emit(notes)
  out
}

# --- Design ----------------------------------------------------------------

#' Validate `method`, mapping the legacy "ancombc" onto "ancombc2"
#' @noRd
.mp_da_methods <- function(method) {
  choices <- c("wilcox", "lm", "ancombc2", "deseq2", "aldex2", "ancombc")
  if (!is.character(method) || !length(method)) {
    rlang::abort("`method` must be a character vector.",
                 class = "micropubr_bad_method")
  }
  bad <- setdiff(method, choices)
  if (length(bad)) {
    rlang::abort(
      c(sprintf("Unknown `method`: %s.", paste(bad, collapse = ", ")),
        i = sprintf("Choose from: %s.",
                    paste(setdiff(choices, "ancombc"), collapse = ", "))),
      class = "micropubr_bad_method")
  }
  unique(replace(method, method == "ancombc", "ancombc2"))
}

#' Pick the engine.args that belong to method `m`
#' @noRd
.mp_engine_args <- function(m, engine.args, methods) {
  if (!length(engine.args)) return(list())
  if (any(names(engine.args) %in% c(methods, "ancombc"))) {
    return(engine.args[[m]] %||% list())
  }
  if (length(methods) > 1) {
    rlang::abort(
      c("With several methods, nest `engine.args` by method.",
        i = "e.g. engine.args = list(ancombc2 = list(n_cl = 4))."),
      class = "micropubr_bad_engine_args")
  }
  engine.args
}

#' Resolve groups, reference, contrasts and a model frame with safe names
#'
#' Levels and columns are recoded to tokens (`L1`, `mp_cov1`, ...) so that any
#' level or column name survives model.matrix() and ANCOM-BC2's formula
#' strings; results are mapped back to the original names afterwards.
#' @noRd
.mp_da_design <- function(object, group.by, groups, reference, adjust,
                          random, pairwise) {
  meta <- .mp_sample_data(object)
  g <- meta[[group.by]]
  lv_all <- if (is.factor(g)) levels(droplevels(g)) else
    sort(unique(as.character(g[!is.na(g)])))

  lv <- lv_all
  if (!is.null(groups)) {
    groups <- as.character(groups)
    if (length(groups) < 2 || anyDuplicated(groups) ||
        !all(groups %in% lv_all)) {
      rlang::abort(
        c("`groups` must name at least two distinct levels of `group.by`.",
          i = sprintf("Levels: %s.", paste(lv_all, collapse = ", "))),
        class = "micropubr_bad_groups")
    }
    lv <- groups
  }
  reference <- reference %||% lv[1]
  if (!reference %in% lv) {
    rlang::abort(
      sprintf("`reference = \"%s\"` is not among the compared levels: %s.",
              reference, paste(lv, collapse = ", ")),
      class = "micropubr_bad_groups")
  }
  lv <- c(reference, setdiff(lv, reference))

  vars <- c(group.by, adjust, random)
  in_lv <- as.character(g) %in% lv
  complete <- stats::complete.cases(meta[, vars, drop = FALSE])
  notes <- character(0)
  if (any(in_lv & !complete)) {
    notes <- c(notes, sprintf(
      "%d sample(s) with missing values in %s were excluded.",
      sum(in_lv & !complete), paste(vars, collapse = ", ")))
  }
  meta <- meta[in_lv & complete, , drop = FALSE]

  n <- table(factor(meta[[group.by]], levels = lv))
  if (n[[reference]] < 2) {
    rlang::abort(
      sprintf("The reference level \"%s\" needs at least two samples.",
              reference),
      class = "micropubr_need_groups")
  }
  tiny <- names(n)[n < 2]
  if (length(tiny)) {
    notes <- c(notes, sprintf(
      "Level(s) %s dropped: fewer than two samples.",
      paste(tiny, collapse = ", ")))
    lv <- setdiff(lv, tiny)
    meta <- meta[meta[[group.by]] %in% lv, , drop = FALSE]
    n <- n[lv]
  }
  if (length(lv) < 2) {
    rlang::abort(
      "`group.by` needs at least two levels with two or more samples each.",
      class = "micropubr_need_groups")
  }
  low <- names(n)[n < 5]
  if (length(low)) {
    notes <- c(notes, sprintf(paste(
      "Group(s) %s have fewer than 5 samples; differential-abundance tests",
      "have little power and unstable variance estimates at this size."),
      paste(sprintf("%s (n = %d)", low, as.integer(n[low])),
            collapse = ", ")))
  }

  tokens <- stats::setNames(paste0("L", seq_along(lv)), lv)
  mf <- data.frame(
    mp_group = factor(unname(tokens[as.character(meta[[group.by]])]),
                      levels = unname(tokens)),
    row.names = rownames(meta))
  cov_names <- character(0)
  for (i in seq_along(adjust)) {
    v <- meta[[adjust[i]]]
    if (is.character(v) || is.logical(v)) v <- factor(v)
    nm <- paste0("mp_cov", i)
    mf[[nm]] <- v
    cov_names <- c(cov_names, nm)
  }
  if (!is.null(random)) mf$mp_rand <- factor(meta[[random]])
  notes <- c(notes, .mp_guard_collinear(mf, adjust))

  contrasts <- if (isTRUE(pairwise)) {
    pr <- utils::combn(seq_along(lv), 2)
    data.frame(level = lv[pr[2, ]], ref = lv[pr[1, ]],
               stringsAsFactors = FALSE)
  } else {
    data.frame(level = lv[-1], ref = lv[1], stringsAsFactors = FALSE)
  }
  contrasts$label <- sprintf("%s vs %s", contrasts$level, contrasts$ref)

  list(object = phyloseq::prune_samples(rownames(meta), object),
       mf = mf, group.by = group.by, levels = lv, tokens = tokens,
       reference = lv[1], adjust = adjust, cov_names = cov_names,
       random = random, contrasts = contrasts, pairwise = isTRUE(pairwise),
       n = stats::setNames(as.integer(n), names(n)), notes = notes)
}

#' Covariates so strongly tied to the grouping that adjusting for them is
#' unstable
#'
#' Cramer's V for a categorical covariate, R-squared on the group for a numeric
#' one. In a confounding simulation on real counts, ANCOM-BC2 adjusted for a
#' covariate that matched the group in 90% of samples (V ~ 0.8) produced many
#' false positives; the linear model's standard errors widen instead. Either
#' way, the data cannot separate the two effects well.
#' @noRd
.mp_guard_collinear <- function(mf, adjust) {
  if (!length(adjust)) return(character(0))
  g <- mf$mp_group
  out <- character(0)
  for (i in seq_along(adjust)) {
    v <- mf[[paste0("mp_cov", i)]]
    strength <- if (is.factor(v)) {
      tab <- table(g, droplevels(v))
      k <- min(dim(tab)) - 1
      if (k < 1) next
      chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
      c(value = sqrt(unname(chi) / (sum(tab) * k)), cut = 0.7)
    } else {
      r2 <- summary(stats::lm(v ~ g))$r.squared
      c(value = r2, cut = 0.5)
    }
    if (is.finite(strength[["value"]]) &&
        strength[["value"]] >= strength[["cut"]]) {
      out <- c(out, sprintf(paste(
        "Covariate %s is strongly associated with the grouping (%s = %.2f);",
        "the data can barely separate their effects, so adjusted estimates are",
        "unstable and adjusted ANCOM-BC2 in particular can produce false",
        "positives."), adjust[i],
        if (is.factor(v)) "Cramer's V" else "R-squared",
        strength[["value"]]))
    }
  }
  out
}

#' Count matrix (taxa x samples) after the prevalence filter
#' @noRd
.mp_da_counts <- function(object, prevalence) {
  counts <- .mp_otu_matrix(object, samples_as_rows = FALSE)
  keep <- rowMeans(counts > 0) >= prevalence & rowSums(counts) > 0
  if (sum(keep) < 2) {
    rlang::abort(
      c(sprintf("Fewer than two taxa are detected in >= %g%% of samples.",
                100 * prevalence),
        i = "Lower `prevalence` or agglomerate less (`tax.level`)."),
      class = "micropubr_da_empty")
  }
  notes <- character(0)
  if (any(abs(counts - round(counts)) > 1e-8)) {
    notes <- paste(
      "The abundance table is not integer counts (relative abundances?):",
      "DESeq2 rounds them and the CLR pseudocount of 1 assumes counts.")
  }
  out <- counts[keep, , drop = FALSE]
  attr(out, "n_filtered") <- sum(!keep)
  attr(out, "notes") <- notes
  out
}

#' CLR transform of a taxa x samples count matrix, returned samples x taxa
#' @noRd
.mp_clr_counts <- function(counts, pseudocount = 1) {
  logm <- log(t(counts) + pseudocount)
  sweep(logm, 1, rowMeans(logm), "-")
}

# --- Engines ---------------------------------------------------------------

#' Dispatch one method; every engine returns feature_id, contrast, log2FC,
#' lfcSE, pvalue, padj, robust and structural_zero
#' @noRd
.mp_da_engine <- function(m, counts, design, p.adjust.method, args) {
  switch(
    m,
    wilcox   = .mp_da_wilcox(counts, design, p.adjust.method),
    lm       = .mp_da_lm(counts, design, p.adjust.method),
    ancombc2 = .mp_da_ancombc2(counts, design, p.adjust.method, args),
    deseq2   = .mp_da_deseq2(counts, design, p.adjust.method, args),
    aldex2   = .mp_da_aldex2(counts, design, p.adjust.method, args)
  )
}

#' Assemble one engine's result for one contrast
#' @noRd
.mp_da_rows <- function(feature_id, contrast, log2FC, lfcSE, pvalue, padj,
                        robust = NA, structural_zero = FALSE) {
  n <- length(feature_id)
  data.frame(feature_id = as.character(feature_id),
             contrast = rep(contrast, n),
             log2FC = as.numeric(log2FC),
             lfcSE = rep_len(as.numeric(lfcSE), n),
             pvalue = as.numeric(pvalue),
             padj = as.numeric(padj),
             robust = rep_len(as.logical(robust), n),
             structural_zero = rep_len(as.logical(structural_zero), n),
             stringsAsFactors = FALSE)
}

#' Note for methods that cannot use `adjust`/`random`
#' @noRd
.mp_da_unadjusted_note <- function(design, label) {
  if (!length(design$adjust) && is.null(design$random)) return(character(0))
  sprintf(paste("%s cannot adjust for covariates or repeated measures;",
                "`adjust`/`random` were ignored and its results are",
                "unadjusted."), label)
}

#' Wilcoxon rank-sum test on CLR values
#' @noRd
.mp_da_wilcox <- function(counts, design, p.adjust.method) {
  clr <- .mp_clr_counts(counts)
  g <- design$mf$mp_group
  tok <- design$tokens
  out <- lapply(seq_len(nrow(design$contrasts)), function(i) {
    ct <- design$contrasts[i, ]
    a <- clr[g == tok[[ct$ref]], , drop = FALSE]
    b <- clr[g == tok[[ct$level]], , drop = FALSE]
    d <- colMeans(b) - colMeans(a)
    se <- sqrt(apply(b, 2, stats::var) / nrow(b) +
                 apply(a, 2, stats::var) / nrow(a))
    pv <- vapply(seq_len(ncol(clr)), function(j) tryCatch(
      suppressWarnings(stats::wilcox.test(b[, j], a[, j])$p.value),
      error = function(e) NA_real_), numeric(1))
    # natural-log CLR difference -> log2 fold change
    .mp_da_rows(colnames(clr), ct$label, d / log(2), se / log(2), pv,
                stats::p.adjust(pv, p.adjust.method))
  })
  out <- do.call(rbind, out)
  attr(out, "notes") <- .mp_da_unadjusted_note(design, "The Wilcoxon test")
  out
}

#' Contrast vector over model coefficients for `level` vs `ref`
#' @noRd
.mp_contrast_vector <- function(coef_names, design, level, ref) {
  L <- stats::setNames(numeric(length(coef_names)), coef_names)
  tok <- design$tokens
  if (level != design$reference) L[paste0("mp_group", tok[[level]])] <- 1
  if (ref != design$reference) L[paste0("mp_group", tok[[ref]])] <- -1
  L
}

#' Linear (mixed) model on CLR values
#' @noRd
.mp_da_lm <- function(counts, design, p.adjust.method) {
  clr <- .mp_clr_counts(counts)
  rhs <- paste(c("mp_group", design$cov_names), collapse = " + ")
  if (!is.null(design$random)) {
    return(.mp_da_lmer(clr, design, rhs, p.adjust.method))
  }
  X <- stats::model.matrix(stats::as.formula(paste("~", rhs)),
                           data = design$mf)
  if (qr(X)$rank < ncol(X)) {
    rlang::abort(
      c("The design is rank-deficient.",
        i = "A covariate in `adjust` is confounded with `group.by` or with another covariate."),
      class = "micropubr_da_design")
  }
  df <- nrow(X) - ncol(X)
  if (df < 1) {
    rlang::abort("Too few samples for the number of model terms.",
                 class = "micropubr_da_design")
  }
  # One least-squares solve for every taxon at once.
  XtXi <- chol2inv(chol(crossprod(X)))
  coef <- XtXi %*% crossprod(X, clr)
  rownames(coef) <- colnames(X)
  sigma2 <- colSums((clr - X %*% coef)^2) / df

  out <- lapply(seq_len(nrow(design$contrasts)), function(i) {
    ct <- design$contrasts[i, ]
    L <- .mp_contrast_vector(colnames(X), design, ct$level, ct$ref)
    est <- drop(crossprod(L, coef))
    se <- sqrt(drop(crossprod(L, XtXi %*% L)) * sigma2)
    pv <- 2 * stats::pt(-abs(est / se), df)
    .mp_da_rows(colnames(clr), ct$label, est / log(2), se / log(2), pv,
                stats::p.adjust(pv, p.adjust.method))
  })
  do.call(rbind, out)
}

#' Linear mixed model per taxon (random intercept), Satterthwaite df
#' @noRd
.mp_da_lmer <- function(clr, design, rhs, p.adjust.method) {
  .mp_require("lmerTest", "Linear mixed models (`random`)")
  form <- stats::as.formula(paste("y ~", rhs, "+ (1 | mp_rand)"))
  d <- design$mf
  nct <- nrow(design$contrasts)
  est <- se <- pv <- matrix(NA_real_, ncol(clr), nct)
  for (j in seq_len(ncol(clr))) {
    d$y <- clr[, j]
    fit <- tryCatch(
      suppressMessages(suppressWarnings(lmerTest::lmer(form, data = d))),
      error = function(e) NULL)
    if (is.null(fit)) next
    cn <- names(lme4::fixef(fit))
    L <- do.call(rbind, lapply(seq_len(nct), function(i) {
      .mp_contrast_vector(cn, design, design$contrasts$level[i],
                          design$contrasts$ref[i])
    }))
    ctab <- tryCatch(lmerTest::contest(fit, L = L, joint = FALSE),
                     error = function(e) NULL)
    if (is.null(ctab)) next
    est[j, ] <- ctab$Estimate
    se[j, ] <- ctab[["Std. Error"]]
    pv[j, ] <- ctab[["Pr(>|t|)"]]
  }
  out <- lapply(seq_len(nct), function(i) {
    .mp_da_rows(colnames(clr), design$contrasts$label[i], est[, i] / log(2),
                se[, i] / log(2), pv[, i],
                stats::p.adjust(pv[, i], p.adjust.method))
  })
  do.call(rbind, out)
}

#' ANCOM-BC2, called with a count matrix so no extra converters are needed
#' @noRd
.mp_da_ancombc2 <- function(counts, design, p.adjust.method, args) {
  .mp_require("ANCOMBC", "ANCOM-BC2 differential abundance", bioc = TRUE)
  nlev <- nlevels(design$mf$mp_group)
  pair <- design$pairwise && nlev > 2
  # ANCOMBC drops a one-column metadata table to a vector; the sample id keeps
  # it a data frame whatever the design.
  meta <- design$mf
  meta$mp_sample <- rownames(meta)
  base <- list(
    data = counts, taxa_are_rows = TRUE, meta_data = meta,
    fix_formula = paste(c("mp_group", design$cov_names), collapse = " + "),
    rand_formula = if (!is.null(design$random)) "(1 | mp_rand)" else NULL,
    p_adj_method = p.adjust.method, pseudo_sens = TRUE,
    # Prevalence was already filtered, identically for every method.
    prv_cut = 0, lib_cut = 0, group = "mp_group", struc_zero = TRUE,
    alpha = 0.05, verbose = FALSE, global = FALSE, pairwise = pair,
    dunnet = FALSE, trend = FALSE)
  call_args <- utils::modifyList(base, args, keep.null = TRUE)

  warns <- character(0)
  fit <- NULL
  utils::capture.output(fit <- withCallingHandlers(
    suppressMessages(do.call(ANCOMBC::ancombc2, call_args)),
    warning = function(w) {
      warns <<- c(warns, gsub("\\s+", " ", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }))
  # Our own small-group note already says this more specifically.
  warns <- warns[!grepl("Small sample size|NaNs produced|multi-group comparisons", warns)]

  res <- if (isTRUE(call_args$pairwise)) fit$res_pair else fit$res
  zeros <- .mp_ancombc_zero_matrix(fit$zero_ind, design)
  tag <- paste0("mp_group", design$tokens)
  names(tag) <- names(design$tokens)
  col <- function(prefix, t) {
    nm <- paste0(prefix, t)
    if (nm %in% names(res)) res[[nm]] else rep(NA, nrow(res))
  }

  out <- lapply(seq_len(nrow(design$contrasts)), function(i) {
    ct <- design$contrasts[i, ]
    t <- if (ct$ref == design$reference) tag[[ct$level]] else
      paste0(tag[[ct$level]], "_", tag[[ct$ref]])
    rows <- .mp_da_rows(res$taxon, ct$label, col("lfc_", t) / log(2),
                        col("se_", t) / log(2), col("p_", t), col("q_", t),
                        robust = col("passed_ss_", t))
    # ANCOM-BC2 leaves structural zeros out of `res`; report them explicitly
    # for the contrasts in which a taxon is present in one group only.
    if (!is.null(zeros)) {
      z <- xor(zeros[, ct$level], zeros[, ct$ref])
      ids <- setdiff(rownames(zeros)[z], rows$feature_id)
      if (length(ids)) {
        rows <- rbind(rows, .mp_da_rows(ids, ct$label, NA, NA, NA, NA,
                                        structural_zero = TRUE))
      }
    }
    rows
  })
  out <- do.call(rbind, out)
  attr(out, "zeros") <- zeros
  attr(out, "notes") <- unique(warns)
  out
}

#' ANCOM-BC2's zero_ind table as a logical taxa x level matrix
#' @noRd
.mp_ancombc_zero_matrix <- function(zero_ind, design) {
  if (is.null(zero_ind) || !nrow(zero_ind)) return(NULL)
  cols <- grep("^structural_zero \\(", names(zero_ind), value = TRUE)
  tok <- sub("^structural_zero \\(mp_group = (.+)\\)$", "\\1", cols)
  lev <- names(design$tokens)[match(tok, design$tokens)]
  if (anyNA(lev)) return(NULL)
  m <- as.matrix(zero_ind[, cols, drop = FALSE])
  dimnames(m) <- list(as.character(zero_ind$taxon), lev)
  m[, design$levels[design$levels %in% lev], drop = FALSE]
}

#' DESeq2 Wald test with poscounts size factors
#' @noRd
.mp_da_deseq2 <- function(counts, design, p.adjust.method, args) {
  .mp_require("DESeq2", "DESeq2 differential abundance", bioc = TRUE)
  notes <- if (!is.null(design$random)) {
    "DESeq2 cannot model repeated measures; `random` was ignored."
  } else character(0)
  mf <- design$mf
  mf$mp_rand <- NULL
  cnt <- round(counts)
  storage.mode(cnt) <- "integer"
  f <- stats::as.formula(paste(
    "~", paste(c(design$cov_names, "mp_group"), collapse = " + ")))
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(
    countData = cnt, colData = mf, design = f))
  dds <- suppressMessages(do.call(
    DESeq2::DESeq,
    utils::modifyList(list(object = dds, sfType = "poscounts",
                           quiet = TRUE), args)))
  tok <- design$tokens
  out <- lapply(seq_len(nrow(design$contrasts)), function(i) {
    ct <- design$contrasts[i, ]
    r <- DESeq2::results(dds, contrast = c("mp_group", tok[[ct$level]],
                                           tok[[ct$ref]]),
                         pAdjustMethod = p.adjust.method)
    .mp_da_rows(rownames(r), ct$label, r$log2FoldChange, r$lfcSE,
                r$pvalue, r$padj)
  })
  out <- do.call(rbind, out)
  attr(out, "notes") <- notes
  out
}

#' ALDEx2 Welch's t-test, one two-group run per contrast
#' @noRd
.mp_da_aldex2 <- function(counts, design, p.adjust.method, args) {
  .mp_require("ALDEx2", "ALDEx2 differential abundance", bioc = TRUE)
  g <- design$mf$mp_group
  tok <- design$tokens
  out <- lapply(seq_len(nrow(design$contrasts)), function(i) {
    ct <- design$contrasts[i, ]
    idx <- g %in% c(tok[[ct$ref]], tok[[ct$level]])
    # ALDEx2 orders conditions alphabetically; these labels make
    # diff.btw = comparison - reference.
    conds <- ifelse(g[idx] == tok[[ct$ref]], "a_reference", "b_comparison")
    ax <- suppressMessages(do.call(
      ALDEx2::aldex,
      utils::modifyList(list(reads = round(counts[, idx, drop = FALSE]),
                             conditions = conds, test = "t", effect = TRUE,
                             denom = "all", verbose = FALSE), args)))
    padj <- if (p.adjust.method == "BH") ax$we.eBH else
      stats::p.adjust(ax$we.ep, p.adjust.method)
    # ALDEx2 works in log2 CLR space, so diff.btw is already a log2 ratio.
    .mp_da_rows(rownames(ax), ct$label, ax$diff.btw, NA, ax$we.ep, padj)
  })
  out <- do.call(rbind, out)
  attr(out, "notes") <- .mp_da_unadjusted_note(design, "ALDEx2 (t-test)")
  out
}

# --- Output ----------------------------------------------------------------

#' Add labels and mean abundance, order columns, drop untestable rows
#' @noRd
.mp_da_finalise <- function(out, counts, labels) {
  out <- out[!is.na(out$log2FC) | out$structural_zero, , drop = FALSE]
  rel <- sweep(counts, 2, pmax(colSums(counts), .Machine$double.eps), "/")
  mean_ab <- rowMeans(rel)
  out$mean <- unname(mean_ab[out$feature_id])
  out$feature <- unname(labels[out$feature_id])
  out$feature[is.na(out$feature)] <- out$feature_id[is.na(out$feature)]
  out <- out[, c("feature", "feature_id", "contrast", "method", "log2FC",
                 "lfcSE", "pvalue", "padj", "mean", "robust",
                 "structural_zero")]
  rownames(out) <- NULL
  out
}

#' Display labels for taxon ids at `tax.level`, made unique
#'
#' Taxa unassigned at `tax.level` are named by their deepest assigned rank, as
#' in the composition plots ("Unassigned Muribaculaceae"). Duplicated names
#' (several "uncultured" genera, say) get their id appended so that labels stay
#' unambiguous on a plot.
#' @noRd
.mp_feature_labels <- function(object, ids, tax.level) {
  out <- stats::setNames(ids, ids)
  if (is.null(tax.level)) return(out)
  tt <- tryCatch(as.data.frame(unclass(phyloseq::tax_table(object)),
                               stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(tt) || !tax.level %in% names(tt)) return(out)
  lab <- .mp_unassigned_labels(tt[ids, , drop = FALSE], tax.level,
                               phyloseq::rank_names(object))
  dup <- lab %in% lab[duplicated(lab)]
  lab[dup] <- sprintf("%s (%s)", lab[dup], ids[dup])
  stats::setNames(lab, ids)
}

#' Packages whose versions describe a DA run
#' @noRd
.mp_da_packages <- function(method, random) {
  pk <- c(wilcox = "stats", lm = "stats", ancombc2 = "ANCOMBC",
          deseq2 = "DESeq2", aldex2 = "ALDEx2")[method]
  if (!is.null(random) && "lm" %in% method) pk <- c(pk, "lmerTest", "lme4")
  unique(unname(pk))
}

#' Installed versions of `pkgs` (NA when absent)
#' @noRd
.mp_versions <- function(pkgs) {
  vapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) {
      as.character(utils::packageVersion(p))
    } else NA_character_
  }, character(1))
}
