# Supervised machine learning -------------------------------------------
# Train a random forest or gradient-boosted model on the feature table of a
# phyloseq object and visualise what it learned: which taxa discriminate the
# groups (importance), how well it separates them (ROC / confusion matrix) or
# how well it predicts a continuous outcome (observed vs predicted).

#' Machine-learning classification and regression on microbiome features
#'
#' Trains a random forest (`"rf"`) or gradient-boosted trees (`"xgb"`) model on
#' the (optionally transformed) feature table of a [phyloseq::phyloseq] object,
#' using a sample-metadata column as the response. The task is detected
#' automatically from the response type  - **classification** for a
#' factor/character column, **regression** for a numeric one  - and can be
#' forced with `task`. Out-of-fold predictions from k-fold cross-validation
#' drive the performance metrics; feature importance comes from a model fit on
#' all samples.
#'
#' @param object A [phyloseq::phyloseq] object with counts.
#' @param group.by Sample-metadata column used as the response/outcome.
#' @param method Learner: `"rf"` (randomForest) or `"xgb"` (xgboost).
#' @param task `"auto"` (default), `"classification"` or `"regression"`.
#' @param tax.level Optional rank to agglomerate features to before modelling.
#' @param top.n Optional number of most-abundant features to keep (dimension
#'   reduction); `NULL` keeps them all.
#' @param transform Feature transform applied before modelling: `"clr"`
#'   (centred log-ratio, the default), `"relative"` or `"none"`.
#' @param plot.type What to draw: `"importance"` (ranked feature importance,
#'   the default), `"roc"` (ROC curve, two-class classification only),
#'   `"confusion"` (confusion-matrix heatmap, classification) or `"prediction"`
#'   (observed vs predicted scatter, regression).
#' @param n.features Number of top features shown in the importance plot.
#' @param cv.folds Number of cross-validation folds.
#' @param ntree Number of trees for the random forest.
#' @param nrounds Number of boosting rounds for xgboost.
#' @param seed Optional integer seed for reproducible folds/models.
#' @param palette Palette name or colour vector.
#' @param theme A `ggplot2` theme; defaults to the active theme (see
#'   [mp_set_defaults()], otherwise [mp_theme_publication()]).
#' @param ... Additional arguments passed to the primary geom of the chosen
#'   plot (e.g. [ggplot2::geom_col()] for importance).
#' @param plot.title,plot.subtitle,plot.caption Optional labels.
#' @param font.size Base font size.
#' @param interactive If `TRUE`, return an interactive `plotly` widget instead
#'   of a static `ggplot` (requires `plotly`). See [mp_alpha()].
#' @param return.data If `TRUE`, return the tidy feature-importance data frame
#'   instead of a plot, with the fitted model (`attr(., "model")`), the
#'   cross-validation predictions (`attr(., "cv")`) and a named list of
#'   performance metrics (`attr(., "metrics")`) attached as attributes.
#' @param verbose Print the method and performance summary.
#'
#' @return A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
#'   data frame when `return.data = TRUE`.
#' @export
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' gp <- phyloseq::tax_glom(GlobalPatterns, "Genus")
#' # Classification: which genera discriminate sample types?
#' mp_ml(gp, group.by = "SampleType", method = "rf")
#' # Two-class ROC curve.
#' soil_vs_feces <- phyloseq::subset_samples(
#'   gp, SampleType %in% c("Soil", "Feces"))
#' mp_ml(soil_vs_feces, group.by = "SampleType", plot.type = "roc")
#' }
mp_ml <- function(object,
                  group.by = NULL,
                  method = c("rf", "xgb"),
                  task = c("auto", "classification", "regression"),
                  tax.level = NULL,
                  top.n = NULL,
                  transform = c("clr", "relative", "none"),
                  plot.type = c("importance", "roc", "confusion", "prediction"),
                  n.features = 20,
                  cv.folds = 5,
                  ntree = 500,
                  nrounds = 150,
                  seed = NULL,
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
  object <- .mp_check_phyloseq(object)
  method <- match.arg(method)
  task <- match.arg(task)
  transform <- match.arg(transform)
  plot.type <- match.arg(plot.type)
  group.by <- .mp_check_var(object, group.by, "group.by", required = TRUE)
  theme <- .mp_theme(theme, font.size)
  .mp_require(if (method == "rf") "randomForest" else "xgboost",
              sprintf("mp_ml(method = \"%s\")", method))
  if (!is.null(seed)) set.seed(seed)

  # Build the sample x feature matrix and the aligned response.
  obj <- if (!is.null(tax.level)) .mp_agglomerate(object, tax.level) else object
  if (!is.null(top.n)) {
    obj <- phyloseq::prune_taxa(.mp_top_taxa(obj, top.n), obj)
  }
  X <- .mp_ml_features(obj, transform)
  meta <- .mp_sample_data(obj)
  y <- meta[rownames(X), group.by]
  keep <- !is.na(y)
  X <- X[keep, , drop = FALSE]
  y <- y[keep]
  if (length(y) < 4) {
    rlang::abort("mp_ml() needs at least 4 samples with a non-missing response.",
                 class = "micropubr_ml_too_few")
  }

  # Resolve the learning task from the response type.
  task <- .mp_ml_task(task, y)
  if (task == "classification") {
    y <- factor(y)
    if (nlevels(y) < 2) {
      rlang::abort("Classification needs at least two classes in `group.by`.",
                   class = "micropubr_ml_one_class")
    }
  } else {
    y <- as.numeric(y)
  }
  .mp_ml_check_plot(plot.type, task, y)
  classes <- if (task == "classification") levels(y) else NULL

  if (verbose) {
    rlang::inform(sprintf(
      "mp_ml(): %s %s on %d samples x %d features (%s transform).",
      method, task, nrow(X), ncol(X), transform))
  }

  # Importance from a full-data fit; performance from cross-validation.
  full <- .mp_ml_fit(X, y, task, method, classes, ntree, nrounds)
  imp <- .mp_ml_importance(full, method, colnames(X))
  imp$importance <- imp$importance / max(imp$importance[imp$importance > 0], 1)
  imp <- imp[order(-imp$importance), , drop = FALSE]
  imp <- .mp_relabel_features(imp, obj, tax.level)
  rownames(imp) <- NULL

  cv <- .mp_ml_cv(X, y, task, method, cv.folds, classes, ntree, nrounds)
  metrics <- .mp_ml_metrics(cv, y, task, cv.folds)
  if (verbose) rlang::inform(metrics$caption)

  if (isTRUE(return.data)) {
    attr(imp, "model") <- full
    attr(imp, "cv") <- cv
    attr(imp, "metrics") <- metrics
    attr(imp, "task") <- task
    return(imp)
  }

  caption <- plot.caption %||% metrics$caption
  dots <- list(...)
  p <- switch(
    plot.type,
    importance = .mp_ml_importance_plot(imp, n.features, palette, theme, dots),
    roc        = .mp_ml_roc_plot(cv, y, palette, theme),
    confusion  = .mp_ml_confusion_plot(cv, y, palette, theme, font.size),
    prediction = .mp_ml_prediction_plot(cv, y, palette, theme, dots)
  )
  p <- p + ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                         caption = caption)
  p <- .mp_tag(p, "mp_ml",
               params = list(method = method, task = task,
                             group.by = group.by, tax.level = tax.level,
                             transform = transform, plot.type = plot.type,
                             cv.folds = cv.folds, ntree = ntree,
                             nrounds = nrounds),
               data = .mp_data_info(object, if (task == "classification")
                 group.by),
               packages = c(if (method == "rf") "randomForest" else "xgboost",
                            if (plot.type == "roc") "pROC"))
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

# --- Feature matrix and task resolution ---------------------------------

#' Build the transformed sample x feature matrix for modelling
#' @noRd
.mp_ml_features <- function(object, transform) {
  switch(
    transform,
    clr      = .mp_clr(object),
    relative = .mp_otu_matrix(.mp_to_relative(object), samples_as_rows = TRUE),
    none     = .mp_otu_matrix(object, samples_as_rows = TRUE)
  )
}

#' Decide classification vs regression from the response
#' @noRd
.mp_ml_task <- function(task, y) {
  if (task != "auto") return(task)
  if (is.numeric(y)) "regression" else "classification"
}

#' Validate the requested plot against the resolved task
#' @noRd
.mp_ml_check_plot <- function(plot.type, task, y) {
  if (task == "regression" && plot.type %in% c("roc", "confusion")) {
    rlang::abort(
      sprintf("plot.type = \"%s\" is only available for classification.",
              plot.type),
      class = "micropubr_ml_bad_plot")
  }
  if (task == "classification" && plot.type == "prediction") {
    rlang::abort(
      "plot.type = \"prediction\" is only available for regression.",
      class = "micropubr_ml_bad_plot")
  }
  if (plot.type == "roc" && nlevels(factor(y)) != 2) {
    rlang::abort(
      c("plot.type = \"roc\" needs a two-class response.",
        i = "Subset `group.by` to two groups, or use plot.type = \"confusion\"."),
      class = "micropubr_ml_bad_plot")
  }
  invisible(TRUE)
}

# --- Model fitting and prediction (method-agnostic wrappers) ------------

#' Fit a random forest or xgboost model
#' @noRd
.mp_ml_fit <- function(X, y, task, method, classes, ntree, nrounds) {
  if (method == "rf") {
    yy <- if (task == "classification") droplevels(factor(y, levels = classes))
          else as.numeric(y)
    return(randomForest::randomForest(x = X, y = yy, ntree = ntree,
                                      importance = TRUE))
  }
  label <- if (task == "classification")
    as.integer(factor(y, levels = classes)) - 1L else as.numeric(y)
  dtrain <- xgboost::xgb.DMatrix(data = X, label = label)
  params <- .mp_ml_xgb_params(task, length(classes))
  xgboost::xgb.train(params = params, data = dtrain, nrounds = nrounds,
                     verbose = 0)
}

#' xgboost objective/params for a given task
#' @noRd
.mp_ml_xgb_params <- function(task, k) {
  base <- list(max_depth = 4, eta = 0.3, nthread = 1)
  if (task == "regression") {
    return(c(base, list(objective = "reg:squarederror")))
  }
  if (k == 2) {
    return(c(base, list(objective = "binary:logistic",
                        eval_metric = "logloss")))
  }
  c(base, list(objective = "multi:softprob", num_class = k,
               eval_metric = "mlogloss"))
}

#' Predict from a fitted model; returns class+prob (classification) or value
#' @noRd
.mp_ml_predict <- function(model, X, task, method, classes) {
  if (method == "rf") {
    if (task == "regression") {
      return(list(value = as.numeric(stats::predict(model, X))))
    }
    prob <- stats::predict(model, X, type = "prob")
    cls <- colnames(prob)[max.col(prob, ties.method = "first")]
    return(list(class = cls, prob = prob))
  }
  dnew <- xgboost::xgb.DMatrix(data = X)
  p <- stats::predict(model, dnew)
  if (task == "regression") return(list(value = as.numeric(p)))
  k <- length(classes)
  if (k == 2) {
    prob <- cbind(1 - p, p)
    colnames(prob) <- classes
  } else {
    if (!is.matrix(p)) p <- matrix(p, ncol = k, byrow = TRUE)
    colnames(p) <- classes
    prob <- p
  }
  cls <- colnames(prob)[max.col(prob, ties.method = "first")]
  list(class = cls, prob = prob)
}

#' Tidy feature-importance table from a fitted model
#' @noRd
.mp_ml_importance <- function(model, method, features) {
  if (method == "rf") {
    imp <- randomForest::importance(model)
    col <- intersect(c("MeanDecreaseGini", "MeanDecreaseAccuracy",
                       "%IncMSE", "IncNodePurity"), colnames(imp))[1]
    if (is.na(col)) col <- colnames(imp)[ncol(imp)]
    return(data.frame(feature = rownames(imp),
                      importance = as.numeric(imp[, col]),
                      stringsAsFactors = FALSE))
  }
  gi <- xgboost::xgb.importance(model = model, feature_names = features)
  df <- data.frame(feature = gi$Feature, importance = gi$Gain,
                   stringsAsFactors = FALSE)
  miss <- setdiff(features, df$feature)
  if (length(miss)) {
    df <- rbind(df, data.frame(feature = miss, importance = 0,
                               stringsAsFactors = FALSE))
  }
  df
}

# --- Cross-validation and metrics ---------------------------------------

#' k-fold cross-validation collecting out-of-fold predictions
#' @noRd
.mp_ml_cv <- function(X, y, task, method, folds, classes, ntree, nrounds) {
  n <- nrow(X)
  folds <- max(2L, min(as.integer(folds), n))
  fold_id <- sample(rep_len(seq_len(folds), n))
  oof_value <- rep(NA_real_, n)
  oof_class <- rep(NA_character_, n)
  oof_prob <- if (!is.null(classes)) {
    matrix(NA_real_, n, length(classes), dimnames = list(NULL, classes))
  } else NULL

  for (k in seq_len(folds)) {
    tr <- fold_id != k
    te <- fold_id == k
    if (task == "classification" && nlevels(droplevels(factor(y[tr]))) < 2) {
      next
    }
    model <- .mp_ml_fit(X[tr, , drop = FALSE], y[tr], task, method, classes,
                        ntree, nrounds)
    pr <- .mp_ml_predict(model, X[te, , drop = FALSE], task, method, classes)
    if (task == "regression") {
      oof_value[te] <- pr$value
    } else {
      oof_class[te] <- pr$class
      cols <- intersect(colnames(pr$prob), classes)
      oof_prob[te, cols] <- pr$prob[, cols, drop = FALSE]
    }
  }
  list(value = oof_value, class = oof_class, prob = oof_prob, classes = classes)
}

#' Cohen's kappa for two label vectors
#' @noRd
.mp_cohen_kappa <- function(obs, pred) {
  tab <- table(obs, pred)
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  if (pe == 1) return(NA_real_)
  (po - pe) / (1 - pe)
}

#' Performance metrics + caption from out-of-fold predictions
#' @noRd
.mp_ml_metrics <- function(cv, y, task, folds) {
  if (task == "regression") {
    obs <- as.numeric(y)
    pred <- cv$value
    ok <- !is.na(pred)
    rmse <- sqrt(mean((pred[ok] - obs[ok])^2))
    ss_res <- sum((pred[ok] - obs[ok])^2)
    ss_tot <- sum((obs[ok] - mean(obs[ok]))^2)
    r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
    return(list(
      rmse = rmse, r2 = r2,
      caption = sprintf("%d-fold CV: RMSE = %.3g, R\u00b2 = %.2f",
                        folds, rmse, r2)))
  }
  obs <- factor(y)
  pred <- factor(cv$class, levels = levels(obs))
  ok <- !is.na(pred)
  acc <- mean(pred[ok] == obs[ok])
  kappa <- .mp_cohen_kappa(obs[ok], pred[ok])
  auc <- NA_real_
  if (nlevels(obs) == 2 && !is.null(cv$prob) &&
      requireNamespace("pROC", quietly = TRUE)) {
    pos <- levels(obs)[2]
    auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(
        obs[ok], cv$prob[ok, pos], quiet = TRUE,
        levels = levels(obs), direction = "<"))),
      error = function(e) NA_real_)
  }
  caption <- sprintf("%d-fold CV: accuracy = %.1f%%, \u03ba = %.2f",
                     folds, 100 * acc, kappa)
  if (!is.na(auc)) caption <- sprintf("%s, AUC = %.3f", caption, auc)
  list(accuracy = acc, kappa = kappa, auc = auc,
       confusion = table(Actual = obs[ok], Predicted = pred[ok]),
       caption = caption)
}

# --- Plots --------------------------------------------------------------

#' Ranked feature-importance bar chart
#' @noRd
.mp_ml_importance_plot <- function(imp, n.features, palette, theme, dots) {
  top <- utils::head(imp, n.features)
  top$feature <- factor(top$feature, levels = rev(unique(top$feature)))
  cols <- .mp_resolve_palette(palette, 256, type = "sequential")
  ggplot2::ggplot(top, ggplot2::aes(x = .data$feature, y = .data$importance,
                                    fill = .data$importance)) +
    .mp_geom(ggplot2::geom_col, params = list(width = 0.78), dots = dots) +
    ggplot2::scale_fill_gradientn(colours = cols, guide = "none") +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Relative importance") +
    theme
}

#' ROC curve from out-of-fold probabilities (two-class only)
#' @noRd
.mp_ml_roc_plot <- function(cv, y, palette, theme) {
  .mp_require("pROC", "ROC curves")
  obs <- factor(y)
  pos <- levels(obs)[2]
  ok <- !is.na(cv$prob[, pos])
  roc <- pROC::roc(obs[ok], cv$prob[ok, pos], quiet = TRUE,
                   levels = levels(obs), direction = "<")
  rdf <- data.frame(fpr = 1 - roc$specificities, tpr = roc$sensitivities)
  rdf <- rdf[order(rdf$fpr, rdf$tpr), , drop = FALSE]
  col <- .mp_resolve_palette(palette, 1, type = "qualitative")[1]
  ggplot2::ggplot(rdf, ggplot2::aes(x = .data$fpr, y = .data$tpr)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey70") +
    ggplot2::geom_step(colour = col, linewidth = 0.8, direction = "hv") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "False positive rate", y = "True positive rate") +
    theme
}

#' Confusion-matrix heatmap from out-of-fold predictions
#' @noRd
.mp_ml_confusion_plot <- function(cv, y, palette, theme, font.size) {
  obs <- factor(y)
  pred <- factor(cv$class, levels = levels(obs))
  ok <- !is.na(pred)
  cm <- as.data.frame(table(Actual = obs[ok], Predicted = pred[ok]),
                      stringsAsFactors = FALSE)
  cols <- .mp_resolve_palette(palette, 256, type = "sequential")
  ggplot2::ggplot(cm, ggplot2::aes(x = .data$Predicted, y = .data$Actual,
                                   fill = .data$Freq)) +
    ggplot2::geom_tile(colour = "grey90", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Freq),
                       size = font.size / 3.2) +
    ggplot2::scale_fill_gradientn(colours = cols, name = "Count") +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(x = "Predicted", y = "Actual") +
    theme +
    .mp_rotate_x(levels(obs), font.size)
}

#' Observed vs predicted scatter for regression
#' @noRd
.mp_ml_prediction_plot <- function(cv, y, palette, theme, dots) {
  df <- data.frame(observed = as.numeric(y), predicted = cv$value)
  df <- df[!is.na(df$predicted), , drop = FALSE]
  col <- .mp_resolve_palette(palette, 1, type = "qualitative")[1]
  ggplot2::ggplot(df, ggplot2::aes(x = .data$observed, y = .data$predicted)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey70") +
    .mp_geom(ggplot2::geom_point,
             params = list(size = 2, alpha = 0.8, colour = col),
             dots = dots) +
    ggplot2::labs(x = "Observed", y = "Predicted (out-of-fold)") +
    theme
}
