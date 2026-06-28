# Supervised machine learning (mp_ml) ------------------------------------

skip_if_not_installed("phyloseq")

test_that("mp_ml classifies and returns a ggplot (random forest)", {
  skip_if_not_installed("randomForest")
  ps <- make_test_ps(nsmp = 16)
  p <- mp_ml(ps, group.by = "Group", method = "rf", cv.folds = 3,
             ntree = 100, seed = 1)
  expect_s3_class(p, "ggplot")
})

test_that("mp_ml classifies with xgboost", {
  skip_if_not_installed("xgboost")
  ps <- make_test_ps(nsmp = 16)
  p <- mp_ml(ps, group.by = "Group", method = "xgb", cv.folds = 3,
             nrounds = 30, seed = 1)
  expect_s3_class(p, "ggplot")
})

test_that("mp_ml return.data yields an importance table with attributes", {
  skip_if_not_installed("randomForest")
  ps <- make_test_ps(nsmp = 16)
  d <- mp_ml(ps, group.by = "Group", method = "rf", cv.folds = 3,
             ntree = 100, seed = 1, return.data = TRUE)
  expect_s3_class(d, "data.frame")
  expect_true(all(c("feature", "importance") %in% names(d)))
  expect_false(is.null(attr(d, "model")))
  expect_false(is.null(attr(d, "metrics")))
  expect_identical(attr(d, "task"), "classification")
})

test_that("mp_ml regresses on a numeric response", {
  skip_if_not_installed("randomForest")
  ps <- make_test_ps(nsmp = 16)
  p <- mp_ml(ps, group.by = "pH", method = "rf", cv.folds = 3,
             ntree = 100, seed = 1, plot.type = "prediction")
  expect_s3_class(p, "ggplot")
  d <- mp_ml(ps, group.by = "pH", method = "rf", cv.folds = 3,
             ntree = 100, seed = 1, return.data = TRUE)
  expect_identical(attr(d, "task"), "regression")
  expect_true(!is.null(attr(d, "metrics")$rmse))
})

test_that("mp_ml draws an ROC curve and a confusion matrix", {
  skip_if_not_installed("randomForest")
  skip_if_not_installed("pROC")
  ps <- make_test_ps(nsmp = 16)
  expect_s3_class(
    mp_ml(ps, group.by = "Group", plot.type = "roc", cv.folds = 3,
          ntree = 100, seed = 1), "ggplot")
  expect_s3_class(
    mp_ml(ps, group.by = "Group", plot.type = "confusion", cv.folds = 3,
          ntree = 100, seed = 1), "ggplot")
})

test_that("mp_ml rejects mismatched plot/task combinations", {
  skip_if_not_installed("randomForest")
  ps <- make_test_ps(nsmp = 16)
  # prediction is regression-only.
  expect_error(
    mp_ml(ps, group.by = "Group", plot.type = "prediction"),
    class = "micropubr_ml_bad_plot")
  # roc needs a two-class response; pH is numeric (regression).
  expect_error(
    mp_ml(ps, group.by = "pH", plot.type = "roc"),
    class = "micropubr_ml_bad_plot")
})
