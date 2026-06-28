test_that("mp_diff parses a generic data frame into a volcano plot", {
  df <- make_da_df()
  p <- mp_diff(df, plot.type = "volcano")
  expect_s3_class(p, "ggplot")
})

test_that("direction classification flags planted up/down features", {
  df <- make_da_df()
  parsed <- micropubr:::.mp_parse_da(df)
  dir <- micropubr:::.mp_da_direction(parsed, fc.threshold = 1,
                                      p.threshold = 0.05)
  expect_setequal(levels(dir), c("Enriched", "Depleted", "NS"))
  expect_gt(sum(dir != "NS"), 0)
})

test_that("MA, forest and dot plots build for data with signal", {
  df <- make_da_df()
  expect_s3_class(mp_diff(df, plot.type = "ma"), "ggplot")
  expect_s3_class(mp_diff(df, plot.type = "forest"), "ggplot")
  expect_s3_class(mp_diff(df, plot.type = "dot"), "ggplot")
})

test_that("forest plot errors cleanly when nothing is significant", {
  df <- data.frame(feature = paste0("f", 1:5), log2FC = rep(0.1, 5),
                   pvalue = rep(0.9, 5), padj = rep(0.99, 5))
  expect_error(mp_diff(df, plot.type = "forest"),
               class = "micropubr_da_empty")
})

test_that("concordance plot requires a second result", {
  df <- make_da_df()
  expect_error(mp_diff(df, plot.type = "concordance"),
               class = "micropubr_da_no_result2")
  expect_s3_class(
    mp_diff(df, plot.type = "concordance", result2 = make_da_df(seed = 2)),
    "ggplot")
})

test_that("missing fold-change column is reported", {
  df <- data.frame(feature = "a", pvalue = 0.1, padj = 0.2)
  expect_error(micropubr:::.mp_parse_da(df), class = "micropubr_da_parse")
})
