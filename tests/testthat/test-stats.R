# Statistical analysis modules -------------------------------------------

test_that("mp_dispersion runs PERMDISP and optional group tests", {
  ps <- make_test_ps()
  p <- mp_dispersion(ps, group.by = "Group")
  expect_s3_class(p, "ggplot")

  p2 <- mp_dispersion(ps, group.by = "Group",
                      add.test = c("anosim", "mrpp", "permanova"),
                      plot.type = "violin")
  expect_s3_class(p2, "ggplot")

  d <- mp_dispersion(ps, group.by = "Group", return.data = TRUE)
  expect_s3_class(d, "data.frame")
  expect_true("distance" %in% names(d))
  expect_false(is.null(attr(d, "betadisper")))
})

test_that("mp_dispersion errors without a usable grouping", {
  ps <- make_test_ps()
  expect_error(mp_dispersion(ps), class = "micropubr_missing_var")
})

test_that("mp_mantel correlates community and environmental distances", {
  ps <- make_test_ps()
  p <- mp_mantel(ps, vars = c("pH", "Age"))
  expect_s3_class(p, "ggplot")

  df <- mp_mantel(ps, vars = c("pH", "Age"), return.data = TRUE)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("community", "explanatory") %in% names(df)))

  pc <- suppressWarnings(
    mp_mantel(ps, vars = c("pH", "Age"), plot.type = "correlogram"))
  expect_s3_class(pc, "ggplot")
})

test_that("mp_mantel needs vars or dist2", {
  ps <- make_test_ps()
  expect_error(mp_mantel(ps), class = "micropubr_mantel_input")
})

test_that("mp_varpart partitions variance among explanatory sets", {
  ps <- make_test_ps()
  p <- mp_varpart(ps, groups = list(Chem = "pH", Host = "Age"))
  expect_s3_class(p, "ggplot")

  df <- mp_varpart(ps, groups = list(Chem = "pH", Host = "Age"),
                   return.data = TRUE)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("fraction", "adj.r2") %in% names(df)))

  expect_error(mp_varpart(ps, groups = list(OnlyOne = "pH")),
               class = "micropubr_varpart_groups")
})

test_that("mp_prevalence draws scatter and core plots", {
  ps <- make_test_ps()
  p <- mp_prevalence(ps, tax.level = "Phylum")
  expect_s3_class(p, "ggplot")

  pc <- mp_prevalence(ps, tax.level = "Phylum", plot.type = "core",
                      prevalence = 0.3)
  expect_s3_class(pc, "ggplot")

  tab <- mp_prevalence(ps, tax.level = "Phylum", return.data = TRUE)
  expect_s3_class(tab, "data.frame")
  expect_true(all(c("prevalence", "mean_abundance") %in% names(tab)))
})

test_that("mp_indicator runs when indicspecies is available", {
  skip_if_not_installed("indicspecies")
  skip_if_not_installed("permute")
  ps <- make_test_ps()
  out <- mp_indicator(ps, group.by = "Group", return.data = TRUE)
  expect_s3_class(out, "data.frame")
})
