skip_if_not_installed("phyloseq")

test_that("mp_taxa bar plot returns a ggplot", {
  ps <- make_test_ps()
  p <- mp_taxa(ps, tax.level = "Phylum", group.by = "Group", top.n = 3)
  expect_s3_class(p, "ggplot")
})

test_that("top.n collapses remaining taxa into an Other category", {
  ps <- make_test_ps()
  df <- micropubr:::.mp_melt_taxa(ps, tax.level = "Genus", top.n = 3)
  expect_true("Other" %in% levels(df$.taxon))
  # 3 kept taxa + Other
  expect_lte(nlevels(df$.taxon), 4)
})

test_that("heatmap and bubble variants return ggplots", {
  ps <- make_test_ps()
  expect_s3_class(mp_taxa(ps, tax.level = "Phylum", plot.type = "heatmap"),
                  "ggplot")
  expect_s3_class(mp_taxa(ps, tax.level = "Phylum", plot.type = "bubble"),
                  "ggplot")
})

test_that("cluster ordering does not error", {
  skip_if_not_installed("vegan")
  ps <- make_test_ps()
  expect_s3_class(
    mp_taxa(ps, tax.level = "Genus", cluster.samples = TRUE), "ggplot")
})

test_that("invalid tax.level raises a classed error", {
  ps <- make_test_ps()
  expect_error(mp_taxa(ps, tax.level = "NotARank"),
               class = "micropubr_bad_rank")
})
