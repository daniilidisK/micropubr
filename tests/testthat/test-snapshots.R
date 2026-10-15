# Visual regression tests (vdiffr) ------------------------------------------
# Deterministic figures only: no permutation tests, no random layouts. Review
# changed snapshots with testthat::snapshot_review().

test_that("differential-abundance figures are stable", {
  skip_if_not_installed("vdiffr")
  df <- make_da_df()
  vdiffr::expect_doppelganger("diff-volcano",
                              mp_diff(df, plot.type = "volcano"))
  vdiffr::expect_doppelganger("diff-forest", mp_diff(df, plot.type = "forest"))
  long <- make_da_long()
  vdiffr::expect_doppelganger("diff-heatmap",
                              mp_diff(long, plot.type = "heatmap"))
  vdiffr::expect_doppelganger("diff-consensus",
                              mp_diff(long, plot.type = "consensus"))
})

test_that("core diversity and composition figures are stable", {
  skip_if_not_installed("vdiffr")
  ps <- make_test_ps()
  vdiffr::expect_doppelganger(
    "alpha-box",
    mp_alpha(ps, measure = "Shannon", group.by = "Group", jitter = FALSE))
  vdiffr::expect_doppelganger(
    "ordination-pcoa",
    mp_ordination(ps, group.by = "Group", permanova = FALSE))
  vdiffr::expect_doppelganger(
    "taxa-bar", mp_taxa(ps, tax.level = "Phylum", group.by = "Group"))
})
