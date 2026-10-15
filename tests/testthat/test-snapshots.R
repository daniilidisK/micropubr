# Visual regression tests (vdiffr) ------------------------------------------
# Deterministic figures only: no permutation tests, no random layouts (ggrepel
# placements differ in the last digits across platforms). Review changed
# snapshots with testthat::snapshot_review().

test_that("differential-abundance figures are stable", {
  skip_if_not_installed("vdiffr")
  # R 4.6 changed plotmath metrics, which shifts the log2/log10 axis titles
  # by a fraction of a pixel; the snapshots were recorded with R 4.6. The
  # files are announced first so that skipping does not delete them.
  for (f in c("diff-volcano", "diff-forest", "diff-heatmap", "diff-consensus")) {
    announce_snapshot_file(name = paste0(f, ".svg"))
  }
  skip_if(getRversion() < "4.6.0", "DA snapshots were recorded with R 4.6")
  df <- make_da_df()
  vdiffr::expect_doppelganger("diff-volcano",
                              mp_diff(df, plot.type = "volcano",
                                      label.top = 0))
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
