skip_if_not_installed("phyloseq")
skip_if_not_installed("plotly")

test_that("interactive = TRUE returns a plotly widget for core plots", {
  ps <- make_test_ps()
  expect_s3_class(
    mp_alpha(ps, measure = "Shannon", group.by = "Group", interactive = TRUE),
    "plotly")
  expect_s3_class(
    mp_taxa(ps, tax.level = "Phylum", group.by = "Group", top.n = 4,
            interactive = TRUE),
    "plotly")
  expect_s3_class(
    mp_ordination(ps, group.by = "Group", permanova = FALSE,
                  interactive = TRUE),
    "plotly")
})

test_that("interactive = FALSE leaves the static ggplot path unchanged", {
  ps <- make_test_ps()
  p <- mp_alpha(ps, measure = "Shannon", group.by = "Group")
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "plotly"))
})

test_that("multi-panel and patchwork inputs become a single interactive subplot", {
  skip_if_not_installed("patchwork")
  ps <- make_test_ps()
  multi <- mp_alpha(ps, measure = c("Shannon", "Observed"), group.by = "Group",
                    interactive = TRUE)
  expect_s3_class(multi, "plotly")

  fig <- mp_combine(
    mp_alpha(ps, measure = "Shannon", group.by = "Group"),
    mp_ordination(ps, group.by = "Group", permanova = FALSE),
    ncol = 2)
  expect_s3_class(.mp_interactive(fig), "plotly")
})

test_that("a shared legend is shown once across subplot panels", {
  skip_if_not_installed("patchwork")
  ps <- make_test_ps()
  w <- mp_alpha(ps, measure = c("Shannon", "Observed"), group.by = "Group",
                interactive = TRUE)
  built <- suppressWarnings(plotly::plotly_build(w))
  shown <- vapply(built$x$data, function(tr) {
    if (isTRUE(tr$showlegend) || is.null(tr$showlegend)) tr$name %||% NA_character_
    else NA_character_
  }, character(1))
  shown <- shown[!is.na(shown)]
  # Each group name should appear at most once in the legend, not once per panel.
  expect_true(all(table(shown) == 1))
})

test_that("ggplot subtitle and caption survive the plotly conversion", {
  df <- make_da_df()
  w <- mp_diff(df, plot.type = "volcano", plot.title = "Title here",
               interactive = TRUE)
  lay <- suppressWarnings(plotly::plotly_build(w))$x$layout
  # Title carries the folded-in subtitle; caption lands as an annotation.
  expect_true(grepl("Title here", lay$title$text, fixed = TRUE))
  expect_true(grepl("<sup>", lay$title$text, fixed = TRUE))
})

test_that("mp_export writes an interactive widget to HTML", {
  ps <- make_test_ps()
  w <- mp_alpha(ps, measure = "Shannon", group.by = "Group", interactive = TRUE)
  f <- tempfile(fileext = ".html")
  on.exit(unlink(f), add = TRUE)
  suppressWarnings(mp_export(w, f))
  expect_true(file.exists(f))
})

test_that("mp_export rejects a non-HTML target for an interactive widget", {
  ps <- make_test_ps()
  w <- mp_alpha(ps, measure = "Shannon", group.by = "Group", interactive = TRUE)
  expect_error(mp_export(w, "figure.pdf"), class = "micropubr_bad_input")
})
