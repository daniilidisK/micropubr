skip_if_not_installed("phyloseq")

test_that("mp_combine returns a patchwork from several ggplots", {
  skip_if_not_installed("patchwork")
  ps <- make_test_ps()
  a <- mp_alpha(ps, measure = "Shannon", group.by = "Group")
  b <- mp_taxa(ps, tax.level = "Phylum")
  expect_s3_class(mp_combine(a, b, ncol = 2), "patchwork")
})

test_that("mp_combine accepts a single list of plots", {
  skip_if_not_installed("patchwork")
  ps <- make_test_ps()
  plots <- mp_alpha(ps, measure = c("Shannon", "Observed"),
                    group.by = "Group", combine = FALSE)
  expect_s3_class(mp_combine(plots), "patchwork")
})

test_that("mp_combine rejects non-ggplot inputs", {
  expect_error(mp_combine(1, "two"), class = "micropubr_bad_input")
})

test_that("mp_export writes a file and returns the path invisibly", {
  ps <- make_test_ps()
  p <- mp_taxa(ps, tax.level = "Phylum")
  f <- tempfile(fileext = ".png")
  out <- mp_export(p, f, width = 80, height = 60, units = "mm")
  expect_true(file.exists(f))
  expect_identical(out, f)
})

test_that("mp_export column presets set width without error", {
  ps <- make_test_ps()
  p <- mp_taxa(ps, tax.level = "Phylum")
  f <- tempfile(fileext = ".pdf")
  expect_silent(suppressMessages(mp_export(p, f, column = "single")))
  expect_true(file.exists(f))
})

test_that("mp_export writes PDF and EPS without cairo (macOS sans XQuartz)", {
  local_mocked_bindings(.mp_has_cairo = function() FALSE)
  p <- mp_taxa(make_test_ps(), tax.level = "Phylum")
  for (ext in c(".pdf", ".eps")) {
    f <- tempfile(fileext = ext)
    expect_silent(suppressMessages(mp_export(p, f, column = "single")))
    expect_gt(file.size(f), 0)
  }
})
