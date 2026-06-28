# Customizability layer: defaults, return.data and ... passthrough --------

test_that("mp_set_defaults / get / reset round-trip", {
  on.exit(mp_reset_defaults(), add = TRUE)
  mp_set_defaults(font.size = 9, palette = "okabe", verbose = TRUE)
  d <- mp_get_defaults()
  expect_equal(d$font.size, 9)
  expect_equal(d$palette, "okabe")
  expect_true(d$verbose)

  mp_reset_defaults()
  d2 <- mp_get_defaults()
  expect_equal(d2$font.size, 11)
  expect_null(d2$palette)
  expect_false(d2$verbose)
})

test_that("mp_set_defaults validates its arguments", {
  on.exit(mp_reset_defaults(), add = TRUE)
  expect_error(mp_set_defaults(font.size = "big"),
               class = "micropubr_bad_default")
  expect_error(mp_set_defaults(theme = "not a theme"),
               class = "micropubr_bad_default")
})

test_that("a global default palette is honoured by plots", {
  on.exit(mp_reset_defaults(), add = TRUE)
  ps <- make_test_ps()
  mp_set_defaults(palette = "okabe")
  expect_s3_class(mp_alpha(ps, measure = "Shannon", group.by = "Group"),
                  "ggplot")
})

test_that("return.data returns tidy data frames, not plots", {
  ps <- make_test_ps()

  a <- mp_alpha(ps, measure = c("Shannon", "Observed"), group.by = "Group",
                return.data = TRUE)
  expect_s3_class(a, "data.frame")
  expect_true(all(c("Shannon", "Observed") %in% names(a)))

  o <- mp_ordination(ps, method = "PCoA", group.by = "Group",
                     return.data = TRUE)
  expect_s3_class(o, "data.frame")
  expect_false(is.null(attr(o, "permanova")))

  tx <- mp_taxa(ps, tax.level = "Phylum", return.data = TRUE)
  expect_s3_class(tx, "data.frame")
  expect_true(".taxon" %in% names(tx))

  dd <- mp_diff(make_da_df(), return.data = TRUE)
  expect_s3_class(dd, "data.frame")
  expect_true("direction" %in% names(dd))
})

test_that("... is forwarded to the primary geom without error", {
  ps <- make_test_ps()
  expect_s3_class(
    mp_ordination(ps, method = "PCoA", group.by = "Group", stroke = 1.2),
    "ggplot")
  expect_s3_class(
    mp_alpha(ps, measure = "Shannon", group.by = "Group", linewidth = 0.8),
    "ggplot")
  expect_s3_class(
    mp_taxa(ps, tax.level = "Phylum", colour = "white"),
    "ggplot")
})
