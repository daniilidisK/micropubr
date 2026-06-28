skip_if_not_installed("phyloseq")

test_that("mp_alpha returns a ggplot for a single measure", {
  ps <- make_test_ps()
  p <- mp_alpha(ps, measure = "Shannon", group.by = "Group")
  expect_s3_class(p, "ggplot")
})

test_that("mp_alpha combines multiple measures into a patchwork", {
  skip_if_not_installed("patchwork")
  ps <- make_test_ps()
  p <- mp_alpha(ps, measure = c("Shannon", "Observed"), group.by = "Group")
  expect_s3_class(p, "patchwork")
})

test_that("mp_alpha returns a named list when combine = FALSE", {
  ps <- make_test_ps()
  p <- mp_alpha(ps, measure = c("Shannon", "Simpson"), group.by = "Group",
                combine = FALSE)
  expect_type(p, "list")
  expect_named(p, c("Shannon", "Simpson"))
})

test_that("mp_alpha works without grouping", {
  ps <- make_test_ps()
  expect_s3_class(mp_alpha(ps, measure = "Shannon"), "ggplot")
})

test_that("mp_alpha adds statistics when requested", {
  skip_if_not_installed("rstatix")
  skip_if_not_installed("ggpubr")
  ps <- make_test_ps()
  p <- mp_alpha(ps, measure = "Shannon", group.by = "Group",
                test = "wilcoxon", add.stats = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("invalid group.by raises a classed error", {
  ps <- make_test_ps()
  expect_error(mp_alpha(ps, group.by = "Nonexistent"),
               class = "micropubr_bad_var")
})

test_that("non-phyloseq input is rejected", {
  expect_error(mp_alpha(data.frame(x = 1)), class = "micropubr_bad_input")
})
