skip_if_not_installed("phyloseq")
skip_if_not_installed("vegan")

test_that("mp_ordination PCoA returns a ggplot with a PCoA caption", {
  ps <- make_test_ps()
  p <- mp_ordination(ps, method = "PCoA", distance = "bray",
                     group.by = "Group", permanova = FALSE)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$caption, "PCoA", fixed = TRUE)
})

test_that("PERMANOVA result is appended to the caption", {
  ps <- make_test_ps()
  p <- mp_ordination(ps, method = "PCoA", group.by = "Group",
                     permanova = TRUE, permutations = 99)
  expect_match(p$labels$caption, "PERMANOVA")
})

test_that("PCA on CLR data works without a distance", {
  ps <- make_test_ps()
  p <- mp_ordination(ps, method = "PCA", group.by = "Group")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$caption, "Aitchison")
})

test_that("Aitchison distance is accepted for PCoA", {
  ps <- make_test_ps()
  expect_s3_class(
    mp_ordination(ps, method = "PCoA", distance = "aitchison",
                  permanova = FALSE),
    "ggplot")
})

test_that("constrained ordination requires constraints", {
  ps <- make_test_ps()
  expect_error(mp_ordination(ps, method = "RDA"),
               class = "micropubr_missing_constraints")
})
