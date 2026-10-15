# TreeSummarizedExperiment input --------------------------------------------

make_test_tse <- function(tree = TRUE) {
  ps <- make_test_ps(tree = tree)
  args <- list(
    assays = list(counts = as(phyloseq::otu_table(ps), "matrix")),
    rowData = S4Vectors::DataFrame(as(phyloseq::tax_table(ps), "matrix")),
    colData = S4Vectors::DataFrame(
      as(phyloseq::sample_data(ps), "data.frame")))
  if (tree) args$rowTree <- phyloseq::phy_tree(ps)
  list(ps = ps,
       tse = do.call(TreeSummarizedExperiment::TreeSummarizedExperiment, args))
}

test_that("mp_as_phyloseq round-trips counts, taxonomy, metadata and tree", {
  skip_if_not_installed("TreeSummarizedExperiment")
  skip_if_not_installed("ape")
  x <- make_test_tse()
  ps2 <- mp_as_phyloseq(x$tse)
  expect_s4_class(ps2, "phyloseq")
  expect_equal(as(phyloseq::otu_table(ps2), "matrix"),
               as(phyloseq::otu_table(x$ps), "matrix"))
  expect_equal(phyloseq::rank_names(ps2), phyloseq::rank_names(x$ps))
  expect_equal(phyloseq::sample_data(ps2)$Group,
               phyloseq::sample_data(x$ps)$Group)
  expect_setequal(phyloseq::phy_tree(ps2)$tip.label,
                  phyloseq::taxa_names(x$ps))
  expect_identical(mp_as_phyloseq(x$ps), x$ps)
})

test_that("every mp_*() function accepts a TreeSummarizedExperiment", {
  skip_if_not_installed("TreeSummarizedExperiment")
  x <- make_test_tse(tree = FALSE)
  expect_s3_class(mp_alpha(x$tse, measure = "Shannon", group.by = "Group"),
                  "ggplot")
  expect_s3_class(mp_ordination(x$tse, group.by = "Group",
                                permutations = 99), "ggplot")
  res <- mp_run_da(x$tse, group.by = "Group")
  expect_equal(res, mp_run_da(x$ps, group.by = "Group"),
               ignore_attr = TRUE)
})

test_that("unsupported inputs and assays are reported", {
  expect_error(mp_as_phyloseq(data.frame(a = 1)),
               class = "micropubr_bad_input")
  skip_if_not_installed("TreeSummarizedExperiment")
  x <- make_test_tse(tree = FALSE)
  expect_error(mp_as_phyloseq(x$tse, assay.type = "relabundance"),
               class = "micropubr_bad_input")
  expect_error(mp_alpha(data.frame(a = 1)), class = "micropubr_bad_input")
})
