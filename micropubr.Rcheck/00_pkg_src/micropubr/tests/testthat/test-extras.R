skip_if_not_installed("phyloseq")

test_that("mp_tree builds with a tree and ggtree backend", {
  skip_if_not_installed("ggtree")
  skip_if_not_installed("ape")
  ps <- make_test_ps(tree = TRUE)
  p <- mp_tree(ps, layout = "circular", color.by = "Phylum",
               annotate = "taxonomy")
  expect_s3_class(p, "ggplot")
})

test_that("mp_tree errors without a tree", {
  ps <- make_test_ps(tree = FALSE)
  expect_error(mp_tree(ps), class = "micropubr_missing_tree")
})

test_that("mp_network builds a correlation network", {
  skip_if_not_installed("ggraph")
  skip_if_not_installed("igraph")
  # Build data with genuine co-occurrence: the first six taxa share a latent
  # gradient, guaranteeing some edges pass the thresholds.
  set.seed(42)
  n <- 24; k <- 12
  latent <- stats::rnorm(n)
  otu <- t(vapply(seq_len(k), function(i) {
    stats::rpois(n, lambda = exp(2 + (i <= 6) * 0.9 * latent +
                                   stats::rnorm(n, 0, 0.2)))
  }, numeric(n)))
  rownames(otu) <- paste0("OTU", seq_len(k))
  colnames(otu) <- paste0("S", seq_len(n))
  tax <- cbind(Kingdom = "Bacteria",
               Phylum = sample(c("Firmicutes", "Bacteroidetes"), k, TRUE),
               Genus = paste0("Genus", seq_len(k)))
  rownames(tax) <- rownames(otu)
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax))
  p <- mp_network(ps, method = "spearman", cor.threshold = 0.5,
                  p.threshold = 0.05, node.color.by = "Phylum")
  expect_s3_class(p, "ggplot")
})

test_that("mp_longitudinal spaghetti works on a time-course object", {
  ps <- make_test_ps(nsmp = 12)
  sd <- phyloseq::sample_data(ps)
  sd$Week <- rep(c(0, 1, 2), length.out = phyloseq::nsamples(ps))
  sd$Subject <- rep(paste0("subj", 1:4), length.out = phyloseq::nsamples(ps))
  phyloseq::sample_data(ps) <- sd
  p <- mp_longitudinal(ps, time.var = "Week", subject.var = "Subject",
                       plot.type = "spaghetti", measure = "Shannon",
                       smooth = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("mp_function builds a heatmap from a matrix", {
  m <- matrix(stats::rpois(20 * 10, 30), nrow = 20,
              dimnames = list(paste0("path", 1:20), paste0("S", 1:10)))
  meta <- data.frame(Group = rep(c("A", "B"), each = 5),
                     row.names = colnames(m))
  p <- mp_function(m, metadata = meta, plot.type = "heatmap",
                   group.by = "Group")
  expect_s3_class(p, "ggplot")
})

test_that("time-axis helper coerces numeric-like factors", {
  expect_true(micropubr:::.mp_time_axis(c("1", "2", "3"))$numeric)
  expect_false(micropubr:::.mp_time_axis(c("a", "b"))$numeric)
})
