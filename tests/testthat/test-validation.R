# Statistical equivalence with the reference implementations -----------------
# The numbers micropubr prints in captions and returns in tables must be the
# numbers the underlying packages compute. These tests call both sides on the
# same data (with identical seeds for permutation tests) and require equality.

test_that("PERMANOVA in mp_ordination equals vegan::adonis2", {
  ps <- make_da_ps(n_per = 8)
  set.seed(42)
  df <- mp_ordination(ps, group.by = "Group", return.data = TRUE)
  ours <- attr(df, "permanova")

  d <- phyloseq::distance(ps, method = "bray")
  meta <- data.frame(phyloseq::sample_data(ps))[labels(d), , drop = FALSE]
  set.seed(42)
  ref <- vegan::adonis2(d ~ Group, data = meta, permutations = 999)
  expect_equal(ours$R2[1], ref$R2[1])
  expect_equal(ours$F[1], ref$F[1])
  expect_equal(ours$`Pr(>F)`[1], ref$`Pr(>F)`[1])
})

test_that("PERMDISP in mp_dispersion equals betadisper + permutest", {
  ps <- make_da_ps(n_per = 8)
  set.seed(7)
  df <- mp_dispersion(ps, group.by = "Group", return.data = TRUE)
  d <- phyloseq::distance(ps, method = "bray")
  g <- factor(data.frame(phyloseq::sample_data(ps))[labels(d), "Group"])
  set.seed(7)
  bd <- vegan::betadisper(d, g)
  expect_equal(df$distance, unname(bd$distances))
})

test_that("alpha diversity equals phyloseq::estimate_richness", {
  ps <- make_test_ps()
  ours <- mp_alpha(ps, measure = c("Shannon", "Observed"),
                   return.data = TRUE)
  ref <- suppressWarnings(
    phyloseq::estimate_richness(ps, measures = c("Shannon", "Observed")))
  expect_equal(ours$Shannon, ref$Shannon)
  expect_equal(ours$Observed, ref$Observed)
})

test_that("the lm engine equals per-taxon lm() on CLR values", {
  ps <- make_da_ps()
  res <- mp_run_da(ps, group.by = "Group", method = "lm", adjust = "Age",
                   prevalence = 0)
  counts <- as(phyloseq::otu_table(ps), "matrix")
  counts <- counts[rowSums(counts) > 0, ]
  logm <- log(t(counts) + 1)
  clr <- sweep(logm, 1, rowMeans(logm), "-")
  meta <- data.frame(phyloseq::sample_data(ps))[rownames(clr), ]
  for (taxon in c("OTU1", "OTU10", "OTU20")) {
    fit <- summary(stats::lm(clr[, taxon] ~ Group + Age, data = meta))
    co <- fit$coefficients["GroupB", ]
    row <- res[res$feature_id == taxon, ]
    expect_equal(row$log2FC, unname(co["Estimate"]) / log(2))
    expect_equal(row$lfcSE, unname(co["Std. Error"]) / log(2))
    expect_equal(row$pvalue, unname(co["Pr(>|t|)"]))
  }
})

test_that("the wilcox engine equals wilcox.test() on CLR values", {
  ps <- make_da_ps()
  res <- mp_run_da(ps, group.by = "Group", method = "wilcox",
                   prevalence = 0)
  counts <- as(phyloseq::otu_table(ps), "matrix")
  counts <- counts[rowSums(counts) > 0, ]
  logm <- log(t(counts) + 1)
  clr <- sweep(logm, 1, rowMeans(logm), "-")
  g <- data.frame(phyloseq::sample_data(ps))[rownames(clr), "Group"]
  for (taxon in c("OTU1", "OTU15")) {
    ref <- suppressWarnings(stats::wilcox.test(clr[g == "B", taxon],
                                               clr[g == "A", taxon]))
    row <- res[res$feature_id == taxon, ]
    expect_equal(row$pvalue, ref$p.value)
    expect_equal(row$log2FC, (mean(clr[g == "B", taxon]) -
                                mean(clr[g == "A", taxon])) / log(2))
  }
})

test_that("pairwise and omnibus p-values are exact, not rstatix-rounded", {
  set.seed(3)
  df <- data.frame(y = c(rnorm(10), rnorm(10, 0.8), rnorm(10, 0.3)),
                   g = factor(rep(c("a", "b", "c"), each = 10)))
  res <- micropubr:::.mp_pairwise_stats(df, "y", "g", test = "wilcoxon")
  ref <- mapply(function(a, b) stats::wilcox.test(df$y[df$g == a],
                                                  df$y[df$g == b])$p.value,
                as.character(res$group1), as.character(res$group2))
  expect_equal(res$p, unname(ref))
  expect_equal(res$p.adj, stats::p.adjust(unname(ref), "BH"))

  tk <- micropubr:::.mp_pairwise_stats(df, "y", "g", test = "anova")
  ref_tk <- stats::TukeyHSD(stats::aov(y ~ g, df))$g[, "p adj"]
  expect_equal(tk$p.adj, unname(ref_tk[paste(tk$group2, tk$group1,
                                             sep = "-")]))

  glob <- micropubr:::.mp_global_test(df, "y", "g", "kruskal")
  expect_equal(glob, sprintf("Kruskal-Wallis, p = %s",
                             format.pval(stats::kruskal.test(y ~ g, df)$p.value,
                                         digits = 2)))
})

test_that("PCoA axis percentages match ape::pcoa", {
  skip_if_not_installed("ape")
  ps <- make_da_ps()
  df <- mp_ordination(ps, method = "PCoA", distance = "bray",
                      permanova = FALSE, return.data = TRUE)
  ref <- ape::pcoa(phyloseq::distance(ps, method = "bray"))
  pct <- as.numeric(sub(".*\\((.*)%\\).*", "\\1", names(df)[1:2]))
  expect_equal(pct, round(100 * ref$values$Relative_eig[1:2], 1))
})

test_that("taxa unassigned at the plotted rank are not folded into Other", {
  ps <- make_test_ps(ntax = 8)
  tt <- as(phyloseq::tax_table(ps), "matrix")
  top <- names(sort(phyloseq::taxa_sums(ps), decreasing = TRUE))[1:2]
  tt[top, "Genus"] <- NA
  phyloseq::tax_table(ps) <- phyloseq::tax_table(tt)
  df <- mp_taxa(ps, tax.level = "Genus", top.n = 4, return.data = TRUE)
  labs <- levels(df$.taxon)
  expect_true(any(grepl("^Unassigned ", labs)))
  # Four top taxa shown individually, plus "Other".
  expect_equal(length(setdiff(labs, "Other")), 4)
})

test_that("placeholder rank values are not used as taxon names", {
  df <- data.frame(Kingdom = c("Root", "Bacteria", "Unassigned"),
                   Phylum = c(NA, NA, NA), Genus = c(NA, NA, NA))
  labs <- micropubr:::.mp_unassigned_labels(df, "Genus",
                                            c("Kingdom", "Phylum", "Genus"))
  expect_equal(labs, c("Unassigned", "Unassigned Bacteria", "Unassigned"))
})
