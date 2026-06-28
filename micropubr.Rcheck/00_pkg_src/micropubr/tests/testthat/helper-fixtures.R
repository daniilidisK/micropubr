# Shared test fixtures -----------------------------------------------------

# Build a small, deterministic phyloseq object for tests. Optionally attach a
# random phylogenetic tree.
make_test_ps <- function(seed = 1, tree = FALSE, ntax = 8, nsmp = 12) {
  set.seed(seed)
  otu <- matrix(stats::rpois(ntax * nsmp, lambda = 20), nrow = ntax,
                dimnames = list(paste0("OTU", seq_len(ntax)),
                                paste0("Sample", seq_len(nsmp))))
  tax <- cbind(
    Kingdom = "Bacteria",
    Phylum  = sample(c("Firmicutes", "Bacteroidetes", "Proteobacteria"),
                     ntax, replace = TRUE),
    Genus   = paste0("Genus", seq_len(ntax))
  )
  rownames(tax) <- rownames(otu)
  meta <- data.frame(
    Group = rep(c("Control", "Treatment"), each = nsmp / 2),
    Sex   = rep(c("M", "F"), length.out = nsmp),
    pH    = stats::runif(nsmp, 5, 8),
    Age   = stats::runif(nsmp, 20, 60),
    row.names = colnames(otu)
  )
  comps <- list(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax),
    phyloseq::sample_data(meta)
  )
  if (tree) {
    comps <- c(comps, list(phyloseq::phy_tree(
      ape::rtree(ntax, tip.label = rownames(otu)))))
  }
  do.call(phyloseq::phyloseq, comps)
}

# A generic differential-abundance table with planted signal.
make_da_df <- function(seed = 1, n = 60) {
  set.seed(seed)
  df <- data.frame(
    feature = paste0("ASV", seq_len(n)),
    log2FC  = stats::rnorm(n, 0, 2),
    pvalue  = stats::runif(n)^2,
    baseMean = stats::rlnorm(n, 5),
    lfcSE   = abs(stats::rnorm(n, 0.5, 0.2)),
    stringsAsFactors = FALSE
  )
  df$log2FC[1:10] <- stats::rnorm(10, 4, 0.4) * sample(c(-1, 1), 10, TRUE)
  df$pvalue[1:10] <- stats::runif(10, 1e-9, 1e-5)
  df$padj <- stats::p.adjust(df$pvalue, "BH")
  df
}
