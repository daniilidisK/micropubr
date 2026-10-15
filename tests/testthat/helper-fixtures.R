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

# Counts with a planted differential-abundance signal: taxa 1-4 are 8-fold
# more abundant in group "B" (and "C"), taxa 5-6 8-fold less. Includes an Age
# covariate and a Subject id (two samples per subject) for adjusted and
# repeated-measures designs.
make_da_ps <- function(seed = 1, n_per = 8, ntax = 30, levels = c("A", "B")) {
  set.seed(seed)
  g <- rep(levels, each = n_per)
  n <- length(g)
  base <- stats::rlnorm(ntax, 4, 1)
  fold <- matrix(1, ntax, n)
  shifted <- g != levels[1]
  fold[1:4, shifted] <- 8
  fold[5:6, shifted] <- 1 / 8
  mu <- base * fold
  otu <- matrix(stats::rnbinom(ntax * n, mu = mu, size = 5), nrow = ntax,
                dimnames = list(paste0("OTU", seq_len(ntax)),
                                paste0("S", seq_len(n))))
  tax <- cbind(Kingdom = "Bacteria",
               Genus = paste0("Genus", seq_len(ntax)))
  rownames(tax) <- rownames(otu)
  meta <- data.frame(Group = g, Age = stats::runif(n, 20, 60),
                     Subject = paste0("P", ceiling(seq_len(n) / 2)),
                     row.names = colnames(otu))
  phyloseq::phyloseq(phyloseq::otu_table(otu, taxa_are_rows = TRUE),
                     phyloseq::tax_table(tax), phyloseq::sample_data(meta))
}

# A synthetic long DA table with two methods and two contrasts.
make_da_long <- function(seed = 1, n = 20) {
  set.seed(seed)
  grid <- expand.grid(feature = paste0("T", seq_len(n)),
                      contrast = c("B vs A", "C vs A"),
                      method = c("wilcox", "lm"), stringsAsFactors = FALSE)
  grid$feature_id <- grid$feature
  grid$log2FC <- stats::rnorm(nrow(grid), 0, 0.5)
  grid$pvalue <- stats::runif(nrow(grid), 0.2, 1)
  sig <- grid$feature %in% c("T1", "T2", "T3")
  grid$log2FC[sig] <- 3
  grid$pvalue[sig] <- 1e-6
  grid$padj <- grid$pvalue
  grid$structural_zero <- FALSE
  grid
}
