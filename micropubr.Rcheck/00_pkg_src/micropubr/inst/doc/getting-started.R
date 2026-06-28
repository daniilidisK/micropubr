## ----setup, include = FALSE---------------------------------------------------
has_deps <- requireNamespace("phyloseq", quietly = TRUE) &&
  requireNamespace("vegan", quietly = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  dpi = 120,
  eval = has_deps
)

## ----load---------------------------------------------------------------------
library(micropubr)
data(GlobalPatterns, package = "phyloseq")

# Work with a lightly filtered subset for speed.
gp <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(GlobalPatterns) > 5, GlobalPatterns)

## ----alpha--------------------------------------------------------------------
mp_alpha(gp, measure = c("Shannon", "Observed"),
         group.by = "SampleType", test = "kruskal")

## ----alpha-rarefaction, fig.height = 4----------------------------------------
mp_alpha(gp, plot.type = "rarefaction", group.by = "SampleType")

## ----ordination, fig.height = 5-----------------------------------------------
mp_ordination(gp, method = "PCoA", distance = "bray", group.by = "SampleType")

## ----ordination-pca, fig.height = 5-------------------------------------------
mp_ordination(gp, method = "PCA", group.by = "SampleType")

## ----taxa, fig.height = 5-----------------------------------------------------
mp_taxa(gp, tax.level = "Phylum", group.by = "SampleType", top.n = 8)

## ----diff, eval = TRUE, fig.height = 5----------------------------------------
set.seed(1)
res <- data.frame(
  feature = paste0("ASV", 1:200),
  log2FC  = rnorm(200, 0, 2),
  pvalue  = runif(200)^2
)
res$log2FC[1:20] <- rnorm(20, 4, 0.5) * sample(c(-1, 1), 20, TRUE)
res$pvalue[1:20] <- runif(20, 1e-9, 1e-4)
res$padj <- p.adjust(res$pvalue, "BH")

mp_diff(res, plot.type = "volcano", label.top = 8)

## ----combine, fig.width = 9, fig.height = 4-----------------------------------
fig <- mp_combine(
  mp_alpha(gp, measure = "Shannon", group.by = "SampleType"),
  mp_ordination(gp, group.by = "SampleType", permanova = FALSE),
  ncol = 2, tag.levels = "A"
)
fig

## ----export, eval = FALSE-----------------------------------------------------
# mp_export(fig, "figure1.pdf", column = "double", dpi = 600)

## ----theme--------------------------------------------------------------------
mp_taxa(gp, tax.level = "Phylum", top.n = 6,
        palette = "taxa20", theme = mp_theme_nature())

