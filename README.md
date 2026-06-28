# micropubr <img src="man/figures/logo.png" align="right" height="120" alt="" />

> Publication-ready microbiome visualizations with one line of code.

`micropubr` is a unified, opinionated R package for turning a
[`phyloseq`](https://joey711.github.io/phyloseq/) object into a
journal-quality figure with a single call. It brings the design philosophy of
[`SCpubr`](https://enblacar.github.io/SCpubr-book/) — consistent arguments,
colourblind-safe defaults, automatic statistics, reproducible captions — to
microbiome science.

Every `mp_*()` function returns a plain
[`ggplot2`](https://ggplot2.tidyverse.org/) object, so you keep full control:
just add `+ theme(...)` or any other layer. Pass `interactive = TRUE` to get a
[`plotly`](https://plotly.com/r/) widget instead — same call, hover tooltips
and zoom included.

## Why micropubr?

Most microbiome figures today are hand-assembled from `phyloseq`, `vegan`,
`ggplot2`, `ggtree` and a dozen helpers — dozens of lines per panel, redone for
every paper. `micropubr` collapses each of these into one consistent call:

| | base phyloseq + ggplot2 | micropubr |
|---|---|---|
| Ordination with PERMANOVA + ellipses | ~25 lines | `mp_ordination(ps, group.by = "Group")` |
| Alpha diversity with significance brackets | ~20 lines | `mp_alpha(ps, group.by = "Group")` |
| Top-N stacked composition | ~15 lines | `mp_taxa(ps, tax.level = "Genus")` |

## Installation

```r
# Core dependencies (Bioconductor + CRAN)
install.packages("BiocManager")
BiocManager::install("phyloseq")

# micropubr (development version)
# install.packages("remotes")
remotes::install_github("daniilidisK/micropubr")
```

Optional features pull in additional packages on demand: `plotly` /
`htmlwidgets` (interactive plots), `ggtree` / `ggtreeExtra` (trees), `ggraph` /
`igraph` (networks), `ggalluvial` (longitudinal), `SpiecEasi` (network
inference), `qiime2R` (import) and `shiny` (interactive explorer).

## Quick start

```r
library(micropubr)
data(GlobalPatterns, package = "phyloseq")

# Alpha diversity with automatic Kruskal–Wallis + pairwise brackets
mp_alpha(GlobalPatterns, measure = c("Shannon", "Observed"),
         group.by = "SampleType", test = "kruskal")

# Beta diversity: PCoA on Bray–Curtis with ellipses + PERMANOVA caption
mp_ordination(GlobalPatterns, method = "PCoA", distance = "bray",
              group.by = "SampleType")

# Taxonomic composition: top 10 phyla, faceted by sample type
mp_taxa(GlobalPatterns, tax.level = "Phylum", group.by = "SampleType")

# Compose a multi-panel figure and export at journal column width
fig <- mp_combine(
  mp_alpha(GlobalPatterns, group.by = "SampleType"),
  mp_ordination(GlobalPatterns, group.by = "SampleType"),
  ncol = 2, tag.levels = "A"
)
mp_export(fig, "figure1.pdf", column = "double", dpi = 600)
```

## Interactive figures

Add `interactive = TRUE` to any plotting function to get a live `plotly`
widget — hover for tooltips, zoom and pan, and toggle groups via the legend.
Multiple panels are linked into a single `plotly::subplot()`, and `mp_export()`
saves the result to a self-contained HTML page.

```r
# An interactive ordination
mp_ordination(GlobalPatterns, group.by = "SampleType", interactive = TRUE)

# Save any interactive figure to a shareable web page
w <- mp_alpha(GlobalPatterns, group.by = "SampleType", interactive = TRUE)
mp_export(w, "alpha.html")
```

The circular `circlize`-based figures (`mp_chord()`, `mp_circos()`) are not
`ggplot2` objects, so they have no interactive mode.

## The `mp_*()` family

| Function | Purpose |
|---|---|
| `mp_alpha()` | Alpha diversity (box / violin / dot / rarefaction) with tests |
| `mp_ordination()` | PCoA / NMDS / PCA / RDA / CCA / dbRDA with PERMANOVA |
| `mp_taxa()` | Composition: bar / area / heatmap / bubble / treemap |
| `mp_diff()` | Differential abundance: volcano / MA / forest / dot / concordance |
| `mp_tree()` | Annotated phylogenetic trees (abundance + taxonomy rings) |
| `mp_network()` | Co-occurrence / correlation networks |
| `mp_longitudinal()` | Spaghetti / trajectory / alluvial / horizon over time |
| `mp_function()` | Functional pathway profiles (PICRUSt2 / HUMAnN) |
| `mp_explore()` | Interactive Shiny explorer with reproducible code |

Plus infrastructure: `mp_palette()`, `mp_theme_*()`, `mp_combine()`,
`mp_export()` and the `mp_import_*()` helpers.

## Differential abundance, any tool

`mp_diff()` accepts results from `DESeq2`, `ALDEx2` and `ANCOM-BC(2)` directly,
or any generic data frame. Compare two methods at a glance:

```r
mp_diff(deseq_res, result2 = aldex_res, plot.type = "concordance",
        method.names = c("DESeq2", "ALDEx2"))
```

## License

MIT © Konstantinos Daniilidis
