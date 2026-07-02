# micropubr

> Publication-ready microbiome visualizations with one line of code.

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
[![R >= 4.1.0](https://img.shields.io/badge/R-%3E%3D%204.1.0-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.2.0-green.svg)](NEWS.md)
<!-- badges: end -->

`micropubr` is a unified, opinionated R package for turning a
[`phyloseq`](https://joey711.github.io/phyloseq/) object into a
journal-quality figure with a single call. It applies one consistent design
philosophy across the whole microbiome workflow — shared argument names,
colourblind-safe defaults, automatic statistical overlays and reproducible,
auto-generated captions — so that every analysis, from alpha diversity to
machine learning, is produced the same way.

Every `mp_*()` plotting function returns a plain
[`ggplot2`](https://ggplot2.tidyverse.org/) object, so you keep full control:
just add `+ theme(...)` or any other layer. Pass `interactive = TRUE` to get a
[`plotly`](https://plotly.com/r/) widget instead — same call, hover tooltips
and zoom included.

---

## Table of contents

- [micropubr](#micropubr)
  - [Table of contents](#table-of-contents)
  - [Why micropubr?](#why-micropubr)
  - [Features](#features)
  - [Installation](#installation)
  - [Quick start](#quick-start)
  - [Interactive figures](#interactive-figures)
  - [Function reference](#function-reference)
    - [Core visualizations](#core-visualizations)
    - [Statistics and ecology](#statistics-and-ecology)
    - [Genomics and quality control](#genomics-and-quality-control)
    - [Infrastructure](#infrastructure)
  - [Differential abundance, any tool](#differential-abundance-any-tool)
  - [Palettes and theming](#palettes-and-theming)
  - [Package-wide defaults](#package-wide-defaults)
  - [Importing data](#importing-data)
  - [Documentation](#documentation)
  - [Contributing](#contributing)
  - [License](#license)

## Why micropubr?

Most microbiome figures today are hand-assembled from `phyloseq`, `vegan`,
`ggplot2`, `ggtree` and a dozen helpers — dozens of lines per panel, redone for
every paper. `micropubr` collapses each of these into one consistent call:

| Task | base phyloseq + ggplot2 | micropubr |
|---|---|---|
| Ordination with PERMANOVA + ellipses | ~25 lines | `mp_ordination(ps, group.by = "Group")` |
| Alpha diversity with significance brackets | ~20 lines | `mp_alpha(ps, group.by = "Group")` |
| Top-N stacked composition | ~15 lines | `mp_taxa(ps, tax.level = "Genus")` |

## Features

- **One consistent interface.** A single family of `mp_*()` functions with
  shared argument names (`group.by`, `palette`, `theme`, `plot.title`,
  `return.data`, …) spanning diversity, composition, differential abundance,
  phylogenetics, networks, longitudinal data, functional profiles, machine
  learning and genomics.
- **Publication defaults out of the box.** Colourblind-safe palettes, a clean
  journal theme, sensible typography and automatic, reproducible figure
  captions on every plot.
- **Statistics built in.** Kruskal–Wallis / pairwise tests, PERMANOVA,
  PERMDISP, ANOSIM, MRPP, Mantel tests, variation partitioning and indicator
  species analysis are computed and annotated automatically.
- **Static or interactive.** Every ggplot-based function accepts
  `interactive = TRUE` to return a linked `plotly` widget with the same call.
- **Composable and exportable.** Combine panels with `mp_combine()` and write
  figures at journal column widths (or as self-contained HTML) with
  `mp_export()`.
- **`ggplot2` all the way down.** Plots are ordinary `ggplot` objects — extend
  them with any layer, theme or scale you like.
- **`return.data = TRUE` everywhere.** Every function can hand back the tidy
  data frame (and, where relevant, the fitted model) behind the figure.

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
`igraph` (networks), `ggalluvial` (longitudinal), `treemapify` (treemaps),
`circlize` (chord / Circos figures), `randomForest` / `xgboost` / `pROC`
(machine learning), `DESeq2` / `ALDEx2` / `ANCOMBC` (differential abundance),
`SpiecEasi` (network inference), `Rsamtools` / `GenomicAlignments` (coverage),
`ShortRead` (FASTQ quality), `qiime2R` / `biomformat` (import) and `shiny`
(interactive explorer).

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

The circular `circlize`-based figures (`mp_chord()`, `mp_circos()`,
`mp_coverage()`) are not `ggplot2` objects, so they have no interactive mode.

## Function reference

### Core visualizations

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

### Statistics and ecology

| Function | Purpose |
|---|---|
| `mp_dispersion()` | Multivariate dispersion (PERMDISP) with ANOSIM / MRPP / PERMANOVA |
| `mp_mantel()` | Mantel test between community and explanatory distances |
| `mp_varpart()` | Variation partitioning across explanatory matrices |
| `mp_prevalence()` | Prevalence–abundance and core-microbiome plots |
| `mp_indicator()` | Indicator species analysis |
| `mp_ml()` | Machine learning (random forest / xgboost): importance, ROC, confusion |
| `mp_run_da()` | Run a differential-abundance method and return tidy results |

### Genomics and quality control

| Function | Purpose |
|---|---|
| `mp_coverage()` | Genome / contig sequencing-depth (coverage) tracks from BAM/CRAM |
| `mp_circos()` | Circular (Circos-style) genome plot with coverage / GC / link tracks |
| `mp_chord()` | Chord diagram linking top taxa to a grouping variable |
| `mp_qc_profile()` | Per-cycle FASTQ read-quality profiles (DADA2 truncation) |

### Infrastructure

| Function | Purpose |
|---|---|
| `mp_combine()` | Compose multi-panel figures with shared legends and tags |
| `mp_export()` | Save static figures (PDF/PNG/…) or interactive widgets (HTML) |
| `mp_palette()`, `mp_palettes()` | Colourblind-safe qualitative / sequential / diverging palettes |
| `mp_theme_publication()`, `mp_theme_nature()`, `mp_theme_cell()`, `mp_theme_dark()` | Journal-style themes |
| `mp_set_defaults()`, `mp_get_defaults()`, `mp_reset_defaults()` | Package-wide defaults |
| `mp_import_qiime2()`, `mp_import_dada2()`, `mp_import_biom()`, `mp_import_metaphlan()` | Import to `phyloseq` |
| `mp_explore()` | Interactive Shiny explorer that emits reproducible code |

## Differential abundance, any tool

`mp_diff()` accepts results from `DESeq2`, `ALDEx2` and `ANCOM-BC(2)` directly,
or any generic data frame. Compare two methods at a glance:

```r
mp_diff(deseq_res, result2 = aldex_res, plot.type = "concordance",
        method.names = c("DESeq2", "ALDEx2"))
```

Use `mp_run_da()` to run a method end-to-end from a `phyloseq` object and hand
the tidy result straight to `mp_diff()`.

## Palettes and theming

All palettes are colourblind-tolerant and designed for microbiome figures —
qualitative palettes for taxa and groups, plus sequential and diverging ramps
for heatmaps and continuous overlays.

```r
mp_palettes()               # list available palettes
mp_palette("taxa10")        # a qualitative palette
mp_palette("diverging", n = 7)

# Swap the theme on any figure
mp_taxa(GlobalPatterns, tax.level = "Phylum", theme = mp_theme_nature())
```

## Package-wide defaults

Set a house style once and every `mp_*()` function will pick it up:

```r
mp_set_defaults(font.size = 9, palette = "okabe")
mp_get_defaults()
mp_reset_defaults()
```

## Importing data

Bring data in from the common upstream pipelines and get a `phyloseq` object
back, ready for the rest of the family:

```r
ps <- mp_import_qiime2(
  feature.table = "table.qza",
  taxonomy      = "taxonomy.qza",
  metadata      = "metadata.tsv"
)
```

Importers are also provided for DADA2 (`mp_import_dada2()`), BIOM
(`mp_import_biom()`) and MetaPhlAn (`mp_import_metaphlan()`).

## Documentation

Three vignettes walk through the package end to end:

- **Getting started with micropubr** — the core workflow and first figures.
- **Visualization gallery** — every plot type and its variants.
- **Statistics, machine learning and genomics** — tests, `mp_ml()`, coverage
  and Circos figures.

```r
browseVignettes("micropubr")
```

Full per-function documentation is available with `?mp_alpha` (and friends).

<!-- ## Citation

If `micropubr` contributes to your work, please cite it:

```r
citation("micropubr")
``` -->

## Contributing

Bug reports and feature requests are welcome via the
[issue tracker](https://github.com/daniilidisK/micropubr/issues). Pull requests
are appreciated — please keep new plotting functions within the shared `mp_*()`
argument conventions.

## License

MIT © Konstantinos Daniilidis — see [LICENSE.md](LICENSE.md).
