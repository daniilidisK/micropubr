# micropubr

> Publication-ready microbiome visualizations with one line of code.

`micropubr` is a unified, opinionated R package for turning a
[`phyloseq`](https://joey711.github.io/phyloseq/) object into a
journal-quality figure with a single call. It applies one consistent
design philosophy across the whole microbiome workflow — shared argument
names, colourblind-safe defaults, automatic statistical overlays and
reproducible, auto-generated captions — so that every analysis, from
alpha diversity to machine learning, is produced the same way.

Every `mp_*()` plotting function returns a plain
[`ggplot2`](https://ggplot2.tidyverse.org/) object, so you keep full
control: just add `+ theme(...)` or any other layer. Pass
`interactive = TRUE` to get a [`plotly`](https://plotly.com/r/) widget
instead — same call, hover tooltips and zoom included.

------------------------------------------------------------------------

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
  - [Differential abundance with several
    methods](#differential-abundance-with-several-methods)
  - [Methods text, alt text and statistical
    guardrails](#methods-text-alt-text-and-statistical-guardrails)
  - [Palettes and theming](#palettes-and-theming)
  - [Package-wide defaults](#package-wide-defaults)
  - [Importing data](#importing-data)
  - [Documentation](#documentation)
  - [Contributing](#contributing)
  - [License](#license)

## Why micropubr?

Most microbiome figures today are hand-assembled from `phyloseq`,
`vegan`, `ggplot2`, `ggtree` and a dozen helpers — dozens of lines per
panel, redone for every paper. `micropubr` collapses each of these into
one consistent call:

| Task | base phyloseq + ggplot2 | micropubr |
|----|----|----|
| Ordination with PERMANOVA + ellipses | ~25 lines | `mp_ordination(ps, group.by = "Group")` |
| Alpha diversity with significance brackets | ~20 lines | `mp_alpha(ps, group.by = "Group")` |
| Top-N stacked composition | ~15 lines | `mp_taxa(ps, tax.level = "Genus")` |

## Features

- **One consistent interface.** A single family of `mp_*()` functions
  with shared argument names (`group.by`, `palette`, `theme`,
  `plot.title`, `return.data`, …) spanning diversity, composition,
  differential abundance, phylogenetics, networks, longitudinal data,
  functional profiles, machine learning and genomics.
- **Publication defaults out of the box.** Colourblind-safe palettes, a
  clean journal theme, sensible typography and automatic, reproducible
  figure captions on every plot.
- **Statistics built in.** Kruskal–Wallis / pairwise tests, PERMANOVA,
  PERMDISP, ANOSIM, MRPP, Mantel tests, variation partitioning and
  indicator species analysis are computed and annotated automatically.
- **Static or interactive.** Every ggplot-based function accepts
  `interactive = TRUE` to return a linked `plotly` widget with the same
  call.
- **Composable and exportable.** Combine panels with
  [`mp_combine()`](https://daniilidisk.github.io/micropubr/reference/mp_combine.md)
  and write figures at journal column widths (or as self-contained HTML)
  with
  [`mp_export()`](https://daniilidisk.github.io/micropubr/reference/mp_export.md).
- **`ggplot2` all the way down.** Plots are ordinary `ggplot` objects —
  extend them with any layer, theme or scale you like.
- **`return.data = TRUE` everywhere.** Every function can hand back the
  tidy data frame (and, where relevant, the fitted model) behind the
  figure.
- **Differential abundance done carefully.**
  [`mp_run_da()`](https://daniilidisk.github.io/micropubr/reference/mp_run_da.md)
  runs Wilcoxon, covariate-adjusted linear (mixed) models, ANCOM-BC2,
  DESeq2 and ALDEx2 on one shared design, for two or many groups, and
  [`mp_diff()`](https://daniilidisk.github.io/micropubr/reference/mp_diff.md)
  shows where they agree.
- **Figures that write their own methods.** Every figure records what
  was computed;
  [`mp_methods()`](https://daniilidisk.github.io/micropubr/reference/mp_methods.md)
  turns that into a methods paragraph with versions and references, and
  [`mp_alt_text()`](https://daniilidisk.github.io/micropubr/reference/mp_alt_text.md)
  into accessible alt text.
- **Statistical guardrails.** Small groups, richness compared across
  very uneven depths, and a PERMANOVA undermined by unequal dispersions
  are flagged as they happen.
- **phyloseq or TreeSummarizedExperiment.** Every function accepts
  either container, so the `mia` / OMA ecosystem works out of the box.

## Installation

``` r

# Core dependencies (Bioconductor + CRAN)
install.packages("BiocManager")
BiocManager::install("phyloseq")

# micropubr (development version)
# install.packages("remotes")
remotes::install_github("daniilidisK/micropubr")
```

Optional features pull in additional packages on demand: `plotly` /
`htmlwidgets` (interactive plots), `ggtree` / `ggtreeExtra` (trees),
`ggraph` / `igraph` (networks), `ggalluvial` (longitudinal),
`treemapify` (treemaps), `circlize` (chord / Circos figures),
`randomForest` / `xgboost` / `pROC` (machine learning), `DESeq2` /
`ALDEx2` / `ANCOMBC` (differential abundance), `SpiecEasi` (network
inference), `Rsamtools` / `GenomicAlignments` (coverage), `ShortRead`
(FASTQ quality), `qiime2R` / `biomformat` (import) and `shiny`
(interactive explorer).

## Quick start

``` r

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
widget — hover for tooltips, zoom and pan, and toggle groups via the
legend. Multiple panels are linked into a single
[`plotly::subplot()`](https://rdrr.io/pkg/plotly/man/subplot.html), and
[`mp_export()`](https://daniilidisk.github.io/micropubr/reference/mp_export.md)
saves the result to a self-contained HTML page.

``` r

# An interactive ordination
mp_ordination(GlobalPatterns, group.by = "SampleType", interactive = TRUE)

# Save any interactive figure to a shareable web page
w <- mp_alpha(GlobalPatterns, group.by = "SampleType", interactive = TRUE)
mp_export(w, "alpha.html")
```

The circular `circlize`-based figures
([`mp_chord()`](https://daniilidisk.github.io/micropubr/reference/mp_chord.md),
[`mp_circos()`](https://daniilidisk.github.io/micropubr/reference/mp_circos.md),
[`mp_coverage()`](https://daniilidisk.github.io/micropubr/reference/mp_coverage.md))
are not `ggplot2` objects, so they have no interactive mode.

## Function reference

### Core visualizations

| Function | Purpose |
|----|----|
| [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md) | Alpha diversity (box / violin / dot / rarefaction) with tests |
| [`mp_ordination()`](https://daniilidisk.github.io/micropubr/reference/mp_ordination.md) | PCoA / NMDS / PCA / RDA / CCA / dbRDA with PERMANOVA |
| [`mp_taxa()`](https://daniilidisk.github.io/micropubr/reference/mp_taxa.md) | Composition: bar / area / heatmap / bubble / treemap |
| [`mp_diff()`](https://daniilidisk.github.io/micropubr/reference/mp_diff.md) | Differential abundance: volcano / MA / forest / dot / heatmap / consensus / structural zeros / concordance |
| [`mp_tree()`](https://daniilidisk.github.io/micropubr/reference/mp_tree.md) | Annotated phylogenetic trees (abundance + taxonomy rings) |
| [`mp_network()`](https://daniilidisk.github.io/micropubr/reference/mp_network.md) | Co-occurrence / correlation networks |
| [`mp_longitudinal()`](https://daniilidisk.github.io/micropubr/reference/mp_longitudinal.md) | Spaghetti / trajectory / alluvial / horizon over time |
| [`mp_function()`](https://daniilidisk.github.io/micropubr/reference/mp_function.md) | Functional pathway profiles (PICRUSt2 / HUMAnN) |

### Statistics and ecology

| Function | Purpose |
|----|----|
| [`mp_dispersion()`](https://daniilidisk.github.io/micropubr/reference/mp_dispersion.md) | Multivariate dispersion (PERMDISP) with ANOSIM / MRPP / PERMANOVA |
| [`mp_mantel()`](https://daniilidisk.github.io/micropubr/reference/mp_mantel.md) | Mantel test between community and explanatory distances |
| [`mp_varpart()`](https://daniilidisk.github.io/micropubr/reference/mp_varpart.md) | Variation partitioning across explanatory matrices |
| [`mp_prevalence()`](https://daniilidisk.github.io/micropubr/reference/mp_prevalence.md) | Prevalence–abundance and core-microbiome plots |
| [`mp_indicator()`](https://daniilidisk.github.io/micropubr/reference/mp_indicator.md) | Indicator species analysis |
| [`mp_ml()`](https://daniilidisk.github.io/micropubr/reference/mp_ml.md) | Machine learning (random forest / xgboost): importance, ROC, confusion |
| [`mp_run_da()`](https://daniilidisk.github.io/micropubr/reference/mp_run_da.md) | Run one or several differential-abundance methods (with covariates, repeated measures and multi-group contrasts) and return one tidy table |

### Genomics and quality control

| Function | Purpose |
|----|----|
| [`mp_coverage()`](https://daniilidisk.github.io/micropubr/reference/mp_coverage.md) | Genome / contig sequencing-depth (coverage) tracks from BAM/CRAM |
| [`mp_circos()`](https://daniilidisk.github.io/micropubr/reference/mp_circos.md) | Circular (Circos-style) genome plot with coverage / GC / link tracks |
| [`mp_chord()`](https://daniilidisk.github.io/micropubr/reference/mp_chord.md) | Chord diagram linking top taxa to a grouping variable |
| [`mp_qc_profile()`](https://daniilidisk.github.io/micropubr/reference/mp_qc_profile.md) | Per-cycle FASTQ read-quality profiles (DADA2 truncation) |

### Reporting

| Function | Purpose |
|----|----|
| [`mp_methods()`](https://daniilidisk.github.io/micropubr/reference/mp_methods.md) | Methods paragraph, software versions, references and statistical notes for a set of figures |
| [`mp_alt_text()`](https://daniilidisk.github.io/micropubr/reference/mp_alt_text.md) | Accessible alt text describing a figure and its key result |

### Infrastructure

| Function | Purpose |
|----|----|
| [`mp_combine()`](https://daniilidisk.github.io/micropubr/reference/mp_combine.md) | Compose multi-panel figures with shared legends and tags |
| [`mp_export()`](https://daniilidisk.github.io/micropubr/reference/mp_export.md) | Save static figures (PDF/PNG/…) or interactive widgets (HTML) |
| [`mp_palette()`](https://daniilidisk.github.io/micropubr/reference/mp_palette.md), [`mp_palettes()`](https://daniilidisk.github.io/micropubr/reference/mp_palettes.md) | Colourblind-safe qualitative / sequential / diverging palettes |
| [`mp_theme_publication()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_publication.md), [`mp_theme_nature()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_nature.md), [`mp_theme_cell()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_cell.md), [`mp_theme_dark()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_dark.md) | Journal-style themes |
| [`mp_darken()`](https://daniilidisk.github.io/micropubr/reference/mp_darken.md) | Recolour a finished figure for a dark background, keeping its layout |
| [`mp_set_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md), [`mp_get_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md), [`mp_reset_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md) | Package-wide defaults |
| [`mp_import_qiime2()`](https://daniilidisk.github.io/micropubr/reference/mp_import_qiime2.md), [`mp_import_dada2()`](https://daniilidisk.github.io/micropubr/reference/mp_import_dada2.md), [`mp_import_biom()`](https://daniilidisk.github.io/micropubr/reference/mp_import_biom.md), [`mp_import_metaphlan()`](https://daniilidisk.github.io/micropubr/reference/mp_import_metaphlan.md) | Import to `phyloseq` |
| [`mp_as_phyloseq()`](https://daniilidisk.github.io/micropubr/reference/mp_as_phyloseq.md) | Convert a `TreeSummarizedExperiment` (accepted directly by every function) |
| [`mp_explore()`](https://daniilidisk.github.io/micropubr/reference/mp_explore.md) | Interactive Shiny explorer that emits reproducible code |

## Differential abundance with several methods

Differential-abundance methods are known to disagree, and unadjusted
tests are easily confounded.
[`mp_run_da()`](https://daniilidisk.github.io/micropubr/reference/mp_run_da.md)
therefore runs one or several methods on a single shared design and
returns one long table (`feature`, `contrast`, `method`, `log2FC`,
`lfcSE`, `pvalue`, `padj`, `robust`, `structural_zero`):

| `method` | Model | Covariates (`adjust`) | Repeated measures (`random`) |
|----|----|:--:|:--:|
| `"wilcox"` | Wilcoxon rank-sum on CLR values (default, no extra packages) | – | – |
| `"lm"` | Linear model on CLR values; mixed model with `random` | ✓ | ✓ (`lmerTest`) |
| `"ancombc2"` | ANCOM-BC2: bias correction, structural zeros, pseudo-count sensitivity | ✓ | ✓ |
| `"deseq2"` | DESeq2 Wald test, `poscounts` size factors | ✓ | – |
| `"aldex2"` | ALDEx2 Welch’s t-test on Monte-Carlo CLR instances | – | – |

``` r

gp <- phyloseq::subset_samples(GlobalPatterns,
                               SampleType %in% c("Feces", "Skin", "Tongue"))

# All groups against a reference, three methods on one design
da <- mp_run_da(gp, group.by = "SampleType", reference = "Skin",
                method = c("wilcox", "lm", "ancombc2"), tax.level = "Genus")

# With your own metadata: adjust for confounders, model repeated measures
# mp_run_da(ps, group.by = "Treatment", method = c("lm", "ancombc2"),
#           adjust = c("Age", "Sex"), random = "Subject")

mp_diff(da, plot.type = "consensus")   # which methods agree, taxon by taxon
mp_diff(da, plot.type = "heatmap")     # fold changes across contrasts
mp_diff(da, plot.type = "volcano", contrast = "Feces vs Skin")
mp_diff(da, plot.type = "zeros")       # ANCOM-BC2 structural zeros
```

Fold changes are always log2 (ANCOM-BC2’s natural-log estimates are
converted), and ANCOM-BC2 results that fail its pseudo-count sensitivity
analysis are not counted as significant unless you pass
`robust.only = FALSE`.
[`mp_diff()`](https://daniilidisk.github.io/micropubr/reference/mp_diff.md)
also reads results produced elsewhere — `DESeq2`, `ALDEx2`,
`ANCOM-BC`/`ANCOM-BC2` objects, any data frame, or a named list of them
for a consensus or concordance plot:

``` r

mp_diff(list(DESeq2 = deseq_res, ALDEx2 = aldex_res), plot.type = "consensus")
```

## Methods text, alt text and statistical guardrails

Every figure records the analysis behind it: the tests and their
settings, transformations, filters, multiple-testing correction, group
sizes and package versions.
[`mp_methods()`](https://daniilidisk.github.io/micropubr/reference/mp_methods.md)
turns those records into a methods paragraph with references, ready to
edit into a manuscript;
[`mp_alt_text()`](https://daniilidisk.github.io/micropubr/reference/mp_alt_text.md)
writes alt text for the figure.

``` r

p1 <- mp_alpha(GlobalPatterns, group.by = "SampleType", test = "kruskal")
p2 <- mp_ordination(GlobalPatterns, group.by = "SampleType")
mp_methods(p1, p2)                   # or mp_methods() for every figure so far
mp_methods(p1, p2, format = "markdown")
mp_alt_text(p2)
```

While figures are made, micropubr also checks for conditions under which
their statistics are easy to over-read — groups smaller than five
samples, richness compared across a more than ten-fold range of
sequencing depth, a significant PERMANOVA alongside unequal dispersions
(PERMDISP), a method that had to ignore the covariates you asked for, a
covariate so entangled with the grouping that adjusting for it is
unstable — and warns. The same notes are listed by
[`mp_methods()`](https://daniilidisk.github.io/micropubr/reference/mp_methods.md);
`mp_set_defaults(guardrails = "caption")` writes them into the figure
caption instead, and `"off"` silences them.

## Palettes and theming

All palettes are colourblind-tolerant and designed for microbiome
figures — qualitative palettes for taxa and groups, plus sequential and
diverging ramps for heatmaps and continuous overlays.

``` r

mp_palettes()               # list available palettes
mp_palette("taxa10")        # a qualitative palette
mp_palette("diverging", n = 7)

# Swap the theme on any figure
mp_taxa(GlobalPatterns, tax.level = "Phylum", theme = mp_theme_nature())
```

## Package-wide defaults

Set a house style once and every `mp_*()` function will pick it up:

``` r

mp_set_defaults(font.size = 9, palette = "okabe")
mp_get_defaults()
mp_reset_defaults()
```

## Importing data

Bring data in from the common upstream pipelines and get a `phyloseq`
object back, ready for the rest of the family:

``` r

ps <- mp_import_qiime2(
  features = "table.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)
```

Importers are also provided for DADA2
([`mp_import_dada2()`](https://daniilidisk.github.io/micropubr/reference/mp_import_dada2.md)),
BIOM
([`mp_import_biom()`](https://daniilidisk.github.io/micropubr/reference/mp_import_biom.md))
and MetaPhlAn
([`mp_import_metaphlan()`](https://daniilidisk.github.io/micropubr/reference/mp_import_metaphlan.md)).

Data already in Bioconductor’s `TreeSummarizedExperiment` (for example
from `mia`) needs no conversion: every `mp_*()` function accepts it
directly, and
[`mp_as_phyloseq()`](https://daniilidisk.github.io/micropubr/reference/mp_as_phyloseq.md)
converts it explicitly when you want the `phyloseq` object.

## Documentation

Three vignettes walk through the package end to end:

- **Getting started with micropubr** — the core workflow and first
  figures.
- **Visualization gallery** — every plot type and its variants.
- **Statistics, machine learning and genomics** — tests,
  [`mp_ml()`](https://daniilidisk.github.io/micropubr/reference/mp_ml.md),
  coverage and Circos figures.

``` r

browseVignettes("micropubr")
```

Full per-function documentation is available with
[`?mp_alpha`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md)
(and friends).

## Contributing

Bug reports and feature requests are welcome via the [issue
tracker](https://github.com/daniilidisK/micropubr/issues). Pull requests
are appreciated — please keep new plotting functions within the shared
`mp_*()` argument conventions.

## License

MIT © Konstantinos Daniilidis — see
[LICENSE.md](https://daniilidisk.github.io/micropubr/LICENSE.md).
