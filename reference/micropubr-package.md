# micropubr: Publication-Ready Microbiome Visualizations

micropubr provides a unified, opinionated family of `mp_*()` functions
that turn a
[phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
object into a publication-ready
[ggplot2](https://ggplot2.tidyverse.org/reference/ggplot.html) figure
with a single call. It applies one consistent design philosophy across
the whole microbiome workflow: consistent arguments, colourblind-safe
defaults, automatic statistical overlays and reproducible,
auto-generated figure captions, while always returning a plain `ggplot`
object you can extend with `+`.

## Core plotting functions

- [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md):

  Alpha (within-sample) diversity with statistics.

- [`mp_ordination()`](https://daniilidisk.github.io/micropubr/reference/mp_ordination.md):

  Beta diversity / ordination with PERMANOVA.

- [`mp_taxa()`](https://daniilidisk.github.io/micropubr/reference/mp_taxa.md):

  Taxonomic composition (bar, area, heatmap, bubble).

- [`mp_diff()`](https://daniilidisk.github.io/micropubr/reference/mp_diff.md):

  Differential-abundance result visualisation.

- [`mp_tree()`](https://daniilidisk.github.io/micropubr/reference/mp_tree.md):

  Annotated phylogenetic trees.

- [`mp_network()`](https://daniilidisk.github.io/micropubr/reference/mp_network.md):

  Co-occurrence / correlation networks.

- [`mp_longitudinal()`](https://daniilidisk.github.io/micropubr/reference/mp_longitudinal.md):

  Time-course and paired-sample plots.

- [`mp_function()`](https://daniilidisk.github.io/micropubr/reference/mp_function.md):

  Functional pathway profiles (PICRUSt2/HUMAnN).

## Statistical analyses

- [`mp_dispersion()`](https://daniilidisk.github.io/micropubr/reference/mp_dispersion.md):

  Multivariate dispersion (PERMDISP) with optional ANOSIM / MRPP /
  PERMANOVA.

- [`mp_mantel()`](https://daniilidisk.github.io/micropubr/reference/mp_mantel.md):

  Mantel test and correlogram.

- [`mp_varpart()`](https://daniilidisk.github.io/micropubr/reference/mp_varpart.md):

  Variance partitioning among explanatory sets.

- [`mp_prevalence()`](https://daniilidisk.github.io/micropubr/reference/mp_prevalence.md):

  Prevalence-abundance and core-microbiome plots.

- [`mp_indicator()`](https://daniilidisk.github.io/micropubr/reference/mp_indicator.md):

  Indicator species analysis.

- [`mp_run_da()`](https://daniilidisk.github.io/micropubr/reference/mp_run_da.md):

  Run a differential-abundance method for
  [`mp_diff()`](https://daniilidisk.github.io/micropubr/reference/mp_diff.md).

## Microbial genomics

- [`mp_coverage()`](https://daniilidisk.github.io/micropubr/reference/mp_coverage.md):

  Genome / contig sequencing-depth tracks.

- [`mp_circos()`](https://daniilidisk.github.io/micropubr/reference/mp_circos.md):

  Circular (Circos-style) genome plots.

## Infrastructure

[`mp_palette()`](https://daniilidisk.github.io/micropubr/reference/mp_palette.md),
the `mp_theme_*()` family,
[`mp_set_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md),
[`mp_combine()`](https://daniilidisk.github.io/micropubr/reference/mp_combine.md),
[`mp_export()`](https://daniilidisk.github.io/micropubr/reference/mp_export.md)
and the `mp_import_*()` helpers. Every plotting function accepts
`return.data = TRUE` to return the tidy data frame behind the plot, and
`...` to forward arguments to its primary geom.

## See also

Useful links:

- <https://github.com/daniilidisK/micropubr>

- <https://daniilidisk.github.io/micropubr/>

- Report bugs at <https://github.com/daniilidisK/micropubr/issues>

## Author

**Maintainer**: Konstantinos Daniilidis <daniilidis@ieee.org>
