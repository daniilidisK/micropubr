# Plot taxonomic composition

Visualises community composition at a chosen taxonomic rank as stacked
bars, smoothed area charts, heatmaps or bubble plots. The top `top.n`
taxa are shown individually and the remainder are collapsed into an
"Other" category with a muted colour.

## Usage

``` r
mp_taxa(
  object,
  tax.level = "Phylum",
  plot.type = c("bar", "area", "heatmap", "bubble", "treemap"),
  group.by = NULL,
  top.n = 10,
  relative = TRUE,
  order.by = "abundance",
  cluster.samples = FALSE,
  other.label = "Other",
  palette = NULL,
  theme = NULL,
  ...,
  plot.title = NULL,
  plot.subtitle = NULL,
  plot.caption = NULL,
  font.size = getOption("micropubr.font.size", 11),
  interactive = FALSE,
  return.data = FALSE,
  verbose = getOption("micropubr.verbose", FALSE)
)
```

## Arguments

- object:

  A [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
  object.

- tax.level:

  Taxonomic rank to agglomerate to (e.g. `"Phylum"`, `"Genus"`). Must be
  one of
  [`phyloseq::rank_names()`](https://rdrr.io/pkg/phyloseq/man/rank_names.html).

- plot.type:

  One of `"bar"`, `"area"`, `"heatmap"`, `"bubble"` or `"treemap"` (the
  latter requires the `treemapify` package).

- group.by:

  Optional metadata column to facet samples by.

- top.n:

  Number of most-abundant taxa to display individually.

- relative:

  Show relative abundance (`TRUE`) or raw counts.

- order.by:

  How to order samples on the x-axis: a metadata column name,
  `"abundance"` (by the most dominant taxon, the default), or
  `"cluster"` (hierarchical clustering on Bray-Curtis distances).
  Ignored when `cluster.samples = TRUE`.

- cluster.samples:

  Order samples by hierarchical clustering; shorthand for
  `order.by = "cluster"`.

- other.label:

  Label for the collapsed low-abundance category.

- palette:

  Palette name or colour vector; defaults to a taxon palette.

- theme:

  A `ggplot2` theme; defaults to the active theme (see
  [`mp_set_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md),
  otherwise
  [`mp_theme_publication()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_publication.md)).

- ...:

  Additional arguments passed to the primary geom (e.g.
  [`ggplot2::geom_col()`](https://ggplot2.tidyverse.org/reference/geom_bar.html)
  for bars).

- plot.title, plot.subtitle, plot.caption:

  Optional labels.

- font.size:

  Base font size.

- interactive:

  If `TRUE`, return an interactive `plotly` widget instead of a static
  `ggplot` (requires `plotly`). See
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md).

- return.data:

  If `TRUE`, return the aggregated taxon-by-sample tidy data frame
  instead of a plot.

- verbose:

  Print the methods used.

## Value

A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
data frame when `return.data = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(GlobalPatterns, package = "phyloseq")
mp_taxa(GlobalPatterns, tax.level = "Phylum", group.by = "SampleType")
} # }
```
