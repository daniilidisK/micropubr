# Plot functional pathway profiles

Visualises functional-prediction tables such as PICRUSt2 pathway
abundances or HUMAnN pathway/gene-family tables. Accepts either a
[phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
object (whose "taxa" are pathways) or a plain features-by-samples matrix
/ data frame plus optional sample metadata.

## Usage

``` r
mp_function(
  object,
  metadata = NULL,
  plot.type = c("heatmap", "bar", "bubble"),
  group.by = NULL,
  top.n = 25,
  category = NULL,
  relative = TRUE,
  scale = c("none", "row"),
  palette = NULL,
  theme = NULL,
  plot.title = NULL,
  plot.subtitle = NULL,
  plot.caption = NULL,
  font.size = getOption("micropubr.font.size", 11),
  interactive = FALSE,
  verbose = getOption("micropubr.verbose", FALSE)
)
```

## Arguments

- object:

  A [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
  object, or a matrix / data frame with pathways in rows and samples in
  columns.

- metadata:

  Optional sample metadata data frame (rows = samples) used when
  `object` is a matrix.

- plot.type:

  One of `"heatmap"`, `"bar"` or `"bubble"`.

- group.by:

  Sample metadata column used to facet / summarise.

- top.n:

  Number of most-abundant pathways to display.

- category:

  Optional column (when `object` carries feature metadata, or a named
  vector) grouping pathways into higher-level categories for the bar
  plot (e.g. KEGG level 1/2).

- relative:

  Convert to per-sample relative abundance.

- scale:

  Row scaling for the heatmap: `"none"` or `"row"` (z-score).

- palette:

  Palette name or colour vector.

- theme:

  A `ggplot2` theme; defaults to
  [`mp_theme_publication()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_publication.md).

- plot.title, plot.subtitle, plot.caption:

  Optional labels.

- font.size:

  Base font size.

- interactive:

  If `TRUE`, return an interactive `plotly` widget instead of a static
  `ggplot` (requires `plotly`). See
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md).

- verbose:

  Print the steps used.

## Value

A `ggplot` object, or a `plotly` widget when `interactive = TRUE`.

## Examples

``` r
set.seed(1)
m <- matrix(rpois(20 * 12, 50), nrow = 20,
            dimnames = list(paste0("path", 1:20), paste0("S", 1:12)))
meta <- data.frame(Group = rep(c("A", "B"), each = 6),
                   row.names = colnames(m))
mp_function(m, metadata = meta, plot.type = "heatmap", group.by = "Group")
```
