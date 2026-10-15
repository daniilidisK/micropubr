# Mantel test between community and explanatory distances

Correlates the community dissimilarity matrix with a second distance
matrix – typically an environmental/spatial distance derived from
numeric metadata – via
[`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html),
and visualises the relationship as a distance-distance scatter or a
Mantel correlogram.

## Usage

``` r
mp_mantel(
  object,
  distance = "bray",
  vars = NULL,
  dist2 = NULL,
  method = c("pearson", "spearman"),
  permutations = 999,
  plot.type = c("scatter", "correlogram"),
  point.size = 1.2,
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

- distance:

  Community distance metric (default `"bray"`).

- vars:

  Numeric metadata column(s) used to build the second (Euclidean,
  standardised) distance matrix. Ignored if `dist2` is given.

- dist2:

  An explicit second distance (a `dist` object or square matrix).

- method:

  Correlation method: `"pearson"` or `"spearman"`.

- permutations:

  Number of permutations.

- plot.type:

  `"scatter"` or `"correlogram"`.

- point.size:

  Point size for the scatter plot.

- palette, theme, plot.title, plot.subtitle, plot.caption, font.size:

  See
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md).

- ...:

  Passed to the primary geom.

- interactive:

  If `TRUE`, return an interactive `plotly` widget instead of a static
  `ggplot` (requires `plotly`). See
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md).

- return.data:

  If `TRUE`, return the underlying data frame instead of a plot
  (pairwise distances for `"scatter"`, the correlogram table for
  `"correlogram"`).

- verbose:

  Print the methods used.

## Value

A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
data frame when `return.data = TRUE`.
