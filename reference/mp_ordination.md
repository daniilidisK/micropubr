# Ordination plot with ellipses and PERMANOVA

Computes an ordination from a
[phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
object and renders a publication-ready scatter plot. Distance-based
methods (PCoA, NMDS, dbRDA) accept any
[`phyloseq::distance()`](https://rdrr.io/pkg/phyloseq/man/distance.html)
metric; PCA operates on centred log-ratio (Aitchison) coordinates; RDA
and CCA support constrained ordination against metadata. Group ellipses
and a PERMANOVA result are added automatically.

## Usage

``` r
mp_ordination(
  object,
  method = c("PCoA", "NMDS", "PCA", "RDA", "CCA", "dbRDA"),
  distance = "bray",
  group.by = NULL,
  shape.by = NULL,
  ellipse = TRUE,
  ellipse.type = c("norm", "t", "euclid"),
  ellipse.level = 0.95,
  permanova = TRUE,
  permutations = 999,
  constraints = NULL,
  add.vectors = FALSE,
  envfit = NULL,
  axes = c(1, 2),
  label.samples = FALSE,
  point.size = 2.5,
  point.alpha = 0.85,
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

- method:

  One of `"PCoA"`, `"NMDS"`, `"PCA"`, `"RDA"`, `"CCA"`, `"dbRDA"`.

- distance:

  Distance metric for distance-based methods. Any value accepted by
  [`phyloseq::distance()`](https://rdrr.io/pkg/phyloseq/man/distance.html)
  (e.g. `"bray"`, `"jaccard"`, `"unifrac"`, `"wunifrac"`, `"canberra"`)
  plus `"aitchison"`.

- group.by:

  Metadata column mapped to point colour (and used for ellipses and
  PERMANOVA).

- shape.by:

  Optional metadata column mapped to point shape.

- ellipse:

  Draw group ellipses.

- ellipse.type:

  Ellipse type passed to
  [`ggplot2::stat_ellipse()`](https://ggplot2.tidyverse.org/reference/stat_ellipse.html).

- ellipse.level:

  Confidence level for ellipses.

- permanova:

  Run
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html)
  and append the result to the caption (distance-based methods only).

- permutations:

  Number of PERMANOVA permutations.

- constraints:

  For `"RDA"`, `"CCA"` and `"dbRDA"`, a character vector of metadata
  columns to constrain the ordination on.

- add.vectors:

  Overlay constraint/biplot arrows (constrained methods).

- envfit:

  Optional character vector of (numeric) metadata columns to fit onto
  the ordination as vectors via
  [`vegan::envfit()`](https://vegandevs.github.io/vegan/reference/envfit.html).
  Significant vectors are overlaid as arrows and summarised in the
  caption.

- axes:

  Length-2 integer vector of axes to plot.

- label.samples:

  Label individual sample points (via `ggrepel`).

- point.size:

  Point size.

- point.alpha:

  Point opacity.

- palette:

  Palette name, colour vector, or `NULL`. See
  [`mp_palette()`](https://daniilidisk.github.io/micropubr/reference/mp_palette.md).

- theme:

  A `ggplot2` theme; defaults to the active theme (see
  [`mp_set_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md),
  otherwise
  [`mp_theme_publication()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_publication.md)).

- ...:

  Additional arguments passed to the point geom
  ([`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)).

- plot.title, plot.subtitle, plot.caption:

  Optional labels; the caption is auto-generated when `NULL`.

- font.size:

  Base font size.

- interactive:

  If `TRUE`, return an interactive `plotly` widget instead of a static
  `ggplot` (requires `plotly`). See
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md).

- return.data:

  If `TRUE`, return the underlying data frame of sample scores and
  metadata (with the ordination object and any PERMANOVA/envfit results
  attached as attributes) instead of a plot.

- verbose:

  Print the methods used.

## Value

A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
data frame when `return.data = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(GlobalPatterns, package = "phyloseq")
mp_ordination(GlobalPatterns, method = "PCoA", distance = "bray",
              group.by = "SampleType")
} # }
```
