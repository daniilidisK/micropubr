# Multivariate dispersion (PERMDISP) with optional group tests

Runs
[`vegan::betadisper()`](https://vegandevs.github.io/vegan/reference/betadisper.html)
to test homogeneity of multivariate dispersions (PERMDISP) – the
standard companion to PERMANOVA – and plots each sample's distance to
its group centroid. Additional group-difference statistics (ANOSIM,
MRPP, PERMANOVA) can be appended to the caption.

## Usage

``` r
mp_dispersion(
  object,
  group.by = NULL,
  distance = "bray",
  add.test = NULL,
  permutations = 999,
  plot.type = c("box", "violin"),
  point.size = 1.5,
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

- group.by:

  Sample-metadata column defining the groups.

- distance:

  Distance metric (see
  [`mp_ordination()`](https://daniilidisk.github.io/micropubr/reference/mp_ordination.md);
  default `"bray"`).

- add.test:

  Optional character vector of further tests to summarise in the
  caption: any of `"anosim"`, `"mrpp"`, `"permanova"`.

- permutations:

  Number of permutations.

- plot.type:

  One of `"box"` or `"violin"`.

- point.size:

  Size of jittered sample points.

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

  If `TRUE`, return the per-sample distance-to-centroid data frame (with
  the betadisper object attached) instead of a plot.

- verbose:

  Print the methods used.

## Value

A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
data frame when `return.data = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(GlobalPatterns, package = "phyloseq")
mp_dispersion(GlobalPatterns, group.by = "SampleType",
              add.test = "permanova")
} # }
```
