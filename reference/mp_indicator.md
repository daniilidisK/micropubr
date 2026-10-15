# Indicator species analysis

Identifies taxa significantly associated with groups (or group
combinations) via
[`indicspecies::multipatt()`](https://emf-creaf.github.io/indicspecies/reference/multipatt.html)
(IndVal.g), and plots the strongest indicators coloured by their
associated group.

## Usage

``` r
mp_indicator(
  object,
  group.by = NULL,
  tax.level = NULL,
  top.n = 20,
  p.threshold = 0.05,
  permutations = 999,
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

- tax.level:

  Optional rank to agglomerate to first.

- top.n:

  Maximum number of indicator taxa to display.

- p.threshold:

  Significance cut-off on the permutation p-value.

- permutations:

  Number of permutations.

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

  If `TRUE`, return the indicator-statistics table.

- verbose:

  Print the methods used.

## Value

A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
data frame when `return.data = TRUE`.
