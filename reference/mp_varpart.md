# Variance partitioning of community structure

Partitions the variation in community composition
(Hellinger-transformed) among two or three sets of explanatory metadata
variables using
[`vegan::varpart()`](https://vegandevs.github.io/vegan/reference/varpart.html),
and plots the adjusted-R-squared of each fraction.

## Usage

``` r
mp_varpart(
  object,
  groups,
  transform = "hellinger",
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

- groups:

  A **named** list mapping each explanatory set to a character vector of
  metadata columns, e.g.
  `list(Host = c("Age", "Sex"), Diet = "Treatment")`. Length 2 or 3.

- transform:

  Community standardisation passed to
  [`vegan::decostand()`](https://vegandevs.github.io/vegan/reference/decostand.html)
  (default `"hellinger"`).

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

  If `TRUE`, return the fractions table instead of a plot.

- verbose:

  Print the methods used.

## Value

A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
data frame when `return.data = TRUE`.
