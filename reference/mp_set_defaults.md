# Get, set and reset micropubr defaults

`micropubr` reads a small set of package-wide defaults that every
`mp_*()` plotting function falls back to when the corresponding argument
is left at its default (`NULL`). Setting them once with
`mp_set_defaults()` gives every subsequent plot a consistent look
without repeating arguments.

## Usage

``` r
mp_set_defaults(
  font.size = NULL,
  palette = NULL,
  theme = NULL,
  verbose = NULL,
  guardrails = NULL
)

mp_get_defaults()

mp_reset_defaults()
```

## Arguments

- font.size:

  Base font size in points used by the default theme.

- palette:

  A palette name (see
  [`mp_palette()`](https://daniilidisk.github.io/micropubr/reference/mp_palette.md))
  or an explicit vector of colours used whenever a function's `palette`
  argument is `NULL`.

- theme:

  A `ggplot2` theme object used whenever a function's `theme` argument
  is `NULL`. If `NULL`,
  [`mp_theme_publication()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_publication.md)
  is used.

- verbose:

  Logical; print the underlying method and parameters.

- guardrails:

  How statistical caveats are reported: `"warn"` (the default) raises a
  warning of class `micropubr_guardrail`, `"caption"` appends a "Note:"
  line to the figure caption, `"off"` stays silent. The caveats are
  recorded either way and listed by
  [`mp_methods()`](https://daniilidisk.github.io/micropubr/reference/mp_methods.md).
  They flag group sizes below 5, richness compared across more than
  10-fold differences in sequencing depth, a significant PERMANOVA
  accompanied by unequal dispersions (PERMDISP), differential-abundance
  methods that ignored the requested covariates, and covariates so
  strongly associated with the grouping that adjusting for them is
  unstable.

## Value

`mp_get_defaults()` returns the current defaults as a named list;
`mp_set_defaults()` and `mp_reset_defaults()` return the (updated)
defaults invisibly.

## Examples

``` r
old <- mp_get_defaults()
mp_set_defaults(font.size = 9, palette = "okabe")
mp_get_defaults()
#> $font.size
#> [1] 9
#> 
#> $palette
#> [1] "okabe"
#> 
#> $theme
#> NULL
#> 
#> $verbose
#> [1] FALSE
#> 
#> $guardrails
#> [1] "warn"
#> 
mp_reset_defaults()
```
