# Convert a finished figure to dark mode

Recolours a finished figure for a dark background while keeping its
layout. Adding
[`mp_theme_dark()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_dark.md)
to a finished plot replaces its whole theme, which undoes layout that
the plotting function set (rotated or wrapped axis labels, bold facet
strips, the absence of grid lines), and it cannot reach dark text and
lines that the layers draw themselves: feature labels, significance
brackets, tree branches. `mp_darken()` changes only colours: the
backgrounds and the theme's text and line colours, plus any near-black
colour a layer sets that is not mapped to data.

## Usage

``` r
mp_darken(plot, background = "#1B1B1B", foreground = "grey90")
```

## Arguments

- plot:

  A `ggplot` or `patchwork` figure.

- background:

  Background colour.

- foreground:

  Colour for text, axes and formerly dark marks.

## Value

The recoloured figure; `plot` itself is not modified.

## See also

[`mp_theme_dark()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_dark.md)
to build a dark figure from the start.

## Examples

``` r
df <- data.frame(feature = paste0("f", 1:30), log2FC = rnorm(30, sd = 2),
                 padj = runif(30)^3)
mp_darken(mp_diff(df, plot.type = "volcano"))
```
