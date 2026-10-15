# micropubr colour palettes

Colourblind-tolerant palettes designed for microbiome figures:
qualitative palettes for taxa and groups, plus sequential and diverging
ramps for heatmaps and fold-changes.

## Usage

``` r
mp_palette(name = "taxa10", n = NULL, reverse = FALSE)
```

## Arguments

- name:

  Palette name. One of `"taxa10"`, `"taxa20"`, `"taxa40"`, `"okabe"`,
  `"groups2"`..`"groups8"`, `"sequential"`, `"viridis"`, `"diverging"`.

- n:

  Number of colours to return. For qualitative palettes, if `n` exceeds
  the palette size a distinct fallback ramp is generated (with a
  message). For sequential/diverging palettes the colours are
  interpolated to exactly `n`. Defaults to the full palette.

- reverse:

  Reverse the colour order.

## Value

A character vector of hex colours.

## Examples

``` r
mp_palette("taxa10")
#>  [1] "#0072B2" "#E69F00" "#009E73" "#D55E00" "#CC79A7" "#56B4E9" "#F0E442"
#>  [8] "#332288" "#44AA99" "#117733"
mp_palette("diverging", n = 7)
#> [1] "#2166AC" "#4F92C3" "#96C3DC" "#F7F7F7" "#F1AE93" "#DA644F" "#B2182B"
```
