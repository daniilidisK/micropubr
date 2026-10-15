# Export a figure at publication resolution and journal column widths

Wraps
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
with sensible defaults for journal submission: device inferred from the
file extension, millimetre units, 300 dpi, and optional standard
column-width presets so figures are sized correctly for single- or
double-column layouts.

## Usage

``` r
mp_export(
  plot,
  filename,
  width = NULL,
  height = NULL,
  units = c("mm", "cm", "in"),
  dpi = 300,
  column = c("none", "single", "onehalf", "double"),
  format = NULL,
  ...
)
```

## Arguments

- plot:

  A `ggplot` or `patchwork` object, or an interactive `plotly` widget
  (as returned with `interactive = TRUE`), which is saved as a
  self-contained `.html` page instead.

- filename:

  Output path. For static plots the extension selects the device
  (`.pdf`, `.png`, `.tiff`, `.svg`, `.eps`, `.jpg`); for an interactive
  widget it must be `.html`.

- width, height:

  Dimensions in `units`. If `column` is set, `width` is taken from the
  preset; if `height` is `NULL` it defaults to `width * 0.75`.

- units:

  One of `"mm"` (default), `"cm"` or `"in"`.

- dpi:

  Raster resolution. Common choices: 300 or 600.

- column:

  Column-width preset: `"none"` (default), `"single"` (89 mm),
  `"onehalf"` (120 mm) or `"double"` (183 mm).

- format:

  Optional explicit device name overriding the extension.

- ...:

  Further arguments passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Value

The `filename`, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- mp_ordination(ps, group.by = "Group")
mp_export(p, "figure1.pdf", column = "single")
mp_export(p, "figure1.tiff", column = "double", dpi = 600)
} # }
```
