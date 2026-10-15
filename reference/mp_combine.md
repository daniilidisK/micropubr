# Combine micropubr plots into a labelled multi-panel figure

A thin convenience wrapper around
[`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)
that applies a consistent layout, collects duplicate legends and adds
automatic panel tags (A, B, C ...), so figures composed from several
`mp_*()` calls look coherent.

## Usage

``` r
mp_combine(
  ...,
  ncol = NULL,
  nrow = NULL,
  guides = c("collect", "keep", "auto"),
  tag.levels = "A",
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  widths = NULL,
  heights = NULL,
  interactive = FALSE
)
```

## Arguments

- ...:

  `ggplot` objects, or a single named list of them.

- ncol, nrow:

  Grid dimensions (passed to
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)).

- guides:

  How to treat legends: `"collect"` (default), `"keep"` or `"auto"`.

- tag.levels:

  Panel-tag style passed to
  [`patchwork::plot_annotation()`](https://patchwork.data-imaginist.com/reference/plot_annotation.html),
  e.g. `"A"`, `"a"`, `"1"`. Use `NULL` for no tags.

- title, subtitle, caption:

  Optional figure-level annotations.

- widths, heights:

  Optional relative panel sizes.

- interactive:

  If `TRUE`, return an interactive `plotly` subplot instead of a static
  `patchwork` (requires `plotly`). Panel tags and figure-level
  annotations are not carried over by the conversion.

## Value

A `patchwork` object, or a `plotly` widget when `interactive = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
a <- mp_alpha(ps, group.by = "Group")
o <- mp_ordination(ps, group.by = "Group")
mp_combine(a, o, ncol = 2, tag.levels = "A")
} # }
```
