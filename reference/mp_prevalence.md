# Prevalence-abundance and core-microbiome plots

Summarises how widely each taxon is detected. `plot.type = "scatter"`
plots prevalence against mean relative abundance (a standard QC /
core-selection view); `plot.type = "core"` draws a core-microbiome
heatmap of prevalence across a gradient of detection thresholds.

## Usage

``` r
mp_prevalence(
  object,
  tax.level = NULL,
  plot.type = c("scatter", "core"),
  detection = 0.001,
  prevalence = 0.5,
  color.by = NULL,
  label.top = 10,
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
  object (counts).

- tax.level:

  Optional rank to agglomerate to first.

- plot.type:

  `"scatter"` or `"core"`.

- detection:

  Detection threshold (relative abundance) defining presence for the
  scatter plot's prevalence and the core threshold line.

- prevalence:

  Prevalence threshold (fraction of samples) for the core line / core
  membership.

- color.by:

  For `"scatter"`, a taxonomy rank to colour points by.

- label.top:

  Number of most-prevalent taxa to label (scatter).

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

  If `TRUE`, return the prevalence/abundance table.

- verbose:

  Print the methods used.

## Value

A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
data frame when `return.data = TRUE`.
