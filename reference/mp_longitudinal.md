# Plot longitudinal and paired-sample microbiome data

Visualises time-course or before/after designs: per-subject diversity
trajectories (spaghetti), ordination trajectories connecting paired
timepoints, alluvial flows of taxonomic shifts, and small-multiple
"horizon" area bands of taxon dynamics.

## Usage

``` r
mp_longitudinal(
  object,
  time.var,
  subject.var,
  plot.type = c("spaghetti", "trajectory", "alluvial", "horizon"),
  group.by = NULL,
  measure = "Shannon",
  tax.level = "Phylum",
  top.n = 8,
  distance = "bray",
  smooth = TRUE,
  palette = NULL,
  theme = NULL,
  plot.title = NULL,
  plot.subtitle = NULL,
  plot.caption = NULL,
  font.size = getOption("micropubr.font.size", 11),
  interactive = FALSE,
  verbose = getOption("micropubr.verbose", FALSE)
)
```

## Arguments

- object:

  A [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
  object.

- time.var:

  Metadata column giving the timepoint (numeric or ordered).

- subject.var:

  Metadata column identifying the subject / individual.

- plot.type:

  One of `"spaghetti"`, `"trajectory"`, `"alluvial"` or `"horizon"`.

- group.by:

  Optional metadata column for colouring (e.g. treatment arm).

- measure:

  Alpha-diversity index for `"spaghetti"` plots.

- tax.level:

  Taxonomic rank for `"alluvial"` / `"horizon"` plots.

- top.n:

  Number of taxa to display in `"alluvial"` / `"horizon"`.

- distance:

  Distance metric for `"trajectory"` ordination.

- smooth:

  Add a per-group LOESS smooth to `"spaghetti"` plots.

- palette:

  Palette name or colour vector.

- theme:

  A `ggplot2` theme; defaults to
  [`mp_theme_publication()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_publication.md).

- plot.title, plot.subtitle, plot.caption:

  Optional labels.

- font.size:

  Base font size.

- interactive:

  If `TRUE`, return an interactive `plotly` widget instead of a static
  `ggplot` (requires `plotly`). See
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md).

- verbose:

  Print the steps used.

## Value

A `ggplot` object, or a `plotly` widget when `interactive = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
mp_longitudinal(ps, time.var = "Week", subject.var = "SubjectID",
                group.by = "Arm", measure = "Shannon")
} # }
```
