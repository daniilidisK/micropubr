# Circular (Circos-style) genome plot

Draws a circular genome plot with `circlize`: an ideogram of contigs
plus optional coverage, GC, heatmap and annotation tracks and ribbon
links (e.g. synteny or correlations). Unlike the other `mp_*()`
functions this draws with base graphics to the current device (the one
exception to the ggplot contract) and returns the assembled track data
invisibly.

## Usage

``` r
mp_circos(
  genome,
  coverage = NULL,
  gc = NULL,
  annotations = NULL,
  links = NULL,
  heatmap = NULL,
  coverage.type = c("line", "area", "heatmap"),
  palette = NULL,
  coverage.color = "#0072B2",
  gc.color = "#009E73",
  link.by = NULL,
  link.palette = NULL,
  annotation.by = NULL,
  annotation.palette = NULL,
  heatmap.palette = NULL,
  track.heights = NULL,
  gap = 2,
  start.degree = 90,
  legend = TRUE,
  title = NULL,
  verbose = getOption("micropubr.verbose", FALSE)
)
```

## Arguments

- genome:

  A data frame describing contigs with columns `contig` and `length`, or
  a named numeric vector of lengths.

- coverage:

  Optional coverage data frame with `contig`, `position` (or `start`)
  and `value`, plus an optional `sample` column; or a named list of such
  data frames (one per sample).

- gc:

  Optional GC-content/skew data frame with `contig`, `start`, `end` and
  `value`.

- annotations:

  Optional features with `contig`, `start`, `end` and optional `label`.

- links:

  Optional ribbon links with `contig1`, `start1`, `end1`, `contig2`,
  `start2`, `end2` and an optional weight column (see `link.by`).

- heatmap:

  Optional data frame for a heatmap track with `contig`, `position` (or
  `start`), `value` and an optional `row` column (one band per row, e.g.
  per sample).

- coverage.type:

  `"line"`, `"area"` or `"heatmap"` (single-sample only).

- palette:

  Palette name or colour vector for contigs / coverage samples.

- coverage.color:

  Colour for a single coverage track.

- gc.color:

  Colour for the GC track.

- link.by:

  Optional numeric column in `links` used to colour and weight the
  ribbons.

- link.palette:

  Palette for `link.by` (defaults to diverging/sequential depending on
  the sign of the values).

- annotation.by:

  Optional column in `annotations` (a feature type) used to colour the
  annotation rectangles.

- annotation.palette:

  Palette for `annotation.by`.

- heatmap.palette:

  Sequential palette for the heatmap track.

- track.heights:

  Optional named list overriding default track heights (`ideogram`,
  `coverage`, `gc`, `heatmap`, `annotations`).

- gap:

  Gap (degrees) between contigs.

- start.degree:

  Starting angle.

- legend:

  Draw a legend describing the coloured tracks.

- title:

  Optional plot title.

- verbose:

  Print the methods used.

## Value

Invisibly, a list of the track data used (so `return.data`-style
workflows still work). Called for the side effect of drawing the plot.

## Details

Coverage can hold a `sample` column (or be a named list of data frames),
in which case one coloured line is drawn per sample. Links can be
coloured and weighted by a numeric column via `link.by` (a diverging
ramp is used when the values are signed, a sequential ramp otherwise),
and annotations can be coloured by feature type via `annotation.by`.

## Examples

``` r
if (FALSE) { # \dontrun{
genome <- data.frame(contig = c("c1", "c2"), length = c(5000, 3000))
cov <- data.frame(contig = "c1", position = seq(1, 5000, 50),
                  value = runif(100, 10, 100))
links <- data.frame(contig1 = "c1", start1 = 100, end1 = 400,
                    contig2 = "c2", start2 = 800, end2 = 1200,
                    weight = 0.8)
mp_circos(genome, coverage = cov, links = links, link.by = "weight")
} # }
```
