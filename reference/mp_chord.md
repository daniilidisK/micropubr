# Chord diagram of taxon-group relationships

Draws a circular chord diagram with `circlize` connecting the top taxa
at a chosen rank to the levels of a grouping variable, with ribbon
widths proportional to mean (relative) abundance. Like
[`mp_circos()`](https://daniilidisk.github.io/micropubr/reference/mp_circos.md)
this draws to the current graphics device (the exception to the ggplot
contract) and returns the underlying taxon-by-group matrix invisibly.

## Usage

``` r
mp_chord(
  object,
  group.by,
  tax.level = "Phylum",
  top.n = 12,
  relative = TRUE,
  palette = NULL,
  gap = 2,
  start.degree = 90,
  transparency = 0.5,
  directional = 0,
  title = NULL,
  verbose = getOption("micropubr.verbose", FALSE)
)
```

## Arguments

- object:

  A [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
  object.

- group.by:

  Sample-metadata column defining the groups (columns of the diagram).

- tax.level:

  Taxonomic rank to agglomerate to (rows of the diagram).

- top.n:

  Number of most-abundant taxa to show; the rest collapse into an
  "Other" sector.

- relative:

  Use relative abundance (`TRUE`) or raw counts.

- palette:

  Palette name or colour vector for the taxa sectors/ribbons.

- gap:

  Gap (degrees) between sectors.

- start.degree:

  Starting angle.

- transparency:

  Ribbon transparency in `[0, 1]`.

- directional:

  Passed to
  [`circlize::chordDiagram()`](https://rdrr.io/pkg/circlize/man/chordDiagram.html):
  `0` (undirected), `1` (taxa to groups) or `-1`.

- title:

  Optional plot title.

- verbose:

  Print the methods used.

## Value

Invisibly, the taxon-by-group abundance matrix that was drawn.

## Examples

``` r
if (FALSE) { # \dontrun{
data(GlobalPatterns, package = "phyloseq")
mp_chord(GlobalPatterns, group.by = "SampleType", tax.level = "Phylum")
} # }
```
