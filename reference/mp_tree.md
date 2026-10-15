# Plot an annotated phylogenetic tree

Renders the phylogenetic tree stored in a
[phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
object using `ggtree`, with optional tip colouring by taxonomy and
abundance rings drawn with `ggtreeExtra`. To keep large trees legible,
agglomerate with `tax.level` and/or restrict to the most abundant tips
with `top.n`.

## Usage

``` r
mp_tree(
  object,
  layout = c("rectangular", "circular", "fan", "radial"),
  annotate = c("abundance", "taxonomy"),
  color.by = NULL,
  group.by = NULL,
  tax.level = NULL,
  top.n = NULL,
  tip.label = FALSE,
  tip.label.size = 2,
  highlight.clades = NULL,
  open.angle = 0,
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
  object containing a `phy_tree`.

- layout:

  Tree layout: `"rectangular"`, `"circular"`, `"fan"` or `"radial"`.

- annotate:

  Annotation layers to draw: any of `"abundance"` and `"taxonomy"`.

- color.by:

  Taxonomic rank used to colour tip points (e.g. `"Phylum"`).

- group.by:

  Metadata column; when set, the abundance ring becomes a per-group
  heatmap instead of a single bar.

- tax.level:

  Optional rank to agglomerate to before plotting.

- top.n:

  Optional number of most-abundant tips to retain.

- tip.label:

  Show tip labels.

- tip.label.size:

  Tip-label font size.

- highlight.clades:

  Named list mapping a taxonomic rank to one or more values whose clades
  should be highlighted, e.g. `list(Phylum = "Firmicutes")`.

- open.angle:

  Opening angle (degrees) for circular/fan layouts.

- palette:

  Palette name or colour vector.

- theme:

  A `ggtree` theme; defaults to
  [`ggtree::theme_tree()`](https://rdrr.io/pkg/ggtree/man/theme_tree.html).

- plot.title, plot.subtitle, plot.caption:

  Optional labels.

- font.size:

  Base font size.

- interactive:

  If `TRUE`, attempt to return an interactive `plotly` widget (requires
  `plotly`). Experimental for trees: `ggtree` layers often do not
  translate, in which case the static plot is returned with a warning.
  See
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md).

- verbose:

  Print the steps used.

## Value

A `ggtree`/`ggplot` object, or a `plotly` widget when
`interactive = TRUE` and the conversion succeeds.

## Examples

``` r
if (FALSE) { # \dontrun{
data(GlobalPatterns, package = "phyloseq")
gp <- phyloseq::tax_glom(GlobalPatterns, "Phylum")
mp_tree(gp, layout = "circular", color.by = "Phylum",
        group.by = "SampleType")
} # }
```
