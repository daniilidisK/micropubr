# Plot a taxon co-occurrence network

Builds a correlation network between taxa and renders it with `ggraph`.
Correlations are computed on centred log-ratio (CLR) coordinates by
default; edges are filtered by correlation strength and adjusted
p-value. SparCC / SPIEC-EASI inference is available when the `SpiecEasi`
package is installed.

## Usage

``` r
mp_network(
  object,
  method = c("spearman", "pearson", "sparcc", "spiec-easi"),
  tax.level = NULL,
  top.n = 50,
  cor.threshold = 0.6,
  p.threshold = 0.05,
  p.adjust.method = "BH",
  layout = c("fr", "stress", "kk", "circle"),
  node.color.by = NULL,
  node.size.by = c("abundance", "degree"),
  edge.color.by = c("sign", "none"),
  label.nodes = FALSE,
  label.size = 2.5,
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

- method:

  Association method: `"spearman"`, `"pearson"`, `"sparcc"` or
  `"spiec-easi"`. The latter two require the `SpiecEasi` package.

- tax.level:

  Optional rank to agglomerate to before building the network.

- top.n:

  Keep only the `top.n` most-abundant taxa (for tractability).

- cor.threshold:

  Minimum absolute correlation to draw an edge.

- p.threshold:

  Maximum adjusted p-value for an edge (correlation methods only).

- p.adjust.method:

  Multiple-testing correction (default `"BH"`).

- layout:

  `ggraph` layout: `"fr"`, `"stress"`, `"kk"` or `"circle"`.

- node.color.by:

  Taxonomic rank used to colour nodes.

- node.size.by:

  Size nodes by `"abundance"` or `"degree"`.

- edge.color.by:

  Colour edges by `"sign"` (default) or `"none"`.

- label.nodes:

  Label nodes with their taxon name.

- label.size:

  Node-label font size.

- palette:

  Palette name or colour vector.

- theme:

  A `ggplot2` theme; defaults to a clean void theme.

- plot.title, plot.subtitle, plot.caption:

  Optional labels.

- font.size:

  Base font size.

- interactive:

  If `TRUE`, return an interactive `plotly` widget instead of a static
  `ggplot` (requires `plotly`). Experimental for network plots: `ggraph`
  edge geoms translate only partially. See
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md).

- verbose:

  Print the steps used.

## Value

A `ggplot`/`ggraph` object, or a `plotly` widget when
`interactive = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(GlobalPatterns, package = "phyloseq")
gp <- phyloseq::tax_glom(GlobalPatterns, "Genus")
mp_network(gp, top.n = 60, node.color.by = "Phylum")
} # }
```
