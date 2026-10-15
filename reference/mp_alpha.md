# Plot alpha (within-sample) diversity

Computes one or more alpha-diversity indices from a
[phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
object and draws publication-ready box, violin or dot plots, optionally
annotated with pairwise statistical tests. Rarefaction curves are also
supported.

## Usage

``` r
mp_alpha(
  object,
  measure = c("Shannon", "Observed"),
  group.by = NULL,
  plot.type = c("box", "violin", "dot", "rarefaction"),
  test = c("wilcoxon", "kruskal", "t.test", "anova"),
  add.stats = TRUE,
  comparisons = NULL,
  p.adjust.method = "BH",
  facet.by = NULL,
  palette = NULL,
  point.size = 1.5,
  jitter = TRUE,
  theme = NULL,
  ...,
  plot.title = NULL,
  plot.subtitle = NULL,
  plot.caption = NULL,
  font.size = getOption("micropubr.font.size", 11),
  combine = TRUE,
  interactive = FALSE,
  return.data = FALSE,
  verbose = getOption("micropubr.verbose", FALSE)
)
```

## Arguments

- object:

  A [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
  object containing **count** data.

- measure:

  Character vector of indices. Any of `"Observed"`, `"Chao1"`, `"ACE"`,
  `"Shannon"`, `"Simpson"`, `"InvSimpson"`, `"Fisher"`, `"Pielou"`
  (evenness) and `"Faith"` (phylogenetic diversity; requires a tree and
  the `picante` package).

- group.by:

  Sample-metadata column to compare across (x-axis / fill). If `NULL`,
  all samples are pooled into a single group and no test is run.

- plot.type:

  One of `"box"`, `"violin"`, `"dot"` or `"rarefaction"`.

- test:

  Statistical test for group comparisons: `"wilcoxon"`, `"kruskal"`,
  `"t.test"` or `"anova"`. Pairwise post-hoc tests are Dunn (for
  `"kruskal"`) and Tukey HSD (for `"anova"`).

- add.stats:

  Add significance brackets / p-values when `group.by` has at least two
  levels.

- comparisons:

  Optional list of length-2 character vectors specifying which group
  pairs to test. If `NULL`, all pairwise comparisons are used.

- p.adjust.method:

  Multiple-testing correction (default `"BH"`).

- facet.by:

  Optional metadata column to facet each panel by.

- palette:

  Palette name, explicit colour vector, or `NULL` for the default
  colourblind-safe palette. See
  [`mp_palette()`](https://daniilidisk.github.io/micropubr/reference/mp_palette.md).

- point.size:

  Size of jittered sample points.

- jitter:

  Overlay individual sample points.

- theme:

  A `ggplot2` theme; defaults to the active theme (see
  [`mp_set_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md),
  otherwise
  [`mp_theme_publication()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_publication.md)).

- ...:

  Additional arguments passed to the primary geom (boxplot, violin or,
  for `plot.type = "dot"`, the point summary).

- plot.title, plot.subtitle, plot.caption:

  Optional labels. When `NULL`, sensible defaults are generated (e.g.
  the omnibus test result as subtitle).

- font.size:

  Base font size in points.

- combine:

  If multiple `measure`s are requested, combine the panels into a single
  [patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html)
  object (`TRUE`) or return a named list of `ggplot` objects (`FALSE`).

- interactive:

  If `TRUE`, return an interactive `plotly` widget instead of a static
  `ggplot` (requires the `plotly` package). Multiple panels are
  assembled with
  [`plotly::subplot()`](https://rdrr.io/pkg/plotly/man/subplot.html).
  Note that some static-only annotations (e.g. significance brackets)
  are not reproduced by the conversion.

- return.data:

  If `TRUE`, return the tidy data frame of per-sample diversity
  estimates and metadata instead of a plot.

- verbose:

  Print the methods used.

## Value

A `ggplot` object (single measure), a `patchwork` (multiple measures,
`combine = TRUE`) or a named list of `ggplot`s; a `plotly` widget when
`interactive = TRUE`; or a data frame when `return.data = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(GlobalPatterns, package = "phyloseq")
mp_alpha(GlobalPatterns, measure = c("Shannon", "Observed"),
         group.by = "SampleType")
} # }
```
