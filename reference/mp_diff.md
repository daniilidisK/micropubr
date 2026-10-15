# Plot differential-abundance results

Visualises differential-abundance results as a volcano, MA, forest or
ranked-dot plot, as a heatmap of fold changes across contrasts, or –
when several methods were run – as a consensus matrix showing where the
methods agree. Results from
[`mp_run_da()`](https://daniilidisk.github.io/micropubr/reference/mp_run_da.md),
`DESeq2`, `ALDEx2` and `ANCOM-BC`/`ANCOM-BC2` are parsed automatically;
any other tool can be supplied as a data frame with explicit column
mappings, and a named list of results from different tools is stacked
for comparison.

## Usage

``` r
mp_diff(
  result,
  plot.type = c("volcano", "ma", "forest", "dot", "heatmap", "consensus", "zeros",
    "concordance"),
  fc.threshold = 1,
  p.threshold = 0.05,
  y.value = c("padj", "pvalue"),
  label.top = 10,
  color.up = "#B2182B",
  color.down = "#2166AC",
  color.ns = "grey75",
  result2 = NULL,
  method.names = c("Method 1", "Method 2"),
  contrast = NULL,
  robust.only = TRUE,
  top.n = 30,
  feature = NULL,
  log2fc = NULL,
  pvalue = NULL,
  padj = NULL,
  mean = NULL,
  lfcSE = NULL,
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

- result:

  A differential-abundance result: the output of
  [`mp_run_da()`](https://daniilidisk.github.io/micropubr/reference/mp_run_da.md),
  a `DESeqResults` object, an `ALDEx2` data frame, an `ANCOM-BC(2)`
  result list, a generic data frame (see
  `feature`/`log2fc`/`pvalue`/`padj`), or a named list of any of these
  (one per method, e.g. `list(DESeq2 = res1, ALDEx2 = res2)`).

- plot.type:

  One of `"volcano"`, `"ma"`, `"forest"`, `"dot"`, `"heatmap"` (fold
  change of significant taxa across contrasts), `"consensus"` (agreement
  across methods), `"zeros"` (ANCOM-BC2 structural zeros) or
  `"concordance"` (fold changes of two methods against each other).

- fc.threshold:

  Absolute log2 fold-change cut-off for significance.

- p.threshold:

  Adjusted p-value cut-off.

- y.value:

  For volcano plots, plot `-log10()` of `"padj"` (default) or `"pvalue"`
  on the y-axis.

- label.top:

  Number of most-significant features to label per panel.

- color.up, color.down, color.ns:

  Colours for enriched, depleted and non-significant features.

- result2:

  A second result (same parsing rules) for `plot.type = "concordance"`.
  Not needed when `result` already holds exactly two methods.

- method.names:

  Length-2 labels for the two methods in a concordance plot.

- contrast:

  Optional character vector of contrasts (as in the `contrast` column,
  e.g. `"Treatment vs Control"`, or an ANCOM-BC2 term such as
  `"groupTreatment"`) to keep.

- robust.only:

  Treat significant features that fail ANCOM-BC2's pseudo-count
  sensitivity analysis as not significant (the default, which matches
  ANCOM-BC2's `diff_robust`). Results without sensitivity information
  are unaffected.

- top.n:

  Maximum number of features in forest, dot, heatmap and consensus plots
  (the most significant are kept; per panel for forest and dot plots).

- feature, log2fc, pvalue, padj, mean, lfcSE:

  Optional column-name overrides when `result` is a generic data frame.

- theme:

  A `ggplot2` theme; defaults to the active theme (see
  [`mp_set_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md),
  otherwise
  [`mp_theme_publication()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_publication.md)).

- ...:

  Additional arguments passed to the point geom.

- plot.title, plot.subtitle, plot.caption:

  Optional labels.

- font.size:

  Base font size.

- interactive:

  If `TRUE`, return an interactive `plotly` widget instead of a static
  `ggplot` (requires `plotly`). See
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md).

- return.data:

  If `TRUE`, return the parsed differential-abundance data frame (with a
  `direction` column) instead of a plot; for `plot.type = "consensus"`,
  the per-feature consensus table.

- verbose:

  Print the detected result type and counts.

## Value

A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
data frame when `return.data = TRUE`.

## Details

Results with several contrasts (multi-group designs) or several methods
are drawn as facets. ANCOM-BC(2) fold changes, which are natural
logarithms, are converted to log2. Taxa that ANCOM-BC2 classifies as
structural zeros (absent from a whole group) are not tested; they are
counted in the subtitle and shown by `plot.type = "zeros"`. When results
carry ANCOM-BC2's pseudo-count sensitivity analysis, taxa that fail it
are not counted as significant (as in ANCOM-BC2's own `diff_robust`
output); they are counted in the subtitle and marked in heatmaps. With
`robust.only = FALSE` they are kept and drawn as open symbols.

## See also

[`mp_run_da()`](https://daniilidisk.github.io/micropubr/reference/mp_run_da.md)
to run the analysis.

## Examples

``` r
df <- data.frame(
  feature = paste0("ASV", 1:50),
  log2FC = rnorm(50, sd = 2),
  pvalue = runif(50)^2,
  padj = p.adjust(runif(50)^2, "BH")
)
mp_diff(df, plot.type = "volcano")
```
