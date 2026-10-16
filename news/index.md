# Changelog

## micropubr 0.3.0

### Differential abundance

- [`mp_run_da()`](https://daniilidisk.github.io/micropubr/reference/mp_run_da.md)
  is rebuilt around one shared design for every method:
  - several methods in one call
    (`method = c("wilcox", "lm", "ancombc2")`), returned as one long
    table with `method` and `contrast` columns;
  - multi-group designs: every level against `reference`, or all pairs
    with `pairwise = TRUE` (previously only two groups);
  - covariate adjustment with `adjust` and repeated measures with
    `random`;
  - a common prevalence filter (`prevalence = 0.1`) and multiple-testing
    correction (`p.adjust.method`);
  - new engines `"lm"` (linear model on CLR values; linear mixed model
    via `lmerTest` with `random`) and `"ancombc2"` (ANCOM-BC2, with
    structural zeros and the pseudo-count sensitivity flag in the new
    `robust` column). `"ancombc"` is now an alias of `"ancombc2"`;
  - DESeq2 now uses `poscounts` size factors, so sparse tables no longer
    fail.
- [`mp_diff()`](https://daniilidisk.github.io/micropubr/reference/mp_diff.md)
  gains `plot.type = "heatmap"` (fold changes across contrasts),
  `"consensus"` (agreement across methods), and `"zeros"` (structural
  zeros), plus `contrast`, `robust.only` and `top.n` arguments.
  Multi-contrast and multi-method results are faceted. ANCOM-BC2 results
  that fail its pseudo-count sensitivity analysis are not counted as
  significant by default (`robust.only = TRUE`, matching ANCOM-BC2’s
  `diff_robust`): in a confounding simulation on real counts nearly all
  of ANCOM-BC2’s false positives failed that analysis. A named list of
  results is stacked for comparison, and `"concordance"` works on a
  single two-method result.
- Forest and dot plots now draw the `top.n` most significant features
  (default 30; per panel when faceted) instead of every significant one,
  and say in the caption how many are shown; their legends no longer
  list “NS”.

### Bug fixes

- ANCOM-BC2 results were plotted from the intercept column
  (`lfc_(Intercept)`) rather than the group effect, and on a log2 axis
  although ANCOM-BC estimates are natural logarithms. Both are fixed,
  and covariate terms are no longer mistaken for group contrasts.
- Results from the original
  [`ANCOMBC::ancombc()`](https://rdrr.io/pkg/ANCOMBC/man/ancombc.html)
  (whose `res` is a list of tables) could not be parsed;
  `mp_run_da(method = "ancombc")` failed for the same reason and also
  ignored the requested reference group.
- ALDEx2 results were plotted using the standardised `effect` size
  instead of the log2 difference `diff.btw`.
- [`mp_ordination()`](https://daniilidisk.github.io/micropubr/reference/mp_ordination.md)
  drew ellipses for groups of three samples, which
  [`stat_ellipse()`](https://ggplot2.tidyverse.org/reference/stat_ellipse.html)
  cannot estimate; it now needs four.
- Composition plots
  ([`mp_taxa()`](https://daniilidisk.github.io/micropubr/reference/mp_taxa.md),
  [`mp_chord()`](https://daniilidisk.github.io/micropubr/reference/mp_chord.md),
  [`mp_longitudinal()`](https://daniilidisk.github.io/micropubr/reference/mp_longitudinal.md))
  folded abundant taxa that are unassigned at the plotted rank into
  “Other”, so fewer than `top.n` taxa were shown. They are now kept and
  labelled by their deepest assigned rank (e.g. “Unassigned
  Lachnospiraceae”).
- Group-comparison p-values (Wilcoxon, t-test, Kruskal-Wallis, ANOVA,
  Tukey HSD) were taken from `rstatix`, which rounds them to three
  significant digits before micropubr adjusted and printed them. They
  are now computed exactly with `stats`.
- PCoA axis percentages are now relative to the sum of all eigenvalues,
  as in [`ape::pcoa()`](https://rdrr.io/pkg/ape/man/pcoa.html) and
  `phyloseq`, so labels match when negative eigenvalues occur.
- [`mp_dispersion()`](https://daniilidisk.github.io/micropubr/reference/mp_dispersion.md)
  rotates long group labels, and
  [`mp_taxa()`](https://daniilidisk.github.io/micropubr/reference/mp_taxa.md)
  wraps long facet labels instead of truncating them.
- [`mp_export()`](https://daniilidisk.github.io/micropubr/reference/mp_export.md)
  wrote no PDF or EPS file on macOS without XQuartz, because the cairo
  devices it used need XQuartz there. It now falls back to the Quartz
  PDF device on macOS and to
  [`pdf()`](https://rdrr.io/r/grDevices/pdf.html)/[`postscript()`](https://rdrr.io/r/grDevices/postscript.html)
  elsewhere.
- README: the
  [`mp_import_qiime2()`](https://daniilidisk.github.io/micropubr/reference/mp_import_qiime2.md)
  example used a non-existent argument.

### Reporting and reproducibility

- New
  [`mp_methods()`](https://daniilidisk.github.io/micropubr/reference/mp_methods.md):
  a methods paragraph, software versions, references and statistical
  notes, generated from the provenance every figure now records.
- New
  [`mp_alt_text()`](https://daniilidisk.github.io/micropubr/reference/mp_alt_text.md):
  accessible alt text for any micropubr figure.
- Statistical guardrails flag small groups, richness compared across
  very uneven sequencing depths, PERMANOVA with unequal dispersions
  (PERMDISP is now run alongside every PERMANOVA in
  [`mp_ordination()`](https://daniilidisk.github.io/micropubr/reference/mp_ordination.md)),
  methods that ignored the requested covariates, and covariates so
  strongly associated with the grouping (Cramer’s V \>= 0.7, or
  R-squared \>= 0.5) that adjusting for them is unstable. Control them
  with `mp_set_defaults(guardrails = "warn" | "caption" | "off")`.

### Themes

- New
  [`mp_darken()`](https://daniilidisk.github.io/micropubr/reference/mp_darken.md)
  recolours a finished figure for a dark background while keeping its
  layout. Adding
  [`mp_theme_dark()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_dark.md)
  to a finished plot replaces the whole theme, which undid rotated axis
  labels and facet strips and left labels and significance brackets
  drawn in black invisible.
- Composition and differential-abundance labels no longer treat
  placeholders such as MGnify’s “Root” or QIIME’s “Unassigned” as taxon
  names.
- `mp_methods(format = "list")` now returns the paragraph without its
  software sentence, plus `versions` and `citations`, so figures
  described one at a time can share a single software sentence.

### Interoperability

- Every `mp_*()` function accepts a `TreeSummarizedExperiment` (or
  `SummarizedExperiment`), and the new
  [`mp_as_phyloseq()`](https://daniilidisk.github.io/micropubr/reference/mp_as_phyloseq.md)
  converts one explicitly, including its row tree.

### Infrastructure

- Statistical equivalence tests against `vegan`, `phyloseq` and `stats`,
  fixture-based tests for real ANCOM-BC/ANCOM-BC2 output, and `vdiffr`
  visual regression tests.
- GitHub Actions for R CMD check, test coverage, BiocCheck and the
  pkgdown site; `inst/CITATION` and `codemeta.json`.

## micropubr 0.2.0

- New `interactive` argument on every ggplot-based plotting function.
  With `interactive = TRUE` they return a `plotly` widget instead of a
  static `ggplot`: hover tooltips, zoom/pan and legend toggling, with
  the ggplot subtitle and caption preserved. Multi-panel output
  (multiple measures or
  [`mp_combine()`](https://daniilidisk.github.io/micropubr/reference/mp_combine.md)
  figures) is assembled into a single
  [`plotly::subplot()`](https://rdrr.io/pkg/plotly/man/subplot.html)
  with a deduplicated shared legend. Supported by
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md),
  [`mp_ordination()`](https://daniilidisk.github.io/micropubr/reference/mp_ordination.md),
  [`mp_taxa()`](https://daniilidisk.github.io/micropubr/reference/mp_taxa.md),
  [`mp_diff()`](https://daniilidisk.github.io/micropubr/reference/mp_diff.md),
  [`mp_longitudinal()`](https://daniilidisk.github.io/micropubr/reference/mp_longitudinal.md),
  [`mp_function()`](https://daniilidisk.github.io/micropubr/reference/mp_function.md),
  [`mp_network()`](https://daniilidisk.github.io/micropubr/reference/mp_network.md),
  [`mp_tree()`](https://daniilidisk.github.io/micropubr/reference/mp_tree.md),
  [`mp_dispersion()`](https://daniilidisk.github.io/micropubr/reference/mp_dispersion.md),
  [`mp_mantel()`](https://daniilidisk.github.io/micropubr/reference/mp_mantel.md),
  [`mp_varpart()`](https://daniilidisk.github.io/micropubr/reference/mp_varpart.md),
  [`mp_prevalence()`](https://daniilidisk.github.io/micropubr/reference/mp_prevalence.md),
  [`mp_indicator()`](https://daniilidisk.github.io/micropubr/reference/mp_indicator.md),
  [`mp_ml()`](https://daniilidisk.github.io/micropubr/reference/mp_ml.md)
  and
  [`mp_combine()`](https://daniilidisk.github.io/micropubr/reference/mp_combine.md).
  The circular `circlize`-based figures
  ([`mp_chord()`](https://daniilidisk.github.io/micropubr/reference/mp_chord.md),
  [`mp_circos()`](https://daniilidisk.github.io/micropubr/reference/mp_circos.md),
  [`mp_coverage()`](https://daniilidisk.github.io/micropubr/reference/mp_coverage.md))
  are not ggplots and have no interactive mode. Requires the new
  (suggested) `plotly` package.
- New
  [`mp_qc_profile()`](https://daniilidisk.github.io/micropubr/reference/mp_qc_profile.md):
  per-cycle FASTQ read-quality profiles (the diagnostic used to pick
  DADA2/QIIME 2 truncation lengths) as a frequency heatmap or a mean ±
  inter-quartile ribbon, with mean/quartile overlays and a Q-score
  reference line. Reads FASTQ via the lightweight `ShortRead` backend,
  or accepts a pre-computed per-cycle data frame, and supports
  `interactive` and `return.data` like the rest of the family.
- [`mp_export()`](https://daniilidisk.github.io/micropubr/reference/mp_export.md)
  now saves an interactive widget to a self-contained `.html` page (via
  `htmlwidgets`) when given one.
- Both vignettes gained a live “Interactive figures” section.
- New
  [`mp_ml()`](https://daniilidisk.github.io/micropubr/reference/mp_ml.md):
  supervised machine learning (random forest and xgboost) for microbiome
  classification *and* regression, with k-fold cross-validation and
  `"importance"`, `"roc"`, `"confusion"` and `"prediction"` plots. Like
  the other functions it returns a tidy importance table (with the
  fitted model and performance metrics attached) when
  `return.data = TRUE`.
- New
  [`mp_chord()`](https://daniilidisk.github.io/micropubr/reference/mp_chord.md):
  chord diagram linking the top taxa to the levels of a grouping
  variable (via `circlize`).
- [`mp_circos()`](https://daniilidisk.github.io/micropubr/reference/mp_circos.md)
  gained colour/weight-mapped links (`link.by`), multi-sample coverage
  tracks, a real heatmap track, feature-type annotation colours
  (`annotation.by`), configurable track heights and an optional legend.
- x-axis tick labels now rotate automatically when categories are many
  or long
  ([`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md),
  [`mp_longitudinal()`](https://daniilidisk.github.io/micropubr/reference/mp_longitudinal.md)),
  so they no longer overlap.
- [`mp_tree()`](https://daniilidisk.github.io/micropubr/reference/mp_tree.md)
  no longer prints the harmless “more than one class ‘phylo’ in cache”
  diagnostic emitted when `phyloseq` and `tidytree`/`ggtree` are loaded
  together.
- Two new vignettes: “Visualization gallery” and “Statistics, machine
  learning and genomics”.

## micropubr 0.1.0

- Initial release.
- Core plotting family:
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md),
  [`mp_ordination()`](https://daniilidisk.github.io/micropubr/reference/mp_ordination.md),
  [`mp_taxa()`](https://daniilidisk.github.io/micropubr/reference/mp_taxa.md),
  [`mp_diff()`](https://daniilidisk.github.io/micropubr/reference/mp_diff.md),
  [`mp_tree()`](https://daniilidisk.github.io/micropubr/reference/mp_tree.md),
  [`mp_network()`](https://daniilidisk.github.io/micropubr/reference/mp_network.md),
  [`mp_longitudinal()`](https://daniilidisk.github.io/micropubr/reference/mp_longitudinal.md)
  and
  [`mp_function()`](https://daniilidisk.github.io/micropubr/reference/mp_function.md).
- Infrastructure:
  [`mp_palette()`](https://daniilidisk.github.io/micropubr/reference/mp_palette.md),
  the `mp_theme_*()` family,
  [`mp_combine()`](https://daniilidisk.github.io/micropubr/reference/mp_combine.md),
  [`mp_export()`](https://daniilidisk.github.io/micropubr/reference/mp_export.md),
  the `mp_import_*()` helpers and the interactive
  [`mp_explore()`](https://daniilidisk.github.io/micropubr/reference/mp_explore.md)
  Shiny gadget.
- Colourblind-safe default palettes and a clean publication theme
  applied across all functions.
- Automatic statistical overlays (pairwise tests, PERMANOVA) and
  reproducible, auto-generated figure captions.
