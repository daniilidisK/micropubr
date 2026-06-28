# micropubr 0.2.0

* New `interactive` argument on every ggplot-based plotting function. With
  `interactive = TRUE` they return a `plotly` widget instead of a static
  `ggplot`: hover tooltips, zoom/pan and legend toggling, with the ggplot
  subtitle and caption preserved. Multi-panel output (multiple measures or
  `mp_combine()` figures) is assembled into a single `plotly::subplot()` with a
  deduplicated shared legend. Supported by `mp_alpha()`, `mp_ordination()`,
  `mp_taxa()`, `mp_diff()`, `mp_longitudinal()`, `mp_function()`,
  `mp_network()`, `mp_tree()`, `mp_dispersion()`, `mp_mantel()`, `mp_varpart()`,
  `mp_prevalence()`, `mp_indicator()`, `mp_ml()` and `mp_combine()`. The
  circular `circlize`-based figures (`mp_chord()`, `mp_circos()`,
  `mp_coverage()`) are not ggplots and have no interactive mode. Requires the
  new (suggested) `plotly` package.
* `mp_export()` now saves an interactive widget to a self-contained `.html`
  page (via `htmlwidgets`) when given one.
* Both vignettes gained a live "Interactive figures" section.
* New `mp_ml()`: supervised machine learning (random forest and xgboost) for
  microbiome classification *and* regression, with k-fold cross-validation and
  `"importance"`, `"roc"`, `"confusion"` and `"prediction"` plots. Like the
  other functions it returns a tidy importance table (with the fitted model and
  performance metrics attached) when `return.data = TRUE`.
* New `mp_chord()`: chord diagram linking the top taxa to the levels of a
  grouping variable (via `circlize`).
* `mp_circos()` gained colour/weight-mapped links (`link.by`), multi-sample
  coverage tracks, a real heatmap track, feature-type annotation colours
  (`annotation.by`), configurable track heights and an optional legend.
* x-axis tick labels now rotate automatically when categories are many or long
  (`mp_alpha()`, `mp_longitudinal()`), so they no longer overlap.
* `mp_tree()` no longer prints the harmless "more than one class 'phylo' in
  cache" diagnostic emitted when `phyloseq` and `tidytree`/`ggtree` are loaded
  together.
* Two new vignettes: "Visualization gallery" and "Statistics, machine learning
  and genomics".

# micropubr 0.1.0

* Initial release.
* Core plotting family: `mp_alpha()`, `mp_ordination()`, `mp_taxa()`,
  `mp_diff()`, `mp_tree()`, `mp_network()`, `mp_longitudinal()` and
  `mp_function()`.
* Infrastructure: `mp_palette()`, the `mp_theme_*()` family, `mp_combine()`,
  `mp_export()`, the `mp_import_*()` helpers and the interactive `mp_explore()`
  Shiny gadget.
* Colourblind-safe default palettes and a clean publication theme applied
  across all functions.
* Automatic statistical overlays (pairwise tests, PERMANOVA) and
  reproducible, auto-generated figure captions.
