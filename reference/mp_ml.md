# Machine-learning classification and regression on microbiome features

Trains a random forest (`"rf"`) or gradient-boosted trees (`"xgb"`)
model on the (optionally transformed) feature table of a
[phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
object, using a sample-metadata column as the response. The task is
detected automatically from the response type - **classification** for a
factor/character column, **regression** for a numeric one - and can be
forced with `task`. Out-of-fold predictions from k-fold cross-validation
drive the performance metrics; feature importance comes from a model fit
on all samples.

## Usage

``` r
mp_ml(
  object,
  group.by = NULL,
  method = c("rf", "xgb"),
  task = c("auto", "classification", "regression"),
  tax.level = NULL,
  top.n = NULL,
  transform = c("clr", "relative", "none"),
  plot.type = c("importance", "roc", "confusion", "prediction"),
  n.features = 20,
  cv.folds = 5,
  ntree = 500,
  nrounds = 150,
  seed = NULL,
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
  object with counts.

- group.by:

  Sample-metadata column used as the response/outcome.

- method:

  Learner: `"rf"` (randomForest) or `"xgb"` (xgboost).

- task:

  `"auto"` (default), `"classification"` or `"regression"`.

- tax.level:

  Optional rank to agglomerate features to before modelling.

- top.n:

  Optional number of most-abundant features to keep (dimension
  reduction); `NULL` keeps them all.

- transform:

  Feature transform applied before modelling: `"clr"` (centred
  log-ratio, the default), `"relative"` or `"none"`.

- plot.type:

  What to draw: `"importance"` (ranked feature importance, the default),
  `"roc"` (ROC curve, two-class classification only), `"confusion"`
  (confusion-matrix heatmap, classification) or `"prediction"` (observed
  vs predicted scatter, regression).

- n.features:

  Number of top features shown in the importance plot.

- cv.folds:

  Number of cross-validation folds.

- ntree:

  Number of trees for the random forest.

- nrounds:

  Number of boosting rounds for xgboost.

- seed:

  Optional integer seed for reproducible folds/models.

- palette:

  Palette name or colour vector.

- theme:

  A `ggplot2` theme; defaults to the active theme (see
  [`mp_set_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md),
  otherwise
  [`mp_theme_publication()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_publication.md)).

- ...:

  Additional arguments passed to the primary geom of the chosen plot
  (e.g.
  [`ggplot2::geom_col()`](https://ggplot2.tidyverse.org/reference/geom_bar.html)
  for importance).

- plot.title, plot.subtitle, plot.caption:

  Optional labels.

- font.size:

  Base font size.

- interactive:

  If `TRUE`, return an interactive `plotly` widget instead of a static
  `ggplot` (requires `plotly`). See
  [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md).

- return.data:

  If `TRUE`, return the tidy feature-importance data frame instead of a
  plot, with the fitted model (`attr(., "model")`), the cross-validation
  predictions (`attr(., "cv")`) and a named list of performance metrics
  (`attr(., "metrics")`) attached as attributes.

- verbose:

  Print the method and performance summary.

## Value

A `ggplot` object, a `plotly` widget when `interactive = TRUE`, or a
data frame when `return.data = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(GlobalPatterns, package = "phyloseq")
gp <- phyloseq::tax_glom(GlobalPatterns, "Genus")
# Classification: which genera discriminate sample types?
mp_ml(gp, group.by = "SampleType", method = "rf")
# Two-class ROC curve.
soil_vs_feces <- phyloseq::subset_samples(
  gp, SampleType %in% c("Soil", "Feces"))
mp_ml(soil_vs_feces, group.by = "SampleType", plot.type = "roc")
} # }
```
