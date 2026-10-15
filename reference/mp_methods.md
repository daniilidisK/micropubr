# Write a methods paragraph for micropubr figures

Assembles a ready-to-edit methods paragraph describing how the supplied
figures were made: the statistical tests and their settings, the
transformations, filters and multiple-testing corrections, sample sizes,
software versions and the references to cite. It is built from the
provenance each `mp_*()` function attaches to its plot, so it always
matches what was actually computed. Statistical caveats raised while the
figures were made (see the guardrails in
[`mp_set_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md))
are returned as notes.

## Usage

``` r
mp_methods(..., format = c("text", "markdown", "list"), references = TRUE)
```

## Arguments

- ...:

  Figures made by `mp_*()` functions: `ggplot`/`patchwork` objects,
  interactive widgets, or lists of them. With no figures, every figure
  made in the current R session is described.

- format:

  `"text"` (plain text with references), `"markdown"`, or `"list"` (the
  components, for programmatic use).

- references:

  Include the reference list.

## Value

For `"text"` and `"markdown"`, a character string of class `mp_methods`
(printed nicely). For `"list"`, the components: `text` (the methods
paragraph, without the software sentence), `software` (that sentence),
`versions` (named package versions, including R), `citations` (in-text
citation labels, named by reference key), `references` and `notes`
(statistical caveats). Keeping the software apart lets a caller that
describes figures one at a time write a single software sentence for all
of them.

## Details

The paragraph covers the statistical-analysis and software items that
reporting guidelines such as STORMS (Mirzayi et al. 2021, *Nature
Medicine* 27:1885) ask for. It does not know your study design, sampling
or sequencing, so treat it as a draft to complete, not a finished
section.

## See also

[`mp_alt_text()`](https://daniilidisk.github.io/micropubr/reference/mp_alt_text.md)
for figure alt text.

## Examples

``` r
# \donttest{
data(GlobalPatterns, package = "phyloseq")
p1 <- mp_alpha(GlobalPatterns, measure = "Shannon", group.by = "SampleType",
               test = "kruskal")
#> Warning: Group(s) Feces (n = 4), Freshwater (n = 2), Freshwater (creek) (n = 3), Mock (n = 3), Ocean (n = 3), Sediment (estuary) (n = 3), Skin (n = 3), Soil (n = 3), Tongue (n = 2) have fewer than 5 samples; group comparisons have little power and their p-values are unstable at this size.
p2 <- mp_ordination(GlobalPatterns, group.by = "SampleType")
#> Warning: Group(s) Feces (n = 4), Freshwater (n = 2), Freshwater (creek) (n = 3), Mock (n = 3), Ocean (n = 3), Sediment (estuary) (n = 3), Skin (n = 3), Soil (n = 3), Tongue (n = 2) have fewer than 5 samples; group comparisons have little power and their p-values are unstable at this size.
#> Warning: Group dispersions differ (PERMDISP p = 0.002); a significant PERMANOVA may reflect differences in spread rather than in location (centroids).
mp_methods(p1, p2)
#> Alpha diversity (Shannon) was estimated for 26 samples with phyloseq (McMurdie
#>   and Holmes 2013). Differences among the 9 SampleType groups (Feces, n = 4;
#>   Freshwater, n = 2; Freshwater (creek), n = 3; Mock, n = 3; Ocean, n = 3;
#>   Sediment (estuary), n = 3; Skin, n = 3; Soil, n = 3; Tongue, n = 2) were
#>   tested with the Kruskal-Wallis test (Kruskal and Wallis 1952), followed by
#>   pairwise Dunn's tests (Dunn 1964) adjusted with the Benjamini-Hochberg
#>   procedure (Benjamini and Hochberg 1995) (rstatix; Kassambara 2026).
#>   Between-sample (beta) diversity of 26 samples was visualised by principal
#>   coordinates analysis (PCoA) of Bray-Curtis dissimilarity computed with
#>   phyloseq (McMurdie and Holmes 2013). Differences in community composition
#>   among SampleType groups (Feces, n = 4; Freshwater, n = 2; Freshwater (creek),
#>   n = 3; Mock, n = 3; Ocean, n = 3; Sediment (estuary), n = 3; Skin, n = 3;
#>   Soil, n = 3; Tongue, n = 2) were tested by permutational multivariate
#>   analysis of variance (PERMANOVA; vegan::adonis2, 999 permutations; Anderson
#>   2001): R² = 0.676, F = 4.43, p = 0.001. Homogeneity of multivariate
#>   dispersions, which PERMANOVA assumes, was assessed with PERMDISP
#>   (vegan::betadisper; Anderson 2006): p = 0.002. Figures were produced in R
#>   4.6.1 (R Core Team 2026) with micropubr 0.3.0 (Daniilidis 2026) and ggplot2
#>   4.0.3 (Wickham 2016), with phyloseq 1.56.0, rstatix 1.1.0 and vegan 2.7.6.
#> 
#> References
#> Anderson MJ (2001). A new method for non-parametric multivariate analysis of
#>   variance. Austral Ecology 26: 32-46.
#> Anderson MJ (2006). Distance-based tests for homogeneity of multivariate
#>   dispersions. Biometrics 62: 245-253.
#> Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
#>   practical and powerful approach to multiple testing. Journal of the Royal
#>   Statistical Society Series B 57: 289-300.
#> Daniilidis K (2026). micropubr: a unified grammar for publication-ready
#>   microbiome visualization in R. R package version 0.3.0.
#>   https://github.com/daniilidisK/micropubr
#> Dunn OJ (1964). Multiple comparisons using rank sums. Technometrics 6: 241-252.
#> Kassambara A (2026). rstatix: Pipe-Friendly Framework for Basic Statistical
#>   Tests. R package version 1.1.0.
#> Kruskal WH, Wallis WA (1952). Use of ranks in one-criterion variance analysis.
#>   Journal of the American Statistical Association 47: 583-621.
#> McMurdie PJ, Holmes S (2013). phyloseq: an R package for reproducible
#>   interactive analysis and graphics of microbiome census data. PLoS ONE 8:
#>   e61217.
#> R Core Team (2026). R: A Language and Environment for Statistical Computing. R
#>   Foundation for Statistical Computing, Vienna, Austria.
#> Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag,
#>   New York.
#> 
#> Statistical notes
#> * Group(s) Feces (n = 4), Freshwater (n = 2), Freshwater (creek) (n = 3), Mock
#>   (n = 3), Ocean (n = 3), Sediment (estuary) (n = 3), Skin (n = 3), Soil (n =
#>   3), Tongue (n = 2) have fewer than 5 samples; group comparisons have little
#>   power and their p-values are unstable at this size.
#> * Group dispersions differ (PERMDISP p = 0.002); a significant PERMANOVA may
#>   reflect differences in spread rather than in location (centroids).
# }
```
