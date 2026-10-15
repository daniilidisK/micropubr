# Generate alt text for a micropubr figure

Describes a figure for readers who cannot see it: the chart type, what
is plotted against what, how it is grouped, how much data it shows, the
key result (e.g. which group has the highest diversity, how many taxa
are differentially abundant) and the statistics in its subtitle and
caption. Journals increasingly ask for alt text with every figure.

## Usage

``` r
mp_alt_text(plot)
```

## Arguments

- plot:

  A figure made by an `mp_*()` function (any `ggplot` works, with a more
  generic description).

## Value

A character string.

## Examples

``` r
# \donttest{
data(GlobalPatterns, package = "phyloseq")
p <- mp_alpha(GlobalPatterns, measure = "Shannon", group.by = "SampleType",
              test = "kruskal")
#> Warning: Group(s) Feces (n = 4), Freshwater (n = 2), Freshwater (creek) (n = 3), Mock (n = 3), Ocean (n = 3), Sediment (estuary) (n = 3), Skin (n = 3), Soil (n = 3), Tongue (n = 2) have fewer than 5 samples; group comparisons have little power and their p-values are unstable at this size.
mp_alt_text(p)
#> [1] "Box plots of Shannon, coloured by SampleType (9 categories), based on 26 data rows. Median Shannon is highest in Soil (6.58) and lowest in Tongue (3.29). Statistics: Kruskal-Wallis, p = 0.005."
# }
```
