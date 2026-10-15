# Interactively explore a phyloseq object

Launches a Shiny gadget that wraps
[`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md),
[`mp_ordination()`](https://daniilidisk.github.io/micropubr/reference/mp_ordination.md)
and
[`mp_taxa()`](https://daniilidisk.github.io/micropubr/reference/mp_taxa.md)
behind interactive controls, with a live preview and the equivalent
reproducible `micropubr` code. Closing the gadget with **Done** returns
the currently displayed `ggplot` object so it can be saved or further
edited.

## Usage

``` r
mp_explore(object)
```

## Arguments

- object:

  A [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
  object.

## Value

Invisibly, the last `ggplot` shown (when run interactively).

## Examples

``` r
if (FALSE) { # \dontrun{
data(GlobalPatterns, package = "phyloseq")
p <- mp_explore(GlobalPatterns)
} # }
```
