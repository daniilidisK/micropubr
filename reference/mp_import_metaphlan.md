# Import a merged MetaPhlAn abundance table

Parses the pipe-delimited lineage column of a merged MetaPhlAn3/4 table
and builds a phyloseq object at the requested rank.

## Usage

``` r
mp_import_metaphlan(path, metadata = NULL, level = "Species")
```

## Arguments

- path:

  Path to a merged MetaPhlAn table (tab-delimited; first column the
  clade lineage).

- metadata:

  Optional sample metadata (path or data frame).

- level:

  Taxonomic level to keep: one of `"Kingdom"`, `"Phylum"`, `"Class"`,
  `"Order"`, `"Family"`, `"Genus"`, `"Species"` (default).

## Value

A [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
object.
