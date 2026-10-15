# Import QIIME 2 artifacts into a phyloseq object

Import QIIME 2 artifacts into a phyloseq object

## Usage

``` r
mp_import_qiime2(features, taxonomy = NULL, tree = NULL, metadata = NULL)
```

## Arguments

- features:

  Path to a feature-table `.qza` (`FeatureTable[Frequency]`).

- taxonomy:

  Optional path to a taxonomy `.qza`.

- tree:

  Optional path to a rooted-tree `.qza`.

- metadata:

  Optional path to a tab-delimited sample-metadata file, or a data
  frame.

## Value

A [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
object.

## See also

[`mp_import_dada2()`](https://daniilidisk.github.io/micropubr/reference/mp_import_dada2.md),
[`mp_import_biom()`](https://daniilidisk.github.io/micropubr/reference/mp_import_biom.md),
[`mp_import_metaphlan()`](https://daniilidisk.github.io/micropubr/reference/mp_import_metaphlan.md)
