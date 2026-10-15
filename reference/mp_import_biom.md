# Import a BIOM file into a phyloseq object

Import a BIOM file into a phyloseq object

## Usage

``` r
mp_import_biom(
  path,
  metadata = NULL,
  parse.taxonomy = c("default", "greengenes"),
  ...
)
```

## Arguments

- path:

  Path to a `.biom` file.

- metadata:

  Optional sample metadata (path or data frame) to merge in.

- parse.taxonomy:

  Taxonomy parsing style: `"default"` or `"greengenes"`.

- ...:

  Passed to
  [`phyloseq::import_biom()`](https://rdrr.io/pkg/phyloseq/man/import_biom.html).

## Value

A [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
object.
