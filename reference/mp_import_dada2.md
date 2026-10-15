# Build a phyloseq object from DADA2 output

Build a phyloseq object from DADA2 output

## Usage

``` r
mp_import_dada2(
  seqtab,
  taxa = NULL,
  metadata = NULL,
  tree = NULL,
  short.names = TRUE
)
```

## Arguments

- seqtab:

  A DADA2 sequence table (samples x ASVs), e.g. from
  `dada2::makeSequenceTable()`.

- taxa:

  Optional taxonomy matrix (ASVs x ranks) from
  `dada2::assignTaxonomy()`.

- metadata:

  Optional sample metadata data frame (rows = samples).

- tree:

  Optional `phylo` object.

- short.names:

  Replace long sequence identifiers with `ASV1`, `ASV2`, ... and store
  the original sequences as `refseq` (requires `Biostrings`).

## Value

A [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
object.
