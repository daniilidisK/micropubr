# Convert a (Tree)SummarizedExperiment to a phyloseq object

Bioconductor's microbiome stack (`mia`, `miaViz`, the OMA book) stores
data in a
[TreeSummarizedExperiment](https://rdrr.io/pkg/TreeSummarizedExperiment/man/TreeSummarizedExperiment-constructor.html)
rather than a `phyloseq` object. Every `mp_*()` function calls this
conversion on its input, so either container can be passed directly; it
is exported for when you want the `phyloseq` object itself.

## Usage

``` r
mp_as_phyloseq(x, assay.type = NULL)
```

## Arguments

- x:

  A `TreeSummarizedExperiment` or `SummarizedExperiment`; a `phyloseq`
  object is returned unchanged.

- assay.type:

  Name of the assay holding the counts. Defaults to `"counts"` when
  present, otherwise the first assay.

## Value

A [phyloseq::phyloseq](https://rdrr.io/pkg/phyloseq/man/phyloseq.html)
object.

## Details

The count assay becomes the OTU table, `colData()` the sample data, the
taxonomic-rank columns of `rowData()` the taxonomy table and, for a
`TreeSummarizedExperiment`, `rowTree()` the phylogeny (pruned and
relabelled to the rows it links to). Non-atomic `colData()` columns,
which `phyloseq` cannot hold, are dropped.

## Examples

``` r
if (requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
  counts <- matrix(rpois(40, 20), nrow = 5,
                   dimnames = list(paste0("OTU", 1:5), paste0("S", 1:8)))
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = counts),
    rowData = S4Vectors::DataFrame(Phylum = rep(c("A", "B"), c(2, 3))),
    colData = S4Vectors::DataFrame(Group = rep(c("x", "y"), 4)))
  mp_as_phyloseq(tse)
}
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 5 taxa and 8 samples ]
#> sample_data() Sample Data:       [ 8 samples by 1 sample variables ]
#> tax_table()   Taxonomy Table:    [ 5 taxa by 1 taxonomic ranks ]
```
