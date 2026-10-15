# Package index

## Package overview

- [`micropubr`](https://daniilidisk.github.io/micropubr/reference/micropubr-package.md)
  [`micropubr-package`](https://daniilidisk.github.io/micropubr/reference/micropubr-package.md)
  : micropubr: Publication-Ready Microbiome Visualizations

## Import & data

Read DADA2, QIIME 2, BIOM and MetaPhlAn output (or convert other
containers) into a `phyloseq` object.

- [`mp_import_dada2()`](https://daniilidisk.github.io/micropubr/reference/mp_import_dada2.md)
  : Build a phyloseq object from DADA2 output
- [`mp_import_qiime2()`](https://daniilidisk.github.io/micropubr/reference/mp_import_qiime2.md)
  : Import QIIME 2 artifacts into a phyloseq object
- [`mp_import_biom()`](https://daniilidisk.github.io/micropubr/reference/mp_import_biom.md)
  : Import a BIOM file into a phyloseq object
- [`mp_import_metaphlan()`](https://daniilidisk.github.io/micropubr/reference/mp_import_metaphlan.md)
  : Import a merged MetaPhlAn abundance table
- [`mp_as_phyloseq()`](https://daniilidisk.github.io/micropubr/reference/mp_as_phyloseq.md)
  : Convert a (Tree)SummarizedExperiment to a phyloseq object

## Quality control

Per-cycle read-quality profiles from FASTQ files or summary tables.

- [`mp_qc_profile()`](https://daniilidisk.github.io/micropubr/reference/mp_qc_profile.md)
  : Per-cycle sequence quality profile

## Diversity

Alpha diversity, ordination with PERMANOVA, multivariate dispersion,
Mantel tests and variance partitioning.

- [`mp_alpha()`](https://daniilidisk.github.io/micropubr/reference/mp_alpha.md)
  : Plot alpha (within-sample) diversity
- [`mp_ordination()`](https://daniilidisk.github.io/micropubr/reference/mp_ordination.md)
  : Ordination plot with ellipses and PERMANOVA
- [`mp_dispersion()`](https://daniilidisk.github.io/micropubr/reference/mp_dispersion.md)
  : Multivariate dispersion (PERMDISP) with optional group tests
- [`mp_mantel()`](https://daniilidisk.github.io/micropubr/reference/mp_mantel.md)
  : Mantel test between community and explanatory distances
- [`mp_varpart()`](https://daniilidisk.github.io/micropubr/reference/mp_varpart.md)
  : Variance partitioning of community structure

## Composition

Taxonomic composition, prevalence / core microbiome and chord diagrams.

- [`mp_taxa()`](https://daniilidisk.github.io/micropubr/reference/mp_taxa.md)
  : Plot taxonomic composition
- [`mp_prevalence()`](https://daniilidisk.github.io/micropubr/reference/mp_prevalence.md)
  : Prevalence-abundance and core-microbiome plots
- [`mp_chord()`](https://daniilidisk.github.io/micropubr/reference/mp_chord.md)
  : Chord diagram of taxon-group relationships

## Differential abundance & biomarkers

Run and visualise differential-abundance analyses, indicator species and
machine-learning biomarker discovery.

- [`mp_run_da()`](https://daniilidisk.github.io/micropubr/reference/mp_run_da.md)
  : Run a differential-abundance analysis
- [`mp_diff()`](https://daniilidisk.github.io/micropubr/reference/mp_diff.md)
  : Plot differential-abundance results
- [`mp_indicator()`](https://daniilidisk.github.io/micropubr/reference/mp_indicator.md)
  : Indicator species analysis
- [`mp_ml()`](https://daniilidisk.github.io/micropubr/reference/mp_ml.md)
  : Machine-learning classification and regression on microbiome
  features

## Networks & phylogeny

- [`mp_network()`](https://daniilidisk.github.io/micropubr/reference/mp_network.md)
  : Plot a taxon co-occurrence network
- [`mp_tree()`](https://daniilidisk.github.io/micropubr/reference/mp_tree.md)
  : Plot an annotated phylogenetic tree

## Longitudinal & function

- [`mp_longitudinal()`](https://daniilidisk.github.io/micropubr/reference/mp_longitudinal.md)
  : Plot longitudinal and paired-sample microbiome data
- [`mp_function()`](https://daniilidisk.github.io/micropubr/reference/mp_function.md)
  : Plot functional pathway profiles

## Microbial genomics

Genome / contig coverage tracks and circular (Circos-style) genome
plots.

- [`mp_coverage()`](https://daniilidisk.github.io/micropubr/reference/mp_coverage.md)
  : Genome / contig coverage plot
- [`mp_circos()`](https://daniilidisk.github.io/micropubr/reference/mp_circos.md)
  : Circular (Circos-style) genome plot

## Interactive exploration

- [`mp_explore()`](https://daniilidisk.github.io/micropubr/reference/mp_explore.md)
  : Interactively explore a phyloseq object

## Reproducibility & reporting

Methods text, figure alt text, multi-panel assembly and journal-ready
export.

- [`mp_methods()`](https://daniilidisk.github.io/micropubr/reference/mp_methods.md)
  : Write a methods paragraph for micropubr figures
- [`mp_alt_text()`](https://daniilidisk.github.io/micropubr/reference/mp_alt_text.md)
  : Generate alt text for a micropubr figure
- [`mp_combine()`](https://daniilidisk.github.io/micropubr/reference/mp_combine.md)
  : Combine micropubr plots into a labelled multi-panel figure
- [`mp_export()`](https://daniilidisk.github.io/micropubr/reference/mp_export.md)
  : Export a figure at publication resolution and journal column widths

## Themes, palettes & options

Journal themes, colourblind-safe palettes and session-wide defaults
([`mp_set_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md),
[`mp_get_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md),
[`mp_reset_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md)).

- [`mp_theme_publication()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_publication.md)
  : Default publication theme for micropubr figures
- [`mp_theme_nature()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_nature.md)
  : Nature-style theme (compact, 8pt sans serif)
- [`mp_theme_cell()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_cell.md)
  : Cell-style theme (boxed panels, slightly larger type)
- [`mp_theme_dark()`](https://daniilidisk.github.io/micropubr/reference/mp_theme_dark.md)
  : Dark theme for slides and posters
- [`mp_darken()`](https://daniilidisk.github.io/micropubr/reference/mp_darken.md)
  : Convert a finished figure to dark mode
- [`mp_palette()`](https://daniilidisk.github.io/micropubr/reference/mp_palette.md)
  : micropubr colour palettes
- [`mp_palettes()`](https://daniilidisk.github.io/micropubr/reference/mp_palettes.md)
  : List available palettes
- [`mp_set_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md)
  [`mp_get_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md)
  [`mp_reset_defaults()`](https://daniilidisk.github.io/micropubr/reference/mp_set_defaults.md)
  : Get, set and reset micropubr defaults
