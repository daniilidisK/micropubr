pkgname <- "micropubr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('micropubr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("mp_alpha")
### * mp_alpha

flush(stderr()); flush(stdout())

### Name: mp_alpha
### Title: Plot alpha (within-sample) diversity
### Aliases: mp_alpha

### ** Examples

## Not run: 
##D data(GlobalPatterns, package = "phyloseq")
##D mp_alpha(GlobalPatterns, measure = c("Shannon", "Observed"),
##D          group.by = "SampleType")
## End(Not run)



cleanEx()
nameEx("mp_circos")
### * mp_circos

flush(stderr()); flush(stdout())

### Name: mp_circos
### Title: Circular (Circos-style) genome plot
### Aliases: mp_circos

### ** Examples

## Not run: 
##D genome <- data.frame(contig = c("c1", "c2"), length = c(5000, 3000))
##D cov <- data.frame(contig = "c1", position = seq(1, 5000, 50),
##D                   value = runif(100, 10, 100))
##D mp_circos(genome, coverage = cov)
## End(Not run)



cleanEx()
nameEx("mp_combine")
### * mp_combine

flush(stderr()); flush(stdout())

### Name: mp_combine
### Title: Combine micropubr plots into a labelled multi-panel figure
### Aliases: mp_combine

### ** Examples

## Not run: 
##D a <- mp_alpha(ps, group.by = "Group")
##D o <- mp_ordination(ps, group.by = "Group")
##D mp_combine(a, o, ncol = 2, tag.levels = "A")
## End(Not run)



cleanEx()
nameEx("mp_coverage")
### * mp_coverage

flush(stderr()); flush(stdout())

### Name: mp_coverage
### Title: Genome / contig coverage plot
### Aliases: mp_coverage

### ** Examples

cov <- data.frame(
  contig = "chr1", position = 1:500,
  depth = 50 + as.integer(20 * sin(seq(0, 6 * pi, length.out = 500)))
)
mp_coverage(cov, plot.type = "area")



cleanEx()
nameEx("mp_diff")
### * mp_diff

flush(stderr()); flush(stdout())

### Name: mp_diff
### Title: Plot differential-abundance results
### Aliases: mp_diff

### ** Examples

df <- data.frame(
  feature = paste0("ASV", 1:50),
  log2FC = rnorm(50, sd = 2),
  pvalue = runif(50)^2,
  padj = p.adjust(runif(50)^2, "BH")
)
mp_diff(df, plot.type = "volcano")



cleanEx()
nameEx("mp_dispersion")
### * mp_dispersion

flush(stderr()); flush(stdout())

### Name: mp_dispersion
### Title: Multivariate dispersion (PERMDISP) with optional group tests
### Aliases: mp_dispersion

### ** Examples

## Not run: 
##D data(GlobalPatterns, package = "phyloseq")
##D mp_dispersion(GlobalPatterns, group.by = "SampleType",
##D               add.test = "permanova")
## End(Not run)



cleanEx()
nameEx("mp_explore")
### * mp_explore

flush(stderr()); flush(stdout())

### Name: mp_explore
### Title: Interactively explore a phyloseq object
### Aliases: mp_explore

### ** Examples

## Not run: 
##D data(GlobalPatterns, package = "phyloseq")
##D p <- mp_explore(GlobalPatterns)
## End(Not run)



cleanEx()
nameEx("mp_export")
### * mp_export

flush(stderr()); flush(stdout())

### Name: mp_export
### Title: Export a figure at publication resolution and journal column
###   widths
### Aliases: mp_export

### ** Examples

## Not run: 
##D p <- mp_ordination(ps, group.by = "Group")
##D mp_export(p, "figure1.pdf", column = "single")
##D mp_export(p, "figure1.tiff", column = "double", dpi = 600)
## End(Not run)



cleanEx()
nameEx("mp_function")
### * mp_function

flush(stderr()); flush(stdout())

### Name: mp_function
### Title: Plot functional pathway profiles
### Aliases: mp_function

### ** Examples

set.seed(1)
m <- matrix(rpois(20 * 12, 50), nrow = 20,
            dimnames = list(paste0("path", 1:20), paste0("S", 1:12)))
meta <- data.frame(Group = rep(c("A", "B"), each = 6),
                   row.names = colnames(m))
mp_function(m, metadata = meta, plot.type = "heatmap", group.by = "Group")



cleanEx()
nameEx("mp_longitudinal")
### * mp_longitudinal

flush(stderr()); flush(stdout())

### Name: mp_longitudinal
### Title: Plot longitudinal and paired-sample microbiome data
### Aliases: mp_longitudinal

### ** Examples

## Not run: 
##D mp_longitudinal(ps, time.var = "Week", subject.var = "SubjectID",
##D                 group.by = "Arm", measure = "Shannon")
## End(Not run)



cleanEx()
nameEx("mp_network")
### * mp_network

flush(stderr()); flush(stdout())

### Name: mp_network
### Title: Plot a taxon co-occurrence network
### Aliases: mp_network

### ** Examples

## Not run: 
##D data(GlobalPatterns, package = "phyloseq")
##D gp <- phyloseq::tax_glom(GlobalPatterns, "Genus")
##D mp_network(gp, top.n = 60, node.color.by = "Phylum")
## End(Not run)



cleanEx()
nameEx("mp_ordination")
### * mp_ordination

flush(stderr()); flush(stdout())

### Name: mp_ordination
### Title: Ordination plot with ellipses and PERMANOVA
### Aliases: mp_ordination

### ** Examples

## Not run: 
##D data(GlobalPatterns, package = "phyloseq")
##D mp_ordination(GlobalPatterns, method = "PCoA", distance = "bray",
##D               group.by = "SampleType")
## End(Not run)



cleanEx()
nameEx("mp_palette")
### * mp_palette

flush(stderr()); flush(stdout())

### Name: mp_palette
### Title: micropubr colour palettes
### Aliases: mp_palette

### ** Examples

mp_palette("taxa10")
mp_palette("diverging", n = 7)



cleanEx()
nameEx("mp_run_da")
### * mp_run_da

flush(stderr()); flush(stdout())

### Name: mp_run_da
### Title: Run a differential-abundance analysis from a phyloseq object
### Aliases: mp_run_da

### ** Examples

## Not run: 
##D data(GlobalPatterns, package = "phyloseq")
##D res <- mp_run_da(GlobalPatterns, group.by = "SampleType",
##D                  groups = c("Feces", "Soil"), method = "wilcox")
##D mp_diff(res, plot.type = "volcano")
## End(Not run)



cleanEx()
nameEx("mp_set_defaults")
### * mp_set_defaults

flush(stderr()); flush(stdout())

### Name: mp_set_defaults
### Title: Get, set and reset micropubr defaults
### Aliases: mp_set_defaults mp_get_defaults mp_reset_defaults

### ** Examples

old <- mp_get_defaults()
mp_set_defaults(font.size = 9, palette = "okabe")
mp_get_defaults()
mp_reset_defaults()



cleanEx()
nameEx("mp_taxa")
### * mp_taxa

flush(stderr()); flush(stdout())

### Name: mp_taxa
### Title: Plot taxonomic composition
### Aliases: mp_taxa

### ** Examples

## Not run: 
##D data(GlobalPatterns, package = "phyloseq")
##D mp_taxa(GlobalPatterns, tax.level = "Phylum", group.by = "SampleType")
## End(Not run)



cleanEx()
nameEx("mp_theme_publication")
### * mp_theme_publication

flush(stderr()); flush(stdout())

### Name: mp_theme_publication
### Title: Default publication theme for micropubr figures
### Aliases: mp_theme_publication

### ** Examples

library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) +
  geom_point() +
  mp_theme_publication()



cleanEx()
nameEx("mp_tree")
### * mp_tree

flush(stderr()); flush(stdout())

### Name: mp_tree
### Title: Plot an annotated phylogenetic tree
### Aliases: mp_tree

### ** Examples

## Not run: 
##D data(GlobalPatterns, package = "phyloseq")
##D gp <- phyloseq::tax_glom(GlobalPatterns, "Phylum")
##D mp_tree(gp, layout = "circular", color.by = "Phylum",
##D         group.by = "SampleType")
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
