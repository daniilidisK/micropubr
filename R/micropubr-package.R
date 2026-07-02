#' micropubr: Publication-Ready Microbiome Visualizations
#'
#' @description
#' micropubr provides a unified, opinionated family of `mp_*()` functions that
#' turn a [phyloseq::phyloseq] object into a publication-ready
#' [ggplot2][ggplot2::ggplot] figure with a single call. It applies one
#' consistent design philosophy across the whole microbiome workflow: consistent
#' arguments, colourblind-safe defaults, automatic statistical overlays and
#' reproducible,
#' auto-generated figure captions, while always returning a plain `ggplot`
#' object you can extend with `+`.
#'
#' @section Core plotting functions:
#' \describe{
#'   \item{[mp_alpha()]}{Alpha (within-sample) diversity with statistics.}
#'   \item{[mp_ordination()]}{Beta diversity / ordination with PERMANOVA.}
#'   \item{[mp_taxa()]}{Taxonomic composition (bar, area, heatmap, bubble).}
#'   \item{[mp_diff()]}{Differential-abundance result visualisation.}
#'   \item{[mp_tree()]}{Annotated phylogenetic trees.}
#'   \item{[mp_network()]}{Co-occurrence / correlation networks.}
#'   \item{[mp_longitudinal()]}{Time-course and paired-sample plots.}
#'   \item{[mp_function()]}{Functional pathway profiles (PICRUSt2/HUMAnN).}
#' }
#'
#' @section Statistical analyses:
#' \describe{
#'   \item{[mp_dispersion()]}{Multivariate dispersion (PERMDISP) with optional
#'     ANOSIM / MRPP / PERMANOVA.}
#'   \item{[mp_mantel()]}{Mantel test and correlogram.}
#'   \item{[mp_varpart()]}{Variance partitioning among explanatory sets.}
#'   \item{[mp_prevalence()]}{Prevalence-abundance and core-microbiome plots.}
#'   \item{[mp_indicator()]}{Indicator species analysis.}
#'   \item{[mp_run_da()]}{Run a differential-abundance method for [mp_diff()].}
#' }
#'
#' @section Microbial genomics:
#' \describe{
#'   \item{[mp_coverage()]}{Genome / contig sequencing-depth tracks.}
#'   \item{[mp_circos()]}{Circular (Circos-style) genome plots.}
#' }
#'
#' @section Infrastructure:
#' [mp_palette()], the `mp_theme_*()` family, [mp_set_defaults()],
#' [mp_combine()], [mp_export()] and the `mp_import_*()` helpers. Every
#' plotting function accepts `return.data = TRUE` to return the tidy data
#' frame behind the plot, and `...` to forward arguments to its primary geom.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import ggplot2
#' @importFrom rlang .data %||% abort warn inform arg_match
#' @importFrom stats as.formula as.dist reorder sd setNames complete.cases
#' @importFrom methods as is
## usethis namespace: end
NULL

# Quiet R CMD check notes for bare symbols passed to ggtreeExtra::geom_fruit(),
# which requires un-namespaced aesthetics (it calls rlang::as_name() on them).
utils::globalVariables(c(".", "label", "value", "group_value"))
