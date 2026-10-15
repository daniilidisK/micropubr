# Chord diagram of taxon-group relationships -----------------------------
# A circlize chord diagram linking taxa to sample groups by their (relative)
# abundance: a compact way to show which groups are dominated by which taxa.

#' Chord diagram of taxon-group relationships
#'
#' Draws a circular chord diagram with `circlize` connecting the top taxa at a
#' chosen rank to the levels of a grouping variable, with ribbon widths
#' proportional to mean (relative) abundance. Like [mp_circos()] this draws to
#' the current graphics device (the exception to the ggplot contract) and
#' returns the underlying taxon-by-group matrix invisibly.
#'
#' @param object A [phyloseq::phyloseq] object.
#' @param group.by Sample-metadata column defining the groups (columns of the
#'   diagram).
#' @param tax.level Taxonomic rank to agglomerate to (rows of the diagram).
#' @param top.n Number of most-abundant taxa to show; the rest collapse into
#'   an "Other" sector.
#' @param relative Use relative abundance (`TRUE`) or raw counts.
#' @param palette Palette name or colour vector for the taxa sectors/ribbons.
#' @param gap Gap (degrees) between sectors.
#' @param start.degree Starting angle.
#' @param transparency Ribbon transparency in `[0, 1]`.
#' @param directional Passed to [circlize::chordDiagram()]: `0` (undirected),
#'   `1` (taxa to groups) or `-1`.
#' @param title Optional plot title.
#' @param verbose Print the methods used.
#'
#' @return Invisibly, the taxon-by-group abundance matrix that was drawn.
#' @export
#' @examples
#' \dontrun{
#' data(GlobalPatterns, package = "phyloseq")
#' mp_chord(GlobalPatterns, group.by = "SampleType", tax.level = "Phylum")
#' }
mp_chord <- function(object,
                     group.by,
                     tax.level = "Phylum",
                     top.n = 12,
                     relative = TRUE,
                     palette = NULL,
                     gap = 2,
                     start.degree = 90,
                     transparency = 0.5,
                     directional = 0,
                     title = NULL,
                     verbose = getOption("micropubr.verbose", FALSE)) {
  object <- .mp_check_phyloseq(object)
  group.by <- .mp_check_var(object, group.by, "group.by", required = TRUE)
  .mp_require("circlize", "Chord diagrams")

  df <- .mp_melt_taxa(object, tax.level = tax.level, top.n = top.n,
                      relative = relative)
  if (!group.by %in% names(df)) {
    rlang::abort(sprintf("`group.by = \"%s\"` is not available after melting.",
                         group.by),
                 class = "micropubr_bad_var")
  }
  agg <- stats::aggregate(
    stats::as.formula(sprintf("Abundance ~ .taxon + `%s`", group.by)),
    data = df, FUN = mean, na.rm = TRUE)
  mat <- stats::xtabs(
    stats::as.formula(sprintf("Abundance ~ .taxon + `%s`", group.by)),
    data = agg)
  mat <- as.matrix(unclass(mat))

  if (verbose) {
    rlang::inform(sprintf(
      "mp_chord(): %d taxa x %d groups at rank '%s'.",
      nrow(mat), ncol(mat), tax.level))
  }

  taxa_cols <- .mp_taxa_colors(palette, df$.taxon, "Other")[rownames(mat)]
  group_cols <- stats::setNames(
    grDevices::grey.colors(ncol(mat), start = 0.35, end = 0.7), colnames(mat))
  grid_col <- c(taxa_cols, group_cols)

  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  circlize::circos.par(start.degree = start.degree, gap.degree = gap)
  circlize::chordDiagram(
    mat, grid.col = grid_col, transparency = transparency,
    directional = directional, annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.12))
  # Sector labels rotated to read radially.
  circlize::circos.trackPlotRegion(
    track.index = 1, bg.border = NA,
    panel.fun = function(x, y) {
      si <- circlize::CELL_META$sector.index
      xlim <- circlize::CELL_META$xlim
      ylim <- circlize::CELL_META$ylim
      circlize::circos.text(mean(xlim), ylim[1] + 0.1, si,
                            facing = "clockwise", niceFacing = TRUE,
                            adj = c(0, 0.5), cex = 0.6)
    })
  if (!is.null(title)) graphics::title(main = title)

  invisible(mat)
}
