# Microbial genomics visualisations --------------------------------------
# Sequencing-oriented plots that complement the community-ecology functions:
# per-base/windowed coverage tracks and circular (Circos-style) genome plots.
# Coverage accepts tidy data frames and common text formats (samtools depth,
# bedGraph, BED); BAM/CRAM input is supported optionally via Rsamtools.

#' Genome / contig coverage plot
#'
#' Draws sequencing-depth (coverage) tracks along a genome or contig from a
#' tidy data frame, a `samtools depth` / bedGraph / BED file, or -- optionally
#' -- a BAM/CRAM alignment. Multiple samples are overlaid by colour and
#' multiple contigs are faceted. Gene annotations and highlight regions can be
#' overlaid.
#'
#' @param x A data frame (with position/depth columns), a path to a
#'   `samtools depth`, bedGraph or BED file, or a path to a BAM/CRAM file.
#' @param format Input format: `"auto"` (guess from extension/columns),
#'   `"depth"`, `"bedgraph"`, `"bed"`, `"bam"` or `"data.frame"`.
#' @param region Optional region to restrict a BAM import, as
#'   `"contig:start-end"` or `"contig"`.
#' @param window Optional window size (bp) for binning to mean depth -- useful
#'   for whole-genome views.
#' @param contig Optional contig name(s) to keep.
#' @param sample Optional sample label (for single-sample inputs).
#' @param plot.type `"line"` or `"area"`.
#' @param log.y Plot depth on a log10 axis.
#' @param annotations Optional gene track: a data frame with `start`, `end`
#'   and optionally `strand` and `label` columns (and `contig` for multi-contig
#'   views).
#' @param highlight Optional data frame of regions to shade, with `start` and
#'   `end` (and optional `contig`).
#' @param palette,theme,plot.title,plot.subtitle,plot.caption,font.size See
#'   [mp_alpha()].
#' @param ... Passed to the primary geom.
#' @param return.data If `TRUE`, return the tidy coverage data frame.
#' @param verbose Print the methods used.
#'
#' @return A `ggplot` object, or a data frame when `return.data = TRUE`.
#' @export
#' @examples
#' cov <- data.frame(
#'   contig = "chr1", position = 1:500,
#'   depth = 50 + as.integer(20 * sin(seq(0, 6 * pi, length.out = 500)))
#' )
#' mp_coverage(cov, plot.type = "area")
mp_coverage <- function(x,
                        format = c("auto", "depth", "bedgraph", "bed",
                                   "bam", "data.frame"),
                        region = NULL,
                        window = NULL,
                        contig = NULL,
                        sample = NULL,
                        plot.type = c("line", "area"),
                        log.y = FALSE,
                        annotations = NULL,
                        highlight = NULL,
                        palette = NULL,
                        theme = NULL,
                        ...,
                        plot.title = NULL,
                        plot.subtitle = NULL,
                        plot.caption = NULL,
                        font.size = getOption("micropubr.font.size", 11),
                        return.data = FALSE,
                        verbose = getOption("micropubr.verbose", FALSE)) {
  format <- match.arg(format)
  plot.type <- match.arg(plot.type)
  theme <- .mp_theme(theme, font.size)

  df <- .mp_read_coverage(x, format, region, sample)
  if (!is.null(contig)) df <- df[df$contig %in% contig, , drop = FALSE]
  if (nrow(df) == 0) {
    rlang::abort("No coverage rows after filtering.",
                 class = "micropubr_empty_coverage")
  }
  if (!is.null(window)) df <- .mp_window_coverage(df, window)
  df <- df[order(df$contig, df$sample, df$position), , drop = FALSE]

  if (verbose) {
    rlang::inform(sprintf(
      "mp_coverage(): %d positions, %d contig(s), %d sample(s)%s.",
      nrow(df), length(unique(df$contig)), length(unique(df$sample)),
      if (!is.null(window)) sprintf(", %d bp windows", window) else ""))
  }
  if (isTRUE(return.data)) return(df)

  nsamp <- length(unique(df$sample))
  cols <- .mp_resolve_palette(palette, max(nsamp, 1), type = "qualitative")

  mapping <- ggplot2::aes(x = .data$position, y = .data$depth,
                          colour = .data$sample, fill = .data$sample)
  p <- ggplot2::ggplot(df, mapping)

  if (!is.null(highlight)) {
    p <- p + .mp_coverage_highlight(highlight)
  }
  if (plot.type == "area") {
    p <- p + .mp_geom(ggplot2::geom_area,
                      params = list(position = "identity", alpha = 0.4,
                                    linewidth = 0.4),
                      dots = list(...))
  } else {
    p <- p + .mp_geom(ggplot2::geom_line,
                      params = list(linewidth = 0.5), dots = list(...))
  }

  p <- p +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::labs(x = "Position (bp)", y = "Depth",
                  title = plot.title, subtitle = plot.subtitle,
                  caption = plot.caption) +
    theme

  if (log.y) p <- p + ggplot2::scale_y_log10()
  if (length(unique(df$contig)) > 1) {
    p <- p + ggplot2::facet_wrap(~ contig, scales = "free_x")
  }
  if (nsamp < 2) p <- p + ggplot2::theme(legend.position = "none")
  if (!is.null(annotations)) {
    p <- .mp_coverage_annotations(p, df, annotations)
  }
  p
}

#' Dispatch coverage input to a tidy data frame
#'
#' Returns a data frame with columns `contig`, `position`, `depth`, `sample`.
#' @noRd
.mp_read_coverage <- function(x, format, region, sample) {
  if (is.data.frame(x)) return(.mp_coverage_from_df(x, sample))
  if (!is.character(x) || length(x) != 1) {
    rlang::abort(
      "`x` must be a data frame or a single file path.",
      class = "micropubr_bad_coverage_input")
  }
  if (format == "auto") format <- .mp_guess_cov_format(x)
  switch(
    format,
    bam      = .mp_coverage_from_bam(x, region, sample),
    bedgraph = .mp_coverage_from_table(x, "bedgraph", sample),
    bed      = .mp_coverage_from_table(x, "bed", sample),
    depth    = .mp_coverage_from_table(x, "depth", sample),
    rlang::abort(sprintf("Cannot read coverage format '%s'.", format),
                 class = "micropubr_bad_coverage_format"))
}

#' Guess a coverage file format from its extension
#' @noRd
.mp_guess_cov_format <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("bam", "cram")) return("bam")
  if (ext %in% c("bedgraph", "bg")) return("bedgraph")
  if (ext == "bed") return("bed")
  "depth"
}

#' Standardise a user data frame to the coverage schema
#' @noRd
.mp_coverage_from_df <- function(df, sample) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  contig <- .mp_pick_col(df, c("contig", "chrom", "chr", "seqnames",
                               "rname"), default = "contig1")
  pos <- .mp_pick_col(df, c("position", "pos", "start", "loc"))
  depth <- .mp_pick_col(df, c("depth", "coverage", "cov", "score", "value"))
  if (is.null(pos) || is.null(depth)) {
    rlang::abort(
      c("Could not find position and depth columns in the data frame.",
        i = "Expected columns like 'position'/'pos' and 'depth'/'coverage'."),
      class = "micropubr_bad_coverage_cols")
  }
  smp <- .mp_pick_col(df, c("sample", "Sample", "library"),
                      default = sample %||% "sample1")
  data.frame(contig = as.character(contig), position = as.numeric(pos),
             depth = as.numeric(depth), sample = as.character(smp),
             stringsAsFactors = FALSE)
}

#' Read a samtools depth / bedGraph / BED file
#' @noRd
.mp_coverage_from_table <- function(path, type, sample) {
  raw <- utils::read.table(path, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, comment.char = "#")
  smp <- sample %||% tools::file_path_sans_ext(basename(path))
  if (type == "depth") {
    out <- data.frame(contig = as.character(raw[[1]]),
                      position = as.numeric(raw[[2]]),
                      depth = as.numeric(raw[[3]]))
  } else {
    val <- if (ncol(raw) >= 4) as.numeric(raw[[4]]) else 1
    out <- data.frame(contig = as.character(raw[[1]]),
                      position = as.numeric(raw[[2]]) + 1,
                      depth = val)
  }
  out$sample <- smp
  out
}

#' Per-base coverage from a BAM/CRAM file (optional Rsamtools backend)
#' @noRd
.mp_coverage_from_bam <- function(path, region, sample) {
  .mp_require(c("Rsamtools", "GenomicAlignments", "GenomicRanges",
                "IRanges"), "BAM/CRAM coverage", bioc = TRUE)
  param <- Rsamtools::ScanBamParam()
  if (!is.null(region)) {
    reg <- .mp_parse_region(region)
    gr <- GenomicRanges::GRanges(
      reg$contig, IRanges::IRanges(reg$start %||% 1, reg$end %||% .Machine$integer.max))
    param <- Rsamtools::ScanBamParam(which = gr)
  }
  ga <- GenomicAlignments::readGAlignments(path, param = param)
  cov <- GenomicAlignments::coverage(ga)
  smp <- sample %||% tools::file_path_sans_ext(basename(path))
  parts <- lapply(names(cov), function(sn) {
    v <- as.numeric(cov[[sn]])
    if (length(v) == 0) return(NULL)
    data.frame(contig = sn, position = seq_along(v), depth = v,
               sample = smp, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  if (!is.null(region)) {
    reg <- .mp_parse_region(region)
    if (!is.null(reg$start)) out <- out[out$position >= reg$start, ]
    if (!is.null(reg$end)) out <- out[out$position <= reg$end, ]
  }
  out
}

#' Parse a "contig:start-end" region string
#' @noRd
.mp_parse_region <- function(region) {
  if (!grepl(":", region)) return(list(contig = region, start = NULL,
                                       end = NULL))
  contig <- sub(":.*$", "", region)
  rng <- sub("^.*:", "", region)
  parts <- strsplit(gsub(",", "", rng), "-")[[1]]
  list(contig = contig,
       start = suppressWarnings(as.numeric(parts[1])),
       end = if (length(parts) > 1) suppressWarnings(as.numeric(parts[2]))
             else NULL)
}

#' Bin coverage into fixed windows (mean depth)
#' @noRd
.mp_window_coverage <- function(df, window) {
  df$position <- (df$position %/% window) * window + window / 2
  agg <- stats::aggregate(depth ~ contig + sample + position, data = df,
                          FUN = mean)
  agg
}

#' Shaded highlight regions for a coverage plot
#' @noRd
.mp_coverage_highlight <- function(highlight) {
  hl <- as.data.frame(highlight, stringsAsFactors = FALSE)
  ggplot2::geom_rect(
    data = hl,
    ggplot2::aes(xmin = .data$start, xmax = .data$end,
                 ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#FDE725", alpha = 0.25)
}

#' Overlay a gene-annotation track near the top of a coverage plot
#' @noRd
.mp_coverage_annotations <- function(p, df, annotations) {
  ann <- as.data.frame(annotations, stringsAsFactors = FALSE)
  ymax <- max(df$depth, na.rm = TRUE)
  ann$.ymin <- ymax * 1.04
  ann$.ymax <- ymax * 1.10
  ann$.ytext <- ymax * 1.12
  ann$.mid <- (ann$start + ann$end) / 2
  fillv <- if ("strand" %in% names(ann)) {
    ggplot2::aes(xmin = .data$start, xmax = .data$end,
                 ymin = .data$.ymin, ymax = .data$.ymax,
                 fill = .data$strand)
  } else {
    ggplot2::aes(xmin = .data$start, xmax = .data$end,
                 ymin = .data$.ymin, ymax = .data$.ymax)
  }
  p <- p + ggplot2::geom_rect(data = ann, mapping = fillv,
                              inherit.aes = FALSE, colour = "grey20",
                              alpha = 0.9)
  if ("label" %in% names(ann)) {
    p <- p + ggplot2::geom_text(
      data = ann,
      ggplot2::aes(x = .data$.mid, y = .data$.ytext,
                   label = .data$label),
      inherit.aes = FALSE, size = 2.5, vjust = 0)
  }
  p
}

# --- Circos-style circular genome plot ----------------------------------

#' Circular (Circos-style) genome plot
#'
#' Draws a circular genome plot with `circlize`: an ideogram of contigs plus
#' optional coverage, GC and annotation tracks and ribbon links (e.g. synteny
#' or correlations). Unlike the other `mp_*()` functions this draws with base
#' graphics to the current device (the one exception to the ggplot contract)
#' and returns the assembled track data invisibly.
#'
#' @param genome A data frame describing contigs with columns `contig` and
#'   `length`, or a named numeric vector of lengths.
#' @param coverage Optional coverage data frame with `contig`, `position` (or
#'   `start`) and `value`.
#' @param gc Optional GC-content/skew data frame with `contig`, `start`, `end`
#'   and `value`.
#' @param annotations Optional features with `contig`, `start`, `end` and
#'   optional `label`.
#' @param links Optional ribbon links with `contig1`, `start1`, `end1`,
#'   `contig2`, `start2`, `end2`.
#' @param coverage.type `"line"`, `"area"` or `"heatmap"`.
#' @param palette Palette name or colour vector for contigs.
#' @param gap Gap (degrees) between contigs.
#' @param start.degree Starting angle.
#' @param title Optional plot title.
#' @param verbose Print the methods used.
#'
#' @return Invisibly, a list of the track data used (so `return.data`-style
#'   workflows still work). Called for the side effect of drawing the plot.
#' @export
#' @examples
#' \dontrun{
#' genome <- data.frame(contig = c("c1", "c2"), length = c(5000, 3000))
#' cov <- data.frame(contig = "c1", position = seq(1, 5000, 50),
#'                   value = runif(100, 10, 100))
#' mp_circos(genome, coverage = cov)
#' }
mp_circos <- function(genome,
                      coverage = NULL,
                      gc = NULL,
                      annotations = NULL,
                      links = NULL,
                      coverage.type = c("line", "area", "heatmap"),
                      palette = NULL,
                      gap = 2,
                      start.degree = 90,
                      title = NULL,
                      verbose = getOption("micropubr.verbose", FALSE)) {
  .mp_require("circlize", "Circos genome plots")
  coverage.type <- match.arg(coverage.type)
  gdf <- .mp_circos_genome(genome)
  cols <- .mp_resolve_palette(palette, nrow(gdf), type = "qualitative")
  names(cols) <- gdf$contig

  if (verbose) {
    rlang::inform(sprintf("mp_circos(): %d contig(s)%s%s%s.",
                          nrow(gdf),
                          if (!is.null(coverage)) " + coverage" else "",
                          if (!is.null(gc)) " + GC" else "",
                          if (!is.null(links)) " + links" else ""))
  }

  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  circlize::circos.par(start.degree = start.degree, gap.degree = gap,
                       cell.padding = c(0, 0, 0, 0))
  circlize::circos.initialize(factors = gdf$contig,
                              xlim = cbind(0, gdf$length))

  # Ideogram track: coloured band + contig labels.
  circlize::circos.track(
    ylim = c(0, 1), track.height = 0.06, bg.border = NA,
    panel.fun = function(x, y) {
      si <- circlize::CELL_META$sector.index
      circlize::circos.rect(circlize::CELL_META$xlim[1], 0,
                            circlize::CELL_META$xlim[2], 1,
                            col = cols[[si]], border = NA)
      circlize::circos.text(circlize::CELL_META$xcenter, 1.8, si,
                            cex = 0.6, facing = "bending.inside",
                            niceFacing = TRUE)
    })

  if (!is.null(coverage)) {
    .mp_circos_coverage(coverage, coverage.type)
  }
  if (!is.null(gc)) {
    .mp_circos_gc(gc)
  }
  if (!is.null(annotations)) {
    .mp_circos_annotations(annotations)
  }
  if (!is.null(links)) {
    .mp_circos_links(links)
  }
  if (!is.null(title)) graphics::title(main = title)

  invisible(list(genome = gdf, coverage = coverage, gc = gc,
                 annotations = annotations, links = links))
}

#' Normalise the genome argument to a data frame of contig/length
#' @noRd
.mp_circos_genome <- function(genome) {
  if (is.numeric(genome) && !is.null(names(genome))) {
    return(data.frame(contig = names(genome),
                      length = as.numeric(genome),
                      stringsAsFactors = FALSE))
  }
  gdf <- as.data.frame(genome, stringsAsFactors = FALSE)
  contig <- .mp_pick_col(gdf, c("contig", "chrom", "seqnames", "name"))
  len <- .mp_pick_col(gdf, c("length", "len", "end", "size"))
  if (is.null(contig) || is.null(len)) {
    rlang::abort(
      c("`genome` needs `contig` and `length` columns.",
        i = "Or pass a named numeric vector of contig lengths."),
      class = "micropubr_bad_genome")
  }
  data.frame(contig = as.character(contig), length = as.numeric(len),
             stringsAsFactors = FALSE)
}

#' Coverage track for mp_circos
#' @noRd
.mp_circos_coverage <- function(coverage, coverage.type) {
  cov <- as.data.frame(coverage, stringsAsFactors = FALSE)
  cov$contig <- as.character(.mp_pick_col(cov, c("contig", "chrom",
                                                 "seqnames")))
  cov$.pos <- as.numeric(.mp_pick_col(cov, c("position", "pos", "start")))
  cov$.val <- as.numeric(.mp_pick_col(cov, c("value", "depth", "coverage",
                                             "score")))
  area <- coverage.type == "area"
  circlize::circos.track(
    factors = cov$contig, x = cov$.pos, y = cov$.val,
    track.height = 0.18, bg.border = "grey80",
    panel.fun = function(x, y) {
      ord <- order(x)
      if (area) {
        circlize::circos.lines(x[ord], y[ord], area = TRUE,
                               col = "#0072B280", border = "#0072B2")
      } else if (coverage.type == "heatmap") {
        circlize::circos.points(x[ord], rep(0.5, length(x)),
                                pch = 15, cex = 0.4)
      } else {
        circlize::circos.lines(x[ord], y[ord], col = "#0072B2")
      }
    })
}

#' GC track for mp_circos
#' @noRd
.mp_circos_gc <- function(gc) {
  g <- as.data.frame(gc, stringsAsFactors = FALSE)
  g$contig <- as.character(.mp_pick_col(g, c("contig", "chrom", "seqnames")))
  g$.pos <- as.numeric(.mp_pick_col(g, c("start", "position", "pos")))
  g$.val <- as.numeric(.mp_pick_col(g, c("value", "gc", "score")))
  circlize::circos.track(
    factors = g$contig, x = g$.pos, y = g$.val,
    track.height = 0.12, bg.border = "grey80",
    panel.fun = function(x, y) {
      ord <- order(x)
      circlize::circos.lines(x[ord], y[ord], col = "#009E73")
    })
}

#' Annotation track (feature rectangles) for mp_circos
#' @noRd
.mp_circos_annotations <- function(annotations) {
  ann <- as.data.frame(annotations, stringsAsFactors = FALSE)
  ann$contig <- as.character(.mp_pick_col(ann, c("contig", "chrom",
                                                 "seqnames")))
  ann$.start <- as.numeric(.mp_pick_col(ann, c("start", "position", "pos")))
  ann$.end <- as.numeric(.mp_pick_col(ann, c("end", "stop"),
                                      default = ann$.start))
  circlize::circos.track(
    factors = ann$contig, ylim = c(0, 1), track.height = 0.06,
    bg.border = NA,
    panel.fun = function(x, y) {
      si <- circlize::CELL_META$sector.index
      sub <- ann[ann$contig == si, , drop = FALSE]
      if (nrow(sub) > 0) {
        circlize::circos.rect(sub$.start, 0.2, sub$.end, 0.8,
                              col = "#D55E00", border = NA)
      }
    })
}

#' Ribbon links for mp_circos
#' @noRd
.mp_circos_links <- function(links) {
  lk <- as.data.frame(links, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(lk))) {
    circlize::circos.link(
      lk$contig1[i], c(lk$start1[i], lk$end1[i]),
      lk$contig2[i], c(lk$start2[i], lk$end2[i]),
      col = "#88888880", border = NA)
  }
}
