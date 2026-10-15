# Sequence quality control -----------------------------------------------
# A pre-processing companion to the community-ecology functions: per-cycle
# read-quality profiles drawn straight from FASTQ files. This is the figure
# you inspect to choose DADA2 truncation lengths *before* building a phyloseq
# object. The per-cycle quality summary is read with `ShortRead::qa()` (the
# same backend `dada2::plotQualityProfile()` uses), so only the lightweight
# `ShortRead` package is needed rather than all of `dada2`.

#' Per-cycle sequence quality profile
#'
#' Draws the read-quality profile across sequencing cycles for one or more
#' FASTQ files -- the diagnostic used to choose DADA2/QIIME 2 truncation
#' lengths. For each cycle (base position) it shows the distribution of Phred
#' quality scores as a frequency heatmap (`"heatmap"`), or the per-cycle mean
#' with an inter-quartile ribbon (`"line"`), overlaid with mean / quartile
#' summary lines. A horizontal reference line marks a target quality (Q30 by
#' default).
#'
#' Input can be FASTQ file paths (optionally `.gz`), a directory of FASTQ
#' files, or a pre-computed per-cycle data frame with columns `file`, `cycle`,
#' `score` and `count` (the shape returned by `return.data = TRUE`), so a
#' summary can be saved once and re-plotted without re-reading the reads.
#' Reading FASTQ files requires the `ShortRead` package.
#'
#' @param x A character vector of FASTQ file paths (wildcards allowed), a
#'   single directory containing FASTQ files, or a per-cycle quality data frame
#'   (`file`, `cycle`, `score`, `count`).
#' @param aggregate Pool all files into a single combined profile.
#' @param n Maximum number of reads to sample per file when reading FASTQ.
#' @param sample.names Optional display names, one per file (defaults to the
#'   file names with the FASTQ extension stripped).
#' @param plot.type `"heatmap"` (per-cycle score frequencies) or `"line"`
#'   (mean with an inter-quartile ribbon).
#' @param show Summary overlays to draw: any of `"mean"` and `"quartiles"`.
#' @param quality.line Phred score for the horizontal reference line (e.g.
#'   `30`); `NULL` to omit it.
#' @param palette,theme,plot.title,plot.subtitle,plot.caption,font.size See
#'   [mp_alpha()].
#' @param ... Passed to the primary geom.
#' @param interactive If `TRUE`, return an interactive `plotly` widget instead
#'   of a static `ggplot` (requires `plotly`). See [mp_alpha()].
#' @param return.data If `TRUE`, return the tidy per-cycle quality data frame
#'   (`file`, `cycle`, `score`, `count`) instead of a plot.
#' @param verbose Print the methods used.
#'
#' @return A `ggplot` object; a `plotly` widget when `interactive = TRUE`; or a
#'   data frame when `return.data = TRUE`.
#' @export
#' @seealso [mp_import_dada2()] for turning denoised reads into a phyloseq
#'   object once truncation lengths are chosen.
#' @examples
#' # Usually you pass FASTQ paths, e.g.
#' #   mp_qc_profile(c("sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz"))
#' # which needs the 'ShortRead' package. Here we build an equivalent per-cycle
#' # summary by hand so the example runs anywhere:
#' set.seed(1)
#' grid <- expand.grid(file = c("R1.fastq.gz", "R2.fastq.gz"),
#'                     cycle = 1:120, score = 2:40)
#' mu <- 38 - 0.06 * grid$cycle - ifelse(grid$file == "R2.fastq.gz", 5, 0)
#' grid$count <- round(3000 * stats::dnorm(grid$score, mu, 2.5))
#' mp_qc_profile(grid[grid$count > 0, ])
mp_qc_profile <- function(x,
                          aggregate = FALSE,
                          n = 5e5,
                          sample.names = NULL,
                          plot.type = c("heatmap", "line"),
                          show = c("mean", "quartiles"),
                          quality.line = 30,
                          palette = NULL,
                          theme = NULL,
                          ...,
                          plot.title = NULL,
                          plot.subtitle = NULL,
                          plot.caption = NULL,
                          font.size = getOption("micropubr.font.size", 11),
                          interactive = FALSE,
                          return.data = FALSE,
                          verbose = getOption("micropubr.verbose", FALSE)) {
  plot.type <- match.arg(plot.type)
  show <- match.arg(show, several.ok = TRUE)
  theme <- .mp_theme(theme, font.size)

  df <- .mp_qc_input(x, n, sample.names, verbose)
  if (isTRUE(aggregate)) {
    df <- stats::aggregate(count ~ cycle + score, data = df, FUN = sum)
    df$file <- "All reads"
  }
  df <- df[order(df$file, df$cycle, df$score), , drop = FALSE]
  df$file <- factor(df$file, levels = unique(df$file))

  if (isTRUE(return.data)) {
    return(df[, c("file", "cycle", "score", "count")])
  }

  summ <- .mp_qc_summary(df, show)
  files <- levels(df$file)
  reads <- .mp_qc_read_totals(df)

  if (verbose) {
    rlang::inform(sprintf(
      "mp_qc_profile(): %d file(s), %d cycle(s), %s reads profiled.",
      length(files), length(unique(df$cycle)),
      formatC(sum(reads), format = "d", big.mark = ",")))
  }

  p <- switch(
    plot.type,
    heatmap = .mp_qc_heatmap(df, summ, palette, font.size, dots = list(...)),
    line    = .mp_qc_line(df, summ, show, palette, dots = list(...)))

  p <- .mp_qc_quality_line(p, quality.line, df)

  caption <- plot.caption %||% sprintf(
    "%d file(s); up to %s reads sampled per file.",
    length(files), formatC(n, format = "d", big.mark = ","))
  p <- p +
    ggplot2::labs(x = "Cycle (sequencing position)",
                  y = "Quality score (Phred)",
                  title = plot.title, subtitle = plot.subtitle,
                  caption = caption) +
    theme
  if (length(files) < 2) p <- p + ggplot2::theme(legend.position = "right")
  p <- .mp_tag(p, "mp_qc_profile", params = list(plot.type = plot.type))
  if (isTRUE(interactive)) return(.mp_interactive(p))
  p
}

# --- input handling ------------------------------------------------------

#' Normalise mp_qc_profile input to a per-cycle data frame
#'
#' Returns a data frame with columns `file`, `cycle`, `score`, `count`.
#' @noRd
.mp_qc_input <- function(x, n, sample.names, verbose) {
  if (is.data.frame(x)) return(.mp_qc_from_df(x))
  if (!is.character(x)) {
    rlang::abort(
      "`x` must be FASTQ file path(s), a directory, or a per-cycle data frame.",
      class = "micropubr_bad_input")
  }
  files <- .mp_qc_resolve_files(x)
  names <- sample.names %||% .mp_qc_default_names(files)
  if (length(names) != length(files)) {
    rlang::abort(
      sprintf("`sample.names` has length %d but there are %d file(s).",
              length(names), length(files)),
      class = "micropubr_bad_input")
  }
  .mp_qc_from_files(files, n, names, verbose)
}

#' Standardise a user-supplied per-cycle quality data frame
#' @noRd
.mp_qc_from_df <- function(x) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  file <- .mp_pick_col(df, c("file", "File", "sample", "Sample", "lane"),
                       default = "sample1")
  cycle <- .mp_pick_col(df, c("cycle", "Cycle", "position", "pos"))
  score <- .mp_pick_col(df, c("score", "Score", "quality", "Quality"))
  count <- .mp_pick_col(df, c("count", "Count", "n", "reads"), default = 1)
  if (is.null(cycle) || is.null(score)) {
    rlang::abort(
      c("Could not find cycle and score columns in the data frame.",
        i = "Expected columns like 'cycle'/'position' and 'score'/'quality'."),
      class = "micropubr_bad_qc_cols")
  }
  data.frame(file = as.character(file), cycle = as.integer(cycle),
             score = as.numeric(score), count = as.numeric(count),
             stringsAsFactors = FALSE)
}

#' Expand a directory / wildcard / vector of FASTQ paths to existing files
#' @noRd
.mp_qc_resolve_files <- function(x) {
  if (length(x) == 1 && dir.exists(x)) {
    files <- list.files(
      x, pattern = "\\.(fastq|fq)(\\.gz)?$", full.names = TRUE,
      ignore.case = TRUE)
    if (!length(files)) {
      rlang::abort(sprintf("No FASTQ files found in directory '%s'.", x),
                   class = "micropubr_no_fastq")
    }
    return(files)
  }
  files <- unlist(lapply(x, function(p) {
    g <- Sys.glob(p)
    if (length(g)) g else p
  }), use.names = FALSE)
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    rlang::abort(
      c("Some FASTQ file(s) do not exist:",
        stats::setNames(missing, rep("x", length(missing)))),
      class = "micropubr_no_fastq")
  }
  files
}

#' Default display names: file name without the FASTQ extension
#' @noRd
.mp_qc_default_names <- function(files) {
  nm <- basename(files)
  sub("\\.(fastq|fq)(\\.gz)?$", "", nm, ignore.case = TRUE)
}

#' Read per-cycle quality from FASTQ files via ShortRead::qa()
#' @noRd
.mp_qc_from_files <- function(files, n, names, verbose) {
  .mp_require("ShortRead", "Reading FASTQ quality profiles", bioc = TRUE)
  parts <- lapply(seq_along(files), function(i) {
    if (verbose) rlang::inform(sprintf("  reading %s", files[i]))
    srqa <- ShortRead::qa(files[i], n = n)
    q <- srqa[["perCycle"]]$quality
    data.frame(file = names[i],
               cycle = as.integer(as.character(q$Cycle)),
               score = as.numeric(as.character(q$Score)),
               count = as.numeric(q$Count),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, parts)
}

# --- summaries -----------------------------------------------------------

#' Per-cycle mean and quartiles, weighted by read count
#'
#' Returns a long data frame (`file`, `cycle`, `value`, `stat`) restricted to
#' the requested overlays, ready to draw as summary lines.
#' @noRd
.mp_qc_summary <- function(df, show) {
  key <- paste(df$file, df$cycle, sep = "\r")
  parts <- lapply(split(df, key), function(d) {
    tot <- sum(d$count)
    mean <- if (tot > 0) sum(d$score * d$count) / tot else NA_real_
    q <- .mp_qc_wquantile(d$score, d$count, c(0.25, 0.5, 0.75))
    data.frame(file = d$file[1], cycle = d$cycle[1],
               Mean = mean, Q25 = q[1], Median = q[2], Q75 = q[3],
               stringsAsFactors = FALSE)
  })
  wide <- do.call(rbind, parts)
  wide <- wide[order(wide$file, wide$cycle), , drop = FALSE]

  stats_keep <- c(if ("quartiles" %in% show) c("Q25", "Median", "Q75"),
                  if ("mean" %in% show) "Mean")
  long <- do.call(rbind, lapply(stats_keep, function(s) {
    data.frame(file = wide$file, cycle = wide$cycle, value = wide[[s]],
               stat = s, stringsAsFactors = FALSE)
  }))
  long$file <- factor(long$file, levels = levels(df$file))
  long$stat <- factor(long$stat, levels = c("Q25", "Median", "Q75", "Mean"))
  long
}

#' Weighted quantiles from (score, count) pairs
#' @noRd
.mp_qc_wquantile <- function(score, count, probs) {
  o <- order(score)
  score <- score[o]
  cw <- cumsum(count[o])
  tot <- cw[length(cw)]
  if (is.na(tot) || tot <= 0) return(rep(NA_real_, length(probs)))
  vapply(probs, function(p) score[which(cw >= p * tot)[1]], numeric(1))
}

#' Reads per file (the per-cycle count total at the first cycle)
#' @noRd
.mp_qc_read_totals <- function(df) {
  per_cycle <- stats::aggregate(count ~ file + cycle, data = df, FUN = sum)
  tapply(per_cycle$count, per_cycle$file, max, na.rm = TRUE)
}

# --- geometries ----------------------------------------------------------

#' Frequency-heatmap profile with summary lines
#' @noRd
.mp_qc_heatmap <- function(df, summ, palette, font.size, dots = list()) {
  df$.prop <- df$count / stats::ave(df$count, df$file, df$cycle, FUN = sum)
  cols <- .mp_resolve_palette(palette, 256, type = "sequential")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$cycle, y = .data$score)) +
    .mp_geom(ggplot2::geom_tile,
             params = list(mapping = ggplot2::aes(fill = .data$.prop)),
             dots = dots) +
    ggplot2::scale_fill_gradientn(colours = cols, name = "Proportion") +
    ggplot2::geom_line(
      data = summ,
      ggplot2::aes(x = .data$cycle, y = .data$value,
                   colour = .data$stat, linetype = .data$stat),
      linewidth = 0.6, na.rm = TRUE)
  p <- p + .mp_qc_stat_scales(summ)
  if (nlevels(df$file) > 1) p <- p + ggplot2::facet_wrap(~ file)
  p
}

#' Mean + inter-quartile ribbon profile (compact multi-file comparison)
#' @noRd
.mp_qc_line <- function(df, summ, show, palette, dots = list()) {
  wide <- stats::reshape(
    summ, idvar = c("file", "cycle"), timevar = "stat",
    direction = "wide")
  names(wide) <- sub("^value\\.", "", names(wide))
  cols <- .mp_resolve_palette(palette, nlevels(df$file), type = "qualitative")

  p <- ggplot2::ggplot(wide, ggplot2::aes(x = .data$cycle))
  if ("quartiles" %in% show && all(c("Q25", "Q75") %in% names(wide))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$Q25, ymax = .data$Q75, fill = .data$file),
      alpha = 0.25)
  }
  yvar <- if ("Mean" %in% names(wide)) "Mean" else "Median"
  p <- p +
    .mp_geom(ggplot2::geom_line,
             params = list(mapping = ggplot2::aes(y = .data[[yvar]],
                                                  colour = .data$file),
                           linewidth = 0.7),
             dots = dots) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_fill_manual(values = cols, name = NULL)
  p
}

#' Colour/linetype scales for the mean & quartile summary lines
#' @noRd
.mp_qc_stat_scales <- function(summ) {
  cols <- c(Mean = "#009E73", Median = "#D55E00",
            Q25 = "#E69F00", Q75 = "#E69F00")
  ltys <- c(Mean = "solid", Median = "solid", Q25 = "22", Q75 = "22")
  present <- levels(droplevels(summ$stat))
  list(
    ggplot2::scale_colour_manual(values = cols[present], name = NULL,
                                 breaks = present),
    ggplot2::scale_linetype_manual(values = ltys[present], name = NULL,
                                   breaks = present))
}

#' Horizontal reference line at a target Phred score
#' @noRd
.mp_qc_quality_line <- function(p, quality.line, df) {
  if (is.null(quality.line)) return(p)
  p +
    ggplot2::geom_hline(yintercept = quality.line, linetype = "dotted",
                        colour = "grey40", linewidth = 0.4) +
    ggplot2::annotate("text", x = max(df$cycle, na.rm = TRUE),
                      y = quality.line, label = paste0("Q", quality.line),
                      hjust = 1.05, vjust = -0.4, size = 2.8,
                      colour = "grey40")
}
