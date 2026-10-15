# Reporting: provenance, methods text, alt text and guardrails -------------
# Every plotting function records what it did -- the analysis parameters, the
# statistics it computed, the data it saw and any caveats -- on the returned
# plot (attribute "mp_provenance") and in a session log. mp_methods() turns
# those records into a methods paragraph with references, mp_alt_text() into
# an accessible description. The guardrails are the caveats: conditions under
# which a figure's statistics are easy to over-read, reported where the
# analyst will see them rather than left for a reviewer to find.

.mp_session <- new.env(parent = emptyenv())
.mp_session$log <- list()

#' Write a methods paragraph for micropubr figures
#'
#' Assembles a ready-to-edit methods paragraph describing how the supplied
#' figures were made: the statistical tests and their settings, the
#' transformations, filters and multiple-testing corrections, sample sizes,
#' software versions and the references to cite. It is built from the
#' provenance each `mp_*()` function attaches to its plot, so it always
#' matches what was actually computed. Statistical caveats raised while the
#' figures were made (see the guardrails in [mp_set_defaults()]) are returned
#' as notes.
#'
#' The paragraph covers the statistical-analysis and software items that
#' reporting guidelines such as STORMS (Mirzayi et al. 2021, *Nature
#' Medicine* 27:1885) ask for. It does not know your study design, sampling or
#' sequencing, so treat it as a draft to complete, not a finished section.
#'
#' @param ... Figures made by `mp_*()` functions: `ggplot`/`patchwork`
#'   objects, interactive widgets, or lists of them. With no figures, every
#'   figure made in the current R session is described.
#' @param format `"text"` (plain text with references), `"markdown"`, or
#'   `"list"` (the components, for programmatic use).
#' @param references Include the reference list.
#'
#' @return For `"text"` and `"markdown"`, a character string of class
#'   `mp_methods` (printed nicely). For `"list"`, the components: `text` (the
#'   methods paragraph, without the software sentence), `software` (that
#'   sentence), `versions` (named package versions, including R), `citations`
#'   (in-text citation labels, named by reference key), `references` and
#'   `notes` (statistical caveats). Keeping the software apart lets a caller
#'   that describes figures one at a time write a single software sentence for
#'   all of them.
#' @export
#' @seealso [mp_alt_text()] for figure alt text.
#' @examples
#' \donttest{
#' data(GlobalPatterns, package = "phyloseq")
#' p1 <- mp_alpha(GlobalPatterns, measure = "Shannon", group.by = "SampleType",
#'                test = "kruskal")
#' p2 <- mp_ordination(GlobalPatterns, group.by = "SampleType")
#' mp_methods(p1, p2)
#' }
mp_methods <- function(..., format = c("text", "markdown", "list"),
                       references = TRUE) {
  format <- match.arg(format)
  figs <- list(...)
  provs <- if (length(figs)) .mp_collect_provenance(figs) else
    .mp_session$log
  if (!length(provs)) {
    rlang::abort(
      c("No micropubr figures to describe.",
        i = "Pass figures made by mp_*() functions, e.g. mp_methods(p1, p2)."),
      class = "micropubr_no_provenance")
  }

  keys <- character(0)
  cite <- function(key) {
    keys <<- c(keys, key)
    .mp_ref(key)$short
  }
  sentences <- unique(unlist(lapply(provs, .mp_methods_sentences,
                                    cite = cite)))
  versions <- unlist(lapply(provs, function(p) p$versions))
  versions <- versions[!duplicated(names(versions)) & !is.na(versions)]
  software <- .mp_software_sentence(versions, cite)
  text <- paste(c(sentences, software), collapse = " ")
  refs <- vapply(unique(keys), function(k) .mp_ref(k, versions)$full,
                 character(1))
  refs <- sort(unname(refs))
  notes <- unique(unlist(lapply(provs, function(p) p$notes)))

  if (format == "list") {
    citations <- vapply(unique(keys), function(k) .mp_ref(k, versions)$short,
                        character(1))
    return(list(text = paste(sentences, collapse = " "), references = refs,
                notes = notes, software = software,
                versions = c(R = paste0(R.version$major, ".", R.version$minor),
                             versions),
                citations = citations))
  }
  md <- format == "markdown"
  out <- text
  if (references && length(refs)) {
    head <- if (md) "**References**" else "References"
    items <- if (md) paste0("- ", refs) else refs
    out <- paste(c(out, "", head, items), collapse = "\n")
  }
  if (length(notes)) {
    head <- if (md) "**Statistical notes**" else "Statistical notes"
    items <- if (md) paste0("- ", notes) else paste0("* ", notes)
    out <- paste(c(out, "", head, items), collapse = "\n")
  }
  structure(out, class = c("mp_methods", "character"))
}

#' @export
print.mp_methods <- function(x, ...) {
  paras <- strsplit(unclass(x), "\n", fixed = TRUE)[[1]]
  for (para in paras) {
    if (!nzchar(para)) {
      cat("\n")
    } else {
      cat(strwrap(para, width = getOption("width", 80), exdent = 2),
          sep = "\n")
    }
  }
  invisible(x)
}

#' Generate alt text for a micropubr figure
#'
#' Describes a figure for readers who cannot see it: the chart type, what is
#' plotted against what, how it is grouped, how much data it shows, the key
#' result (e.g. which group has the highest diversity, how many taxa are
#' differentially abundant) and the statistics in its subtitle and caption.
#' Journals increasingly ask for alt text with every figure.
#'
#' @param plot A figure made by an `mp_*()` function (any `ggplot` works, with
#'   a more generic description).
#' @return A character string.
#' @export
#' @examples
#' \donttest{
#' data(GlobalPatterns, package = "phyloseq")
#' p <- mp_alpha(GlobalPatterns, measure = "Shannon", group.by = "SampleType",
#'               test = "kruskal")
#' mp_alt_text(p)
#' }
mp_alt_text <- function(plot) {
  if (inherits(plot, "patchwork")) {
    panels <- .mp_patch_plots(plot)
    prov <- attr(plot, "mp_provenance")
    panels <- lapply(panels, function(pn) {
      if (is.null(attr(pn, "mp_provenance"))) attr(pn, "mp_provenance") <- prov
      pn
    })
    if (length(panels) > 1) {
      parts <- vapply(seq_along(panels), function(i) {
        sprintf("Panel %s: %s", LETTERS[i], mp_alt_text(panels[[i]]))
      }, character(1))
      return(paste(c(sprintf("A figure with %d panels.", length(panels)),
                     parts), collapse = " "))
    }
  }
  if (!.mp_is_ggplot(plot)) {
    rlang::abort("`plot` must be a ggplot made by an mp_*() function.",
                 class = "micropubr_bad_input")
  }
  prov <- attr(plot, "mp_provenance")
  generic <- tryCatch(.mp_alt_generic(plot), error = function(e) NULL)
  specific <- if (!is.null(prov)) {
    tryCatch(.mp_alt_specific(prov, plot), error = function(e) NULL)
  }
  stats <- c(.mp_label_text(plot, "subtitle"), .mp_label_text(plot, "caption"))
  stats <- stats[nzchar(stats)]
  stats <- if (length(stats)) {
    lines <- unlist(strsplit(stats, "\n", fixed = TRUE))
    lines <- sub("\\.\\s*$", "", trimws(lines[nzchar(trimws(lines))]))
    paste0("Statistics: ", paste(lines, collapse = "; "), ".")
  }
  out <- paste(c(generic, specific, stats), collapse = " ")
  gsub("\\.\\.", ".", out)
}

# --- Provenance ------------------------------------------------------------

#' Attach a provenance record to a plot and the session log
#'
#' @param emit Raise the guardrail `notes` now (FALSE when they were already
#'   raised upstream, e.g. by mp_run_da()).
#' @noRd
.mp_tag <- function(p, fun, params = list(), stats = list(), data = list(),
                    notes = character(0), packages = character(0),
                    emit = TRUE) {
  notes <- unique(notes[!is.na(notes) & nzchar(notes)])
  prov <- list(
    fun = fun, params = params, stats = stats, data = data, notes = notes,
    caption = .mp_label_text(p, "caption"),
    subtitle = .mp_label_text(p, "subtitle"),
    versions = .mp_versions(unique(c("micropubr", "ggplot2", packages)))
  )
  if (emit) p <- .mp_guard_emit(notes, p)
  .mp_session_log(prov)
  .mp_set_provenance(p, prov)
}

#' @noRd
.mp_set_provenance <- function(p, prov) {
  if (is.list(p) && !.mp_is_ggplot(p) && !inherits(p, "htmlwidget")) {
    return(lapply(p, .mp_set_provenance, prov = prov))
  }
  attr(p, "mp_provenance") <- prov
  p
}

#' @noRd
.mp_session_log <- function(prov) {
  log <- c(.mp_session$log, list(prov))
  .mp_session$log <- utils::tail(log, 200)
  invisible(NULL)
}

#' Provenance records from figures, lists and patchworks
#' @noRd
.mp_collect_provenance <- function(x) {
  if (is.null(x)) return(list())
  own <- attr(x, "mp_provenance")
  if (inherits(x, "patchwork")) {
    inner <- unlist(lapply(.mp_patch_plots(x), .mp_collect_provenance),
                    recursive = FALSE)
    return(c(if (!is.null(own)) list(own), inner))
  }
  if (!is.null(own)) return(list(own))
  if (is.list(x) && !.mp_is_ggplot(x) && !inherits(x, "htmlwidget")) {
    return(unlist(lapply(x, .mp_collect_provenance), recursive = FALSE))
  }
  list()
}

#' The component plots of a patchwork (the last plot is the object itself)
#' @noRd
.mp_patch_plots <- function(x) {
  inner <- tryCatch(x$patches$plots, error = function(e) NULL)
  last <- x
  class(last) <- setdiff(class(last), "patchwork")
  c(inner, list(last))
}

#' A plot label as a single string ("" when absent or not text)
#'
#' @param resolve Also resolve labels ggplot2 only fills in when the plot is
#'   drawn (legend titles taken from the mapping, in ggplot2 >= 4.0), via
#'   ggplot2::get_labs() where it exists.
#' @noRd
.mp_label_text <- function(p, which, resolve = FALSE) {
  lab <- tryCatch({
    if (!.mp_is_ggplot(p)) NULL else {
      labs <- p$labels
      ns <- asNamespace("ggplot2")
      if (resolve && exists("get_labs", ns, inherits = FALSE)) {
        labs <- get("get_labs", envir = ns)(p)
      }
      labs[[which]]
    }
  }, error = function(e) NULL)
  if (is.null(lab)) return("")
  if (is.language(lab) || is.expression(lab)) lab <- .mp_deparse_label(lab)
  paste(as.character(lab), collapse = " ")
}

#' Turn a plotmath label into plain text
#' @noRd
.mp_deparse_label <- function(lab) {
  txt <- paste(deparse(if (is.expression(lab)) lab[[1]] else lab),
               collapse = "")
  txt <- gsub("\\[([^]]+)\\]", "\\1", txt)
  txt <- gsub("italic\\(([^)]+)\\)", "\\1", txt)
  txt <- gsub("~|\\*", " ", txt)
  txt <- gsub(" - ", "-", txt, fixed = TRUE)
  txt <- gsub("\"", "", txt)
  trimws(gsub("\\s+", " ", txt))
}

#' Sample/taxon counts and group sizes for a provenance record
#' @noRd
.mp_data_info <- function(object, group.by = NULL) {
  out <- list(n_samples = phyloseq::nsamples(object),
              n_taxa = phyloseq::ntaxa(object))
  if (!is.null(group.by)) {
    g <- .mp_sample_data(object)[[group.by]]
    tab <- table(g[!is.na(g)])
    out$groups <- stats::setNames(as.integer(tab), names(tab))
  }
  out
}

# --- Guardrails ------------------------------------------------------------

#' Raise guardrail notes as warnings, or add them to the caption
#'
#' Controlled by `mp_set_defaults(guardrails = "warn" | "caption" | "off")`.
#' @noRd
.mp_guard_emit <- function(notes, p = NULL) {
  notes <- unique(notes[!is.na(notes) & nzchar(notes)])
  mode <- getOption("micropubr.guardrails", "warn")
  if (!length(notes) || identical(mode, "off")) return(p)
  if (identical(mode, "caption") && !is.null(p) && .mp_is_ggplot(p) &&
      !inherits(p, "patchwork")) {
    cap <- .mp_label_text(p, "caption")
    cap <- paste(c(if (nzchar(cap)) cap, paste("Note:", notes)),
                 collapse = "\n")
    return(p + ggplot2::labs(caption = cap))
  }
  for (n in notes) rlang::warn(n, class = "micropubr_guardrail")
  p
}

#' Too few samples per group for a test to mean much
#' @noRd
.mp_guard_groups <- function(object, group.by, min.n = 5) {
  if (is.null(group.by)) return(character(0))
  g <- .mp_sample_data(object)[[group.by]]
  n <- table(g[!is.na(g)])
  low <- n[n < min.n]
  if (length(n) < 2 || !length(low)) return(character(0))
  sprintf(paste(
    "Group(s) %s have fewer than %d samples; group comparisons have little",
    "power and their p-values are unstable at this size."),
    paste(sprintf("%s (n = %d)", names(low), as.integer(low)),
          collapse = ", "), min.n)
}

#' Richness compared across very uneven sequencing depths
#' @noRd
.mp_guard_depth <- function(object, measure) {
  rich <- intersect(measure, c("Observed", "Chao1", "ACE", "Fisher"))
  if (!length(rich)) return(character(0))
  d <- phyloseq::sample_sums(object)
  d <- d[d > 0]
  if (length(d) < 2 || max(d) / min(d) <= 10) return(character(0))
  sprintf(paste(
    "Sequencing depth varies %.0f-fold across samples; richness (%s)",
    "depends strongly on depth, so consider rarefying or comparing at equal",
    "coverage."), max(d) / min(d), paste(rich, collapse = ", "))
}

#' PERMDISP check that accompanies a PERMANOVA
#'
#' @return list(p, note) or NULL when it cannot be computed.
#' @noRd
.mp_guard_dispersion <- function(dist, meta, group, permutations) {
  g <- factor(meta[labels(dist), group])
  keep <- !is.na(g)
  if (nlevels(droplevels(g[keep])) < 2) return(NULL)
  d <- if (all(keep)) dist else
    stats::as.dist(as.matrix(dist)[keep, keep])
  p <- tryCatch({
    bd <- vegan::betadisper(d, droplevels(g[keep]))
    vegan::permutest(bd, permutations = permutations)$tab[["Pr(>F)"]][1]
  }, error = function(e) NA_real_)
  if (is.na(p)) return(NULL)
  note <- if (p < 0.05) {
    sprintf(paste(
      "Group dispersions differ (PERMDISP p = %s); a significant PERMANOVA",
      "may reflect differences in spread rather than in location",
      "(centroids)."), format.pval(p, digits = 2))
  } else character(0)
  list(p = p, note = note)
}

# --- Methods sentences ------------------------------------------------------

#' Sentences for one provenance record
#' @noRd
.mp_methods_sentences <- function(prov, cite) {
  f <- switch(
    prov$fun,
    mp_alpha = .mp_m_alpha,
    mp_ordination = .mp_m_ordination,
    mp_dispersion = .mp_m_dispersion,
    mp_taxa = .mp_m_taxa,
    mp_prevalence = .mp_m_prevalence,
    mp_indicator = .mp_m_indicator,
    mp_mantel = .mp_m_mantel,
    mp_network = .mp_m_network,
    mp_ml = .mp_m_ml,
    mp_diff = .mp_m_diff,
    mp_varpart = .mp_m_varpart,
    mp_longitudinal = .mp_m_longitudinal,
    mp_tree = .mp_m_tree,
    .mp_m_generic)
  tryCatch(f(prov, cite), error = function(e) {
    rlang::warn(
      sprintf("Could not describe the %s() figure in detail (%s); using its caption.",
              prov$fun, conditionMessage(e)),
      class = "micropubr_methods_fallback")
    .mp_m_generic(prov, cite)
  })
}

#' "a, b and c"
#' @noRd
.mp_and <- function(x) {
  x <- as.character(x)
  if (length(x) <= 1) return(paste(x, collapse = ""))
  paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
}

#' p-value adjustment as prose, citing the procedure
#' @noRd
.mp_padj_phrase <- function(method, cite) {
  switch(method %||% "BH",
         BH = , fdr = sprintf("the Benjamini-Hochberg procedure (%s)",
                              cite("bh")),
         holm = sprintf("the Holm procedure (%s)", cite("holm")),
         bonferroni = "the Bonferroni correction",
         BY = "the Benjamini-Yekutieli procedure",
         hochberg = "the Hochberg procedure",
         hommel = "the Hommel procedure",
         none = "no correction",
         method)
}

#' @noRd
.mp_group_phrase <- function(data) {
  g <- data$groups
  if (!length(g)) return("")
  sprintf(" (%s)", paste(sprintf("%s, n = %d", names(g), g), collapse = "; "))
}

#' @noRd
.mp_distance_phrase <- function(distance, cite) {
  switch(distance,
         bray = "Bray-Curtis dissimilarity",
         jaccard = "Jaccard distance",
         unifrac = sprintf("unweighted UniFrac distance (%s)",
                           cite("unifrac")),
         wunifrac = sprintf("weighted UniFrac distance (%s)",
                            cite("unifrac")),
         aitchison = sprintf(paste("Aitchison distance (Euclidean distance",
                                   "of centred log-ratio-transformed counts;",
                                   "%s)"), cite("aitchison")),
         sprintf("%s distance", distance))
}

#' @noRd
.mp_m_alpha <- function(prov, cite) {
  pa <- prov$params
  d <- prov$data
  if (identical(pa$plot.type, "rarefaction")) {
    return(sprintf(
      "Rarefaction curves of observed richness were computed for %d samples with vegan (%s).",
      d$n_samples, cite("vegan")))
  }
  s <- sprintf("Alpha diversity (%s) was estimated for %d samples with phyloseq (%s).",
               .mp_and(pa$measure), d$n_samples, cite("phyloseq"))
  k <- length(d$groups)
  if (!isTRUE(pa$add.stats) || k < 2) return(s)
  test <- pa$test
  if (k == 2) {
    glob <- if (test == "t.test") {
      "Welch's t-test"
    } else sprintf("the Wilcoxon rank-sum test (%s)", cite("wilcoxon"))
    return(c(s, sprintf("Differences between %s groups%s were tested with %s.",
                        pa$group.by, .mp_group_phrase(d), glob)))
  }
  glob <- if (test == "anova") "one-way ANOVA" else
    sprintf("the Kruskal-Wallis test (%s)", cite("kruskal"))
  post <- switch(
    test,
    kruskal = sprintf("Dunn's tests (%s)", cite("dunn")),
    anova = "Tukey's honest significant difference test",
    wilcoxon = sprintf("Wilcoxon rank-sum tests (%s)", cite("wilcoxon")),
    t.test = "Welch's t-tests")
  adj <- if (test == "anova") "" else
    sprintf(" adjusted with %s", .mp_padj_phrase(pa$p.adjust.method, cite))
  c(s, sprintf(
    "Differences among the %d %s groups%s were tested with %s, followed by pairwise %s%s (rstatix; %s).",
    k, pa$group.by, .mp_group_phrase(d), glob, post, adj, cite("rstatix")))
}

#' @noRd
.mp_m_ordination <- function(prov, cite) {
  pa <- prov$params
  st <- prov$stats
  d <- prov$data
  uses_dist <- pa$method %in% c("PCoA", "NMDS", "dbRDA")
  meth <- switch(
    pa$method,
    PCoA = "principal coordinates analysis (PCoA)",
    NMDS = sprintf("non-metric multidimensional scaling (NMDS; vegan::metaMDS, two dimensions%s; %s)",
                   if (!is.null(st$stress)) sprintf(", stress = %.3f", st$stress) else "",
                   cite("vegan")),
    PCA = "principal component analysis (PCA) of centred log-ratio-transformed counts",
    RDA = sprintf("redundancy analysis (RDA; vegan; %s)", cite("vegan")),
    CCA = sprintf("canonical correspondence analysis (CCA; vegan; %s)", cite("vegan")),
    dbRDA = sprintf("distance-based redundancy analysis (dbRDA; vegan; %s)", cite("vegan")))
  s <- sprintf("Between-sample (beta) diversity of %d samples was visualised by %s%s%s.",
               d$n_samples, meth,
               if (uses_dist) sprintf(" of %s computed with phyloseq (%s)",
                                      .mp_distance_phrase(pa$distance, cite),
                                      cite("phyloseq")) else "",
               if (length(pa$constraints)) sprintf(", constrained by %s",
                                                   .mp_and(pa$constraints)) else "")
  if (pa$method == "PCA") s <- sub("\\.$", sprintf(" (%s).", cite("aitchison")), s)
  out <- s
  if (!is.null(st$permanova)) {
    out <- c(out, sprintf(
      "Differences in community composition among %s groups%s were tested by permutational multivariate analysis of variance (PERMANOVA; vegan::adonis2, %d permutations; %s): R\u00b2 = %.3f, F = %.2f, p = %s.",
      pa$group.by, .mp_group_phrase(d), pa$permutations, cite("permanova"),
      st$permanova$R2, st$permanova$F, format.pval(st$permanova$p, digits = 2)))
  }
  if (!is.null(st$permdisp_p)) {
    out <- c(out, sprintf(
      "Homogeneity of multivariate dispersions, which PERMANOVA assumes, was assessed with PERMDISP (vegan::betadisper; %s): p = %s.",
      cite("permdisp"), format.pval(st$permdisp_p, digits = 2)))
  }
  if (length(pa$envfit)) {
    out <- c(out, sprintf(
      "Metadata vectors (%s) were fitted onto the ordination with vegan::envfit (%d permutations).",
      .mp_and(pa$envfit), pa$permutations))
  }
  out
}

#' @noRd
.mp_m_dispersion <- function(prov, cite) {
  pa <- prov$params
  st <- prov$stats
  s <- sprintf(
    "Multivariate dispersion (the distance of each sample to its group centroid, based on %s) was compared among %s groups%s with PERMDISP (vegan::betadisper and permutest, %d permutations; %s): F = %.2f, p = %s.",
    .mp_distance_phrase(pa$distance, cite), pa$group.by,
    .mp_group_phrase(prov$data), pa$permutations, cite("permdisp"),
    st$F, format.pval(st$p, digits = 2))
  extra <- c(anosim = "ANOSIM", mrpp = "MRPP", permanova = "PERMANOVA")
  if (length(pa$add.test)) {
    if ("permanova" %in% pa$add.test) cite("permanova")
    s <- c(s, sprintf("Group differences were additionally tested with %s (vegan, %d permutations).",
                      .mp_and(extra[pa$add.test]), pa$permutations))
  }
  s
}

#' @noRd
.mp_m_taxa <- function(prov, cite) {
  pa <- prov$params
  what <- if (isTRUE(pa$relative)) "relative abundance" else "read counts"
  rank <- if (is.null(pa$tax.level)) "feature" else pa$tax.level
  sprintf(
    "Taxonomic composition is shown as %s aggregated at the %s level (phyloseq; %s); the %d most abundant taxa are shown individually%s.",
    what, rank, cite("phyloseq"), pa$top.n,
    if (pa$plot.type %in% c("bar", "area")) {
      sprintf(" and the remainder pooled as \"%s\"", pa$other.label)
    } else "")
}

#' @noRd
.mp_m_prevalence <- function(prov, cite) {
  pa <- prov$params
  sprintf(
    "Prevalence was defined as the fraction of samples in which a taxon's relative abundance exceeds %s%s; taxa with a prevalence of at least %.0f%% were considered part of the core microbiome.",
    format(pa$detection, scientific = FALSE),
    if (!is.null(pa$tax.level)) sprintf(" at the %s level", pa$tax.level) else "",
    100 * pa$prevalence)
}

#' @noRd
.mp_m_indicator <- function(prov, cite) {
  pa <- prov$params
  sprintf(
    "Taxa characteristic of each %s group were identified by indicator species analysis (IndVal.g; indicspecies::multipatt, %d permutations; %s); indicators with p < %g are shown.",
    pa$group.by, pa$permutations, cite("indicspecies"), pa$p.threshold)
}

#' @noRd
.mp_m_mantel <- function(prov, cite) {
  pa <- prov$params
  st <- prov$stats
  sprintf(
    "The correlation between community dissimilarity (%s) and %s was assessed with a Mantel test (%s correlation, %d permutations; vegan; %s, %s): r = %.3f, p = %s.",
    .mp_distance_phrase(pa$distance, cite),
    if (length(pa$vars)) sprintf("Euclidean distance in standardised %s",
                                 .mp_and(pa$vars)) else "a second distance matrix",
    pa$method, pa$permutations, cite("mantel"), cite("vegan"),
    st$r, format.pval(st$p, digits = 2))
}

#' @noRd
.mp_m_network <- function(prov, cite) {
  pa <- prov$params
  meth <- switch(pa$method, spearman = "Spearman", pearson = "Pearson",
                 sparcc = "SparCC", `spiec-easi` = "SPIEC-EASI", pa$method)
  if (pa$method %in% c("spearman", "pearson")) {
    sprintf(
      "A co-occurrence network of the %d most abundant taxa%s was inferred from pairwise %s correlations; edges were kept for |r| >= %g and p < %g after adjustment with %s.",
      pa$top.n, if (!is.null(pa$tax.level)) sprintf(" (%s level)", pa$tax.level) else "",
      meth, pa$cor.threshold, pa$p.threshold,
      .mp_padj_phrase(pa$p.adjust.method, cite))
  } else {
    sprintf("A co-occurrence network of the %d most abundant taxa was inferred with %s.",
            pa$top.n, meth)
  }
}

#' @noRd
.mp_m_ml <- function(prov, cite) {
  pa <- prov$params
  model <- if (pa$method == "rf") {
    sprintf("random forest (randomForest, %d trees; %s)", pa$ntree, cite("randomforest"))
  } else {
    sprintf("gradient-boosted trees (xgboost, %d rounds)", pa$nrounds)
  }
  transform <- switch(pa$transform %||% "clr",
                      clr = "centred log-ratio (CLR)-transformed",
                      relative = "relative", none = "untransformed",
                      pa$transform)
  s <- sprintf(
    "A %s %s model predicting %s from %s%s abundances was evaluated by %d-fold cross-validation.",
    model, pa$task, pa$group.by, transform,
    if (!is.null(pa$tax.level)) sprintf(" %s-level", pa$tax.level) else "",
    pa$cv.folds)
  if (identical(pa$plot.type, "roc")) {
    s <- c(s, sprintf("Receiver operating characteristic curves were computed with pROC (%s).",
                      cite("proc")))
  }
  if (nzchar(prov$caption %||% "")) s <- c(s, paste0(sub("\\.?$", "", prov$caption), "."))
  s
}

#' @noRd
.mp_m_diff <- function(prov, cite) {
  da <- prov$params$da
  pa <- prov$params
  thr <- sprintf(
    "features with an adjusted p-value below %g and an absolute log2 fold change above %g were considered differentially abundant",
    pa$p.threshold, pa$fc.threshold)
  if (is.null(da) || is.null(da$group.by)) {
    tools <- unique(prov$stats$methods %||% "result")
    tools[tools == "result"] <- "a user-supplied results table"
    tool <- .mp_and(tools)
    return(sprintf("Differential-abundance results from %s were visualised; %s.",
                   tool, thr))
  }
  # Versions may be missing for a method (e.g. results stitched together by
  # hand), so look them up softly: `[[` on a missing name would error.
  ver <- function(pkg) {
    v <- da$versions
    if (!is.null(v) && pkg %in% names(v) && !is.na(v[[pkg]])) v[[pkg]] else ""
  }
  phrase <- function(m) {
    switch(
      m,
      wilcox = sprintf("Wilcoxon rank-sum tests (%s) on centred log-ratio (CLR)-transformed counts (pseudocount 1; %s)",
                       cite("wilcoxon"), cite("aitchison")),
      lm = if (is.null(da$random)) {
        sprintf("linear models on CLR-transformed counts (%s)", cite("aitchison"))
      } else {
        sprintf("linear mixed models on CLR-transformed counts with a random intercept for %s (lme4 and lmerTest, Satterthwaite degrees of freedom; %s; %s)",
                da$random, cite("lme4"), cite("lmertest"))
      },
      ancombc2 = sprintf("ANCOM-BC2 (ANCOMBC %s; %s), which corrects for sample-specific sampling fractions, detects structural zeros and assesses each result's sensitivity to the pseudo-count",
                         ver("ANCOMBC"), cite("ancombc2")),
      deseq2 = sprintf("DESeq2 negative-binomial Wald tests with 'poscounts' size factors (DESeq2 %s; %s)",
                       ver("DESeq2"), cite("deseq2")),
      aldex2 = sprintf("ALDEx2 Welch's t-tests on Monte-Carlo CLR instances (ALDEx2 %s; %s)",
                       ver("ALDEx2"), cite("aldex2")),
      m)
  }
  m <- da$method
  rank <- if (is.null(da$tax.level)) "" else
    sprintf(" agglomerated at the %s level", da$tax.level)
  groups <- paste(sprintf("%s, n = %d", names(da$n), da$n), collapse = "; ")
  out <- sprintf(
    "Differential abundance among %s groups (%s) was tested for %s %s on taxa%s detected in at least %.0f%% of samples (%d taxa) using %s.",
    da$group.by, groups,
    if (length(da$contrasts) > 1) "the contrasts" else "the contrast",
    .mp_and(da$contrasts), rank, 100 * da$prevalence,
    da$n_taxa, .mp_and(vapply(m, phrase, character(1))))
  if (length(da$adjust)) {
    adj_m <- intersect(m, c("lm", "ancombc2", "deseq2"))
    un_m <- setdiff(m, adj_m)
    out <- c(out, sprintf("%s adjusted for %s%s.",
                          if (length(m) > 1) "Models were" else "The model was",
                          .mp_and(da$adjust),
                          if (length(un_m)) sprintf("; %s cannot include covariates and %s unadjusted",
                                                    .mp_and(c(wilcox = "the Wilcoxon test", aldex2 = "ALDEx2")[un_m]),
                                                    if (length(un_m) > 1) "were run" else "was run") else ""))
  }
  pair <- if (isTRUE(da$pairwise) && "ancombc2" %in% m) {
    "; ANCOM-BC2 pairwise contrasts control the mixed directional false discovery rate"
  } else ""
  out <- c(out, sprintf("P-values were adjusted within each method and contrast with %s%s; %s%s.",
                        .mp_padj_phrase(da$p.adjust.method, cite), pair, thr,
                        if (isTRUE(pa$robust.only) && "ancombc2" %in% m) ", counting ANCOM-BC2 results as significant only when they passed its pseudo-count sensitivity analysis" else ""))
  if (identical(pa$plot.type, "consensus")) {
    out <- c(out, "Agreement between methods is reported as the number of methods identifying each taxon.")
  }
  out
}

#' @noRd
.mp_m_varpart <- function(prov, cite) {
  g <- prov$params$groups
  sprintf(
    "The variation in community composition (Hellinger-transformed) was partitioned among %s with vegan::varpart (%s), reporting adjusted R\u00b2 for each fraction.",
    .mp_and(sprintf("%s (%s)", names(g), vapply(g, paste, character(1),
                                                 collapse = ", "))),
    cite("vegan"))
}

#' @noRd
.mp_m_longitudinal <- function(prov, cite) {
  pa <- prov$params
  what <- switch(
    pa$plot.type,
    spaghetti = sprintf("%s diversity of each %s", pa$measure, pa$subject.var),
    trajectory = sprintf(paste("the trajectory of each %s through a principal",
                               "coordinates analysis of %s"),
                         pa$subject.var, .mp_distance_phrase(pa$distance, cite)),
    alluvial = sprintf("the relative abundance of the %d most abundant %s taxa",
                       pa$top.n, pa$tax.level),
    horizon = sprintf("the relative abundance of the %d most abundant %s taxa",
                      pa$top.n, pa$tax.level),
    "the community")
  sprintf("Change over %s was shown as %s%s.", pa$time.var, what,
          if (!is.null(pa$group.by)) sprintf(", by %s", pa$group.by) else "")
}

#' @noRd
.mp_m_tree <- function(prov, cite) {
  sprintf(
    "The phylogeny of %d taxa was drawn with ggtree (%s) in a %s layout, annotated with %s.",
    prov$data$n_taxa, cite("ggtree"), prov$params$layout,
    .mp_and(prov$params$annotate))
}

#' @noRd
.mp_m_generic <- function(prov, cite) {
  cap <- trimws(gsub("\\s*\n\\s*", "; ", prov$caption %||% ""))
  if (nzchar(cap)) return(paste0(sub("\\.$", "", cap), "."))
  sprintf("A figure was generated with %s().", prov$fun)
}

#' Closing software sentence with versions
#' @noRd
.mp_software_sentence <- function(versions, cite) {
  pk <- setdiff(names(versions), c("micropubr", "ggplot2", "stats"))
  others <- if (length(pk)) {
    sprintf(", with %s", .mp_and(sprintf("%s %s", pk, versions[pk])))
  } else ""
  sprintf("Figures were produced in R %s.%s (%s) with micropubr %s (%s) and ggplot2 %s (%s)%s.",
          R.version$major, R.version$minor, cite("R"),
          versions[["micropubr"]] %||% "", cite("micropubr"),
          versions[["ggplot2"]] %||% "", cite("ggplot2"), others)
}

# --- References --------------------------------------------------------------

#' One reference: list(short, full)
#' @noRd
.mp_ref <- function(key, versions = character(0)) {
  ver <- function(p) {
    v <- versions[p]
    if (is.na(v)) .mp_versions(p)[[1]] else v
  }
  year <- function(p) {
    d <- tryCatch(utils::packageDate(p), error = function(e) NA)
    if (is.na(d)) format(Sys.Date(), "%Y") else format(d, "%Y")
  }
  switch(
    key,
    micropubr = list(short = "Daniilidis 2026", full = sprintf(
      "Daniilidis K (2026). micropubr: a unified grammar for publication-ready microbiome visualization in R. R package version %s. https://github.com/daniilidisK/micropubr",
      ver("micropubr"))),
    R = list(short = sprintf("R Core Team %s", R.version$year), full = sprintf(
      "R Core Team (%s). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria.",
      R.version$year)),
    ggplot2 = list(short = "Wickham 2016", full =
      "Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag, New York."),
    phyloseq = list(short = "McMurdie and Holmes 2013", full =
      "McMurdie PJ, Holmes S (2013). phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE 8: e61217."),
    vegan = list(short = sprintf("Oksanen et al. %s", year("vegan")), full = sprintf(
      "Oksanen J, Simpson GL, Blanchet FG, et al. (%s). vegan: Community Ecology Package. R package version %s.",
      year("vegan"), ver("vegan"))),
    rstatix = list(short = sprintf("Kassambara %s", year("rstatix")), full = sprintf(
      "Kassambara A (%s). rstatix: Pipe-Friendly Framework for Basic Statistical Tests. R package version %s.",
      year("rstatix"), ver("rstatix"))),
    permanova = list(short = "Anderson 2001", full =
      "Anderson MJ (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecology 26: 32-46."),
    permdisp = list(short = "Anderson 2006", full =
      "Anderson MJ (2006). Distance-based tests for homogeneity of multivariate dispersions. Biometrics 62: 245-253."),
    aitchison = list(short = "Aitchison 1986", full =
      "Aitchison J (1986). The Statistical Analysis of Compositional Data. Chapman and Hall, London."),
    bh = list(short = "Benjamini and Hochberg 1995", full =
      "Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B 57: 289-300."),
    holm = list(short = "Holm 1979", full =
      "Holm S (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6: 65-70."),
    wilcoxon = list(short = "Wilcoxon 1945", full =
      "Wilcoxon F (1945). Individual comparisons by ranking methods. Biometrics Bulletin 1: 80-83."),
    kruskal = list(short = "Kruskal and Wallis 1952", full =
      "Kruskal WH, Wallis WA (1952). Use of ranks in one-criterion variance analysis. Journal of the American Statistical Association 47: 583-621."),
    dunn = list(short = "Dunn 1964", full =
      "Dunn OJ (1964). Multiple comparisons using rank sums. Technometrics 6: 241-252."),
    ancombc2 = list(short = "Lin and Peddada 2024", full =
      "Lin H, Peddada SD (2024). Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. Nature Methods 21: 83-91."),
    deseq2 = list(short = "Love et al. 2014", full =
      "Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15: 550."),
    aldex2 = list(short = "Fernandes et al. 2014", full =
      "Fernandes AD, Reid JN, Macklaim JM, McMurrough TA, Edgell DR, Gloor GB (2014). Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. Microbiome 2: 15."),
    lme4 = list(short = "Bates et al. 2015", full =
      "Bates D, Maechler M, Bolker B, Walker S (2015). Fitting linear mixed-effects models using lme4. Journal of Statistical Software 67(1): 1-48."),
    lmertest = list(short = "Kuznetsova et al. 2017", full =
      "Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package: tests in linear mixed effects models. Journal of Statistical Software 82(13): 1-26."),
    indicspecies = list(short = "De Caceres and Legendre 2009", full =
      "De Caceres M, Legendre P (2009). Associations between species and groups of sites: indices and statistical inference. Ecology 90: 3566-3574."),
    mantel = list(short = "Mantel 1967", full =
      "Mantel N (1967). The detection of disease clustering and a generalized regression approach. Cancer Research 27: 209-220."),
    randomforest = list(short = "Liaw and Wiener 2002", full =
      "Liaw A, Wiener M (2002). Classification and regression by randomForest. R News 2(3): 18-22."),
    proc = list(short = "Robin et al. 2011", full =
      "Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez J-C, Mueller M (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics 12: 77."),
    ggtree = list(short = "Yu et al. 2017", full =
      "Yu G, Smith DK, Zhu H, Guan Y, Lam TT-Y (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution 8: 28-36."),
    unifrac = list(short = "Lozupone and Knight 2005", full =
      "Lozupone C, Knight R (2005). UniFrac: a new phylogenetic method for comparing microbial communities. Applied and Environmental Microbiology 71: 8228-8235."),
    list(short = key, full = key))
}

# --- Alt text ------------------------------------------------------------------

#' Chart type, axes, grouping and size, from the plot object alone
#' @noRd
.mp_alt_generic <- function(p) {
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  # Earlier entries win: ellipses (paths) decorate a scatter plot, points
  # decorate box plots, and so on.
  type <- c(GeomBoxplot = "Box plots", GeomViolin = "Violin plots",
            GeomTile = "Heatmap", GeomRaster = "Heatmap",
            GeomCol = "Bar chart", GeomBar = "Bar chart",
            GeomArea = "Stacked area chart", GeomPoint = "Scatter plot",
            GeomLine = "Line chart", GeomPath = "Line chart")
  hit <- names(type)[names(type) %in% geoms]
  kind <- if (length(hit)) type[[hit[1]]] else "Chart"

  x <- .mp_label_text(p, "x", resolve = TRUE)
  y <- .mp_label_text(p, "y", resolve = TRUE)
  axes <- if (nzchar(y) && nzchar(x)) {
    sprintf(" of %s against %s", y, x)
  } else if (nzchar(y)) sprintf(" of %s", y) else ""

  b <- suppressMessages(suppressWarnings(ggplot2::ggplot_build(p)))
  grp <- ""
  for (aes in c("fill", "colour")) {
    title <- .mp_label_text(p, aes, resolve = TRUE)
    # Count categories from the scale, not the layer data: constant colours
    # (reference lines, outlines) would otherwise be counted as groups.
    sc <- tryCatch(b$plot$scales$get_scales(aes), error = function(e) NULL)
    if (!nzchar(title) || title %in% c(x, y) || is.null(sc)) next
    if (sc$is_discrete()) {
      n <- length(sc$get_limits())
      if (n < 2) next
      grp <- sprintf(", coloured by %s (%d categories)", title, n)
    } else {
      grp <- sprintf(", shaded by %s", title)
    }
    break
  }
  panels <- nrow(b$layout$layout)
  pan <- if (panels > 1) sprintf(", in %d panels", panels) else ""
  n <- if (is.data.frame(p$data)) {
    sprintf(", based on %d data rows", nrow(p$data))
  } else ""
  title <- .mp_label_text(p, "title")
  paste0(if (nzchar(title)) paste0(title, ". ") else "",
         kind, axes, grp, pan, n, ".")
}

#' Result-specific sentence for known micropubr figures
#' @noRd
.mp_alt_specific <- function(prov, p) {
  d <- p$data
  switch(
    prov$fun,
    mp_alpha = {
      y <- rlang::as_label(p$mapping$y)
      x <- rlang::as_label(p$mapping$x)
      if (!y %in% names(d) || !x %in% names(d) || identical(x, ".all")) {
        return(NULL)
      }
      med <- sort(tapply(d[[y]], d[[x]], stats::median, na.rm = TRUE))
      if (length(med) < 2) return(NULL)
      sprintf("Median %s is highest in %s (%.2f) and lowest in %s (%.2f).",
              y, names(med)[length(med)], med[[length(med)]], names(med)[1],
              med[[1]])
    },
    mp_diff = {
      st <- prov$stats
      top <- st$top
      sprintf("%d %s enriched and %d depleted%s.%s",
              st$n_up, if (st$n_up == 1) "feature is" else "features are",
              st$n_down,
              if (length(st$contrasts) == 1) sprintf(" (%s)", st$contrasts) else "",
              if (length(top)) sprintf(" Most significant: %s.", .mp_and(top)) else "")
    },
    mp_taxa = {
      if (!all(c(".taxon", "Abundance") %in% names(d))) return(NULL)
      ab <- sort(tapply(d$Abundance, d$.taxon, mean, na.rm = TRUE),
                 decreasing = TRUE)
      ab <- ab[names(ab) != prov$params$other.label]
      if (!length(ab)) return(NULL)
      sprintf("The most abundant taxa are %s.",
              .mp_and(utils::head(names(ab), 3)))
    },
    NULL)
}
