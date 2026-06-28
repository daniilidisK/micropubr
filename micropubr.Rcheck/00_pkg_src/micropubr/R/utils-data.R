# Internal data-wrangling and validation helpers -------------------------
# None of these are exported. They centralise phyloseq access so the public
# mp_*() functions stay short and consistent.

#' Ensure an optional (Suggests) package is available
#'
#' @param pkg Package name.
#' @param what Human-readable description of the feature requiring it.
#' @param bioc Logical; suggest a Bioconductor install hint instead of CRAN.
#' @noRd
.mp_require <- function(pkg, what = "this feature", bioc = FALSE) {
  ok <- vapply(pkg, requireNamespace, logical(1), quietly = TRUE)
  if (all(ok)) {
    return(invisible(TRUE))
  }
  missing <- pkg[!ok]
  hint <- if (bioc) {
    sprintf('BiocManager::install(c(%s))',
            paste0('"', missing, '"', collapse = ", "))
  } else {
    sprintf('install.packages(c(%s))',
            paste0('"', missing, '"', collapse = ", "))
  }
  rlang::abort(
    c(
      sprintf("%s requires the package(s): %s.",
              what, paste(missing, collapse = ", ")),
      i = paste("Install with:", hint)
    ),
    class = "micropubr_missing_dependency"
  )
}

#' Version-safe ggplot check (is_ggplot() replaced is.ggplot() in 3.5.2)
#' @noRd
.mp_is_ggplot <- function(x) {
  ns <- asNamespace("ggplot2")
  fn <- if (exists("is_ggplot", ns, inherits = FALSE)) "is_ggplot" else
    "is.ggplot"
  get(fn, envir = ns)(x)
}

#' Validate that `object` is a usable phyloseq object
#' @noRd
.mp_check_phyloseq <- function(object, require_tree = FALSE,
                               arg = "object") {
  .mp_require("phyloseq", "micropubr", bioc = TRUE)
  if (!methods::is(object, "phyloseq") &&
      !methods::is(object, "otu_table")) {
    rlang::abort(
      sprintf("`%s` must be a phyloseq object, not <%s>.",
              arg, class(object)[1]),
      class = "micropubr_bad_input"
    )
  }
  if (require_tree && is.null(phyloseq::phy_tree(object, errorIfNULL = FALSE))) {
    rlang::abort(
      c("This plot requires a phylogenetic tree.",
        i = "The supplied phyloseq object has no `phy_tree` slot."),
      class = "micropubr_missing_tree"
    )
  }
  invisible(object)
}

#' Extract sample metadata as a plain data.frame with a `.sample` id column
#' @noRd
.mp_sample_data <- function(object) {
  sd <- phyloseq::sample_data(object, errorIfNULL = FALSE)
  if (is.null(sd)) {
    ids <- phyloseq::sample_names(object)
    return(data.frame(.sample = ids, stringsAsFactors = FALSE,
                      row.names = ids))
  }
  df <- data.frame(sd, check.names = FALSE, stringsAsFactors = FALSE)
  df[[".sample"]] <- rownames(df)
  df
}

#' Validate a metadata column name and return it (or NULL)
#' @noRd
.mp_check_var <- function(object, var, arg = "group.by", required = FALSE) {
  if (is.null(var)) {
    if (required) {
      rlang::abort(sprintf("`%s` is required for this plot.", arg),
                   class = "micropubr_missing_var")
    }
    return(NULL)
  }
  vars <- phyloseq::sample_variables(object)
  bad <- setdiff(var, vars)
  if (length(bad)) {
    rlang::abort(
      c(sprintf("`%s` column(s) not found in sample metadata: %s.",
                arg, paste(bad, collapse = ", ")),
        i = sprintf("Available variables: %s.",
                    paste(vars, collapse = ", "))),
      class = "micropubr_bad_var"
    )
  }
  var
}

#' Convert counts to relative abundance (proportions per sample)
#' @noRd
.mp_to_relative <- function(object) {
  phyloseq::transform_sample_counts(object, function(x) {
    s <- sum(x)
    if (s > 0) x / s else x
  })
}

#' Centred log-ratio transform (with a pseudocount) of the OTU matrix
#' @noRd
.mp_clr <- function(object, pseudocount = 1) {
  mat <- .mp_otu_matrix(object, samples_as_rows = TRUE)
  mat <- mat + pseudocount
  logm <- log(mat)
  clr <- sweep(logm, 1, rowMeans(logm), "-")
  clr
}

#' Extract the OTU table as a base matrix in a known orientation
#'
#' @param samples_as_rows If TRUE rows are samples, columns are taxa.
#' @noRd
.mp_otu_matrix <- function(object, samples_as_rows = TRUE) {
  ot <- phyloseq::otu_table(object)
  mat <- as(ot, "matrix")
  if (phyloseq::taxa_are_rows(object)) {
    # rows are taxa -> transpose to samples x taxa
    if (samples_as_rows) mat <- t(mat)
  } else {
    if (!samples_as_rows) mat <- t(mat)
  }
  mat
}

#' Agglomerate to a taxonomic rank, validating the rank name
#' @noRd
.mp_agglomerate <- function(object, tax.level) {
  if (is.null(tax.level)) {
    return(object)
  }
  ranks <- phyloseq::rank_names(object)
  if (!tax.level %in% ranks) {
    rlang::abort(
      c(sprintf("`tax.level = \"%s\"` is not a valid rank.", tax.level),
        i = sprintf("Available ranks: %s.", paste(ranks, collapse = ", "))),
      class = "micropubr_bad_rank"
    )
  }
  phyloseq::tax_glom(object, taxrank = tax.level, NArm = FALSE)
}

#' Identify the `n` most abundant taxa (by total relative abundance)
#' @noRd
.mp_top_taxa <- function(object, n) {
  rel <- .mp_to_relative(object)
  sums <- phyloseq::taxa_sums(rel)
  ord <- names(sort(sums, decreasing = TRUE))
  utils::head(ord, n)
}

#' Long-format melt with taxon labelling and "Other" collapsing
#'
#' Returns a data.frame from [phyloseq::psmelt()] with an added `.taxon`
#' factor column giving the display label at `tax.level`, where everything
#' outside the top `top.n` taxa is collapsed into `other.label`.
#' @noRd
.mp_melt_taxa <- function(object, tax.level = NULL, top.n = 10,
                          relative = TRUE, other.label = "Other") {
  obj <- .mp_agglomerate(object, tax.level)
  if (relative) obj <- .mp_to_relative(obj)

  keep <- .mp_top_taxa(obj, top.n)
  df <- phyloseq::psmelt(obj)

  label_col <- if (!is.null(tax.level)) tax.level else "OTU"
  raw <- as.character(df[[label_col]])
  raw[is.na(raw)] <- "Unassigned"

  # Map taxa not in the top set to the "Other" bucket.
  keep_labels <- unique(as.character(df[[label_col]][df[["OTU"]] %in% keep]))
  raw[!raw %in% keep_labels] <- other.label

  # Order factor by overall abundance, with Other last.
  lvl_order <- names(.mp_aggregate_abundance(df, label_col))
  lvl_order <- lvl_order[lvl_order %in% raw]
  lvl_order <- c(setdiff(lvl_order, other.label),
                 if (other.label %in% raw) other.label)
  df[[".taxon"]] <- factor(raw, levels = unique(lvl_order))
  df
}

#' Order taxon labels by descending total abundance (helper for factor levels)
#' @noRd
.mp_aggregate_abundance <- function(df, label_col) {
  ab <- tapply(df[["Abundance"]], df[[label_col]], sum, na.rm = TRUE)
  sort(ab, decreasing = TRUE)
}
