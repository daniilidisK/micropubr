# Data import helpers ------------------------------------------------------
# Thin, well-documented wrappers that return a ready-to-plot phyloseq object
# from the common upstream pipeline formats.

#' Import QIIME 2 artifacts into a phyloseq object
#'
#' @param features Path to a feature-table `.qza`
#'   (`FeatureTable[Frequency]`).
#' @param taxonomy Optional path to a taxonomy `.qza`.
#' @param tree Optional path to a rooted-tree `.qza`.
#' @param metadata Optional path to a tab-delimited sample-metadata file, or a
#'   data frame.
#' @return A [phyloseq::phyloseq] object.
#' @export
#' @seealso [mp_import_dada2()], [mp_import_biom()], [mp_import_metaphlan()]
mp_import_qiime2 <- function(features, taxonomy = NULL, tree = NULL,
                             metadata = NULL) {
  .mp_require("qiime2R", "Importing QIIME 2 artifacts")
  qiime2R::qza_to_phyloseq(features = features, taxonomy = taxonomy,
                           tree = tree, metadata = metadata)
}

#' Build a phyloseq object from DADA2 output
#'
#' @param seqtab A DADA2 sequence table (samples x ASVs), e.g. from
#'   `dada2::makeSequenceTable()`.
#' @param taxa Optional taxonomy matrix (ASVs x ranks) from
#'   `dada2::assignTaxonomy()`.
#' @param metadata Optional sample metadata data frame (rows = samples).
#' @param tree Optional `phylo` object.
#' @param short.names Replace long sequence identifiers with `ASV1`, `ASV2`,
#'   ... and store the original sequences as `refseq` (requires `Biostrings`).
#' @return A [phyloseq::phyloseq] object.
#' @export
mp_import_dada2 <- function(seqtab, taxa = NULL, metadata = NULL,
                            tree = NULL, short.names = TRUE) {
  .mp_require("phyloseq", "Importing DADA2 output", bioc = TRUE)
  seqtab <- as.matrix(seqtab)
  comps <- list(phyloseq::otu_table(seqtab, taxa_are_rows = FALSE))

  if (!is.null(taxa)) {
    comps <- c(comps, list(phyloseq::tax_table(as.matrix(taxa))))
  }
  if (!is.null(metadata)) {
    comps <- c(comps, list(phyloseq::sample_data(
      as.data.frame(metadata, stringsAsFactors = FALSE))))
  }
  if (!is.null(tree)) comps <- c(comps, list(phyloseq::phy_tree(tree)))

  ps <- do.call(phyloseq::phyloseq, comps)

  if (short.names) {
    seqs <- phyloseq::taxa_names(ps)
    new_names <- paste0("ASV", seq_along(seqs))
    if (requireNamespace("Biostrings", quietly = TRUE)) {
      refs <- Biostrings::DNAStringSet(seqs)
      names(refs) <- new_names
      phyloseq::taxa_names(ps) <- new_names
      ps <- phyloseq::merge_phyloseq(ps, phyloseq::refseq(refs))
    } else {
      phyloseq::taxa_names(ps) <- new_names
    }
  }
  ps
}

#' Import a BIOM file into a phyloseq object
#'
#' @param path Path to a `.biom` file.
#' @param metadata Optional sample metadata (path or data frame) to merge in.
#' @param parse.taxonomy Taxonomy parsing style: `"default"` or
#'   `"greengenes"`.
#' @param ... Passed to [phyloseq::import_biom()].
#' @return A [phyloseq::phyloseq] object.
#' @export
mp_import_biom <- function(path, metadata = NULL,
                           parse.taxonomy = c("default", "greengenes"),
                           ...) {
  .mp_require("phyloseq", "Importing BIOM files", bioc = TRUE)
  parse.taxonomy <- match.arg(parse.taxonomy)
  pf <- if (parse.taxonomy == "greengenes") {
    phyloseq::parse_taxonomy_greengenes
  } else {
    phyloseq::parse_taxonomy_default
  }
  ps <- phyloseq::import_biom(path, parseFunction = pf, ...)
  if (!is.null(metadata)) ps <- .mp_attach_metadata(ps, metadata)
  ps
}

#' Import a merged MetaPhlAn abundance table
#'
#' Parses the pipe-delimited lineage column of a merged MetaPhlAn3/4 table and
#' builds a phyloseq object at the requested rank.
#'
#' @param path Path to a merged MetaPhlAn table (tab-delimited; first column
#'   the clade lineage).
#' @param metadata Optional sample metadata (path or data frame).
#' @param level Taxonomic level to keep: one of `"Kingdom"`, `"Phylum"`,
#'   `"Class"`, `"Order"`, `"Family"`, `"Genus"`, `"Species"` (default).
#' @return A [phyloseq::phyloseq] object.
#' @export
mp_import_metaphlan <- function(path, metadata = NULL, level = "Species") {
  .mp_require("phyloseq", "Importing MetaPhlAn tables", bioc = TRUE)
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
             "Species")
  level <- rlang::arg_match(level, ranks)
  prefixes <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
  depth <- match(level, ranks)

  raw <- utils::read.delim(path, header = TRUE, check.names = FALSE,
                           comment.char = "", stringsAsFactors = FALSE)
  # Locate the lineage column (first character column) and numeric matrix.
  lineage_col <- names(raw)[1]
  lineage <- raw[[lineage_col]]
  num <- raw[, vapply(raw, is.numeric, logical(1)), drop = FALSE]

  # Keep rows at exactly the requested depth (e.g. ending in s__... for species).
  parts <- strsplit(lineage, "\\|")
  keep <- vapply(parts, length, integer(1)) == depth
  parts <- parts[keep]
  num <- num[keep, , drop = FALSE]
  if (nrow(num) == 0) {
    rlang::abort(sprintf("No clades found at level '%s'.", level),
                 class = "micropubr_metaphlan_empty")
  }

  taxmat <- t(vapply(parts, function(p) {
    vals <- sub("^[kpcofgs]__", "", p)
    length(vals) <- depth
    vals
  }, character(depth)))
  colnames(taxmat) <- ranks[seq_len(depth)]
  ids <- make.unique(taxmat[, depth])
  rownames(taxmat) <- ids
  rownames(num) <- ids

  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(as.matrix(num), taxa_are_rows = TRUE),
    phyloseq::tax_table(taxmat))
  if (!is.null(metadata)) ps <- .mp_attach_metadata(ps, metadata)
  ps
}

#' Convert a (Tree)SummarizedExperiment to a phyloseq object
#'
#' Bioconductor's microbiome stack (`mia`, `miaViz`, the OMA book) stores data
#' in a [TreeSummarizedExperiment][TreeSummarizedExperiment::TreeSummarizedExperiment]
#' rather than a `phyloseq` object. Every `mp_*()` function calls this
#' conversion on its input, so either container can be passed directly; it is
#' exported for when you want the `phyloseq` object itself.
#'
#' The count assay becomes the OTU table, `colData()` the sample data, the
#' taxonomic-rank columns of `rowData()` the taxonomy table and, for a
#' `TreeSummarizedExperiment`, `rowTree()` the phylogeny (pruned and
#' relabelled to the rows it links to). Non-atomic `colData()` columns, which
#' `phyloseq` cannot hold, are dropped.
#'
#' @param x A `TreeSummarizedExperiment` or `SummarizedExperiment`; a
#'   `phyloseq` object is returned unchanged.
#' @param assay.type Name of the assay holding the counts. Defaults to
#'   `"counts"` when present, otherwise the first assay.
#' @return A [phyloseq::phyloseq] object.
#' @export
#' @examples
#' if (requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
#'   counts <- matrix(rpois(40, 20), nrow = 5,
#'                    dimnames = list(paste0("OTU", 1:5), paste0("S", 1:8)))
#'   tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
#'     assays = list(counts = counts),
#'     rowData = S4Vectors::DataFrame(Phylum = rep(c("A", "B"), c(2, 3))),
#'     colData = S4Vectors::DataFrame(Group = rep(c("x", "y"), 4)))
#'   mp_as_phyloseq(tse)
#' }
mp_as_phyloseq <- function(x, assay.type = NULL) {
  .mp_require("phyloseq", "micropubr", bioc = TRUE)
  if (methods::is(x, "phyloseq") || methods::is(x, "otu_table")) return(x)
  if (!methods::is(x, "SummarizedExperiment")) {
    rlang::abort(
      sprintf("Cannot convert <%s>; expected a (Tree)SummarizedExperiment.",
              class(x)[1]),
      class = "micropubr_bad_input")
  }
  .mp_require("SummarizedExperiment", "Converting a SummarizedExperiment",
              bioc = TRUE)
  assays <- SummarizedExperiment::assayNames(x)
  if (is.null(assay.type)) {
    assay.type <- if ("counts" %in% assays) "counts" else assays[1]
  }
  if (!length(assays) || is.na(assay.type) || !assay.type %in% assays) {
    rlang::abort(
      c(sprintf("Assay \"%s\" not found.", assay.type %||% ""),
        i = sprintf("Available assays: %s.",
                    paste(assays, collapse = ", "))),
      class = "micropubr_bad_input")
  }

  counts <- as.matrix(SummarizedExperiment::assay(x, assay.type))
  if (is.null(rownames(counts))) {
    rownames(counts) <- paste0("taxon", seq_len(nrow(counts)))
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
  }
  comps <- list(phyloseq::otu_table(counts, taxa_are_rows = TRUE))

  rd <- as.data.frame(SummarizedExperiment::rowData(x))
  ranks <- .mp_rank_columns(rd)
  if (length(ranks)) {
    tt <- as.matrix(data.frame(lapply(rd[, ranks, drop = FALSE],
                                      as.character),
                               check.names = FALSE))
    rownames(tt) <- rownames(counts)
    comps <- c(comps, list(phyloseq::tax_table(tt)))
  }

  cd <- SummarizedExperiment::colData(x)
  atomic <- vapply(cd, is.atomic, logical(1))
  if (any(atomic)) {
    sd <- as.data.frame(cd[, atomic, drop = FALSE])
    rownames(sd) <- colnames(counts)
    comps <- c(comps, list(phyloseq::sample_data(sd)))
  }

  tree <- .mp_quiet_phylo(.mp_tse_tree(x, rownames(counts)))
  if (!is.null(tree)) comps <- c(comps, list(phyloseq::phy_tree(tree)))

  .mp_quiet_phylo(do.call(phyloseq::phyloseq, comps))
}

#' Taxonomic-rank columns of a rowData table
#'
#' Prefers the conventional rank names (as `mia` writes them); otherwise falls
#' back to every character/factor column.
#' @noRd
.mp_rank_columns <- function(rd) {
  if (!ncol(rd)) return(character(0))
  known <- c("domain", "kingdom", "phylum", "class", "order", "family",
             "genus", "species")
  hit <- names(rd)[tolower(names(rd)) %in% known]
  if (length(hit)) return(hit)
  names(rd)[vapply(rd, function(v) is.character(v) || is.factor(v),
                   logical(1))]
}

#' The row tree of a TreeSummarizedExperiment, tips relabelled to row names
#' @noRd
.mp_tse_tree <- function(x, row_ids) {
  if (!methods::is(x, "TreeSummarizedExperiment") ||
      !requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    return(NULL)
  }
  tree <- tryCatch(TreeSummarizedExperiment::rowTree(x),
                   error = function(e) NULL)
  if (is.null(tree)) return(NULL)
  links <- TreeSummarizedExperiment::rowLinks(x)
  tips <- as.character(links$nodeLab)
  # Only a one-to-one row-to-tip mapping can be expressed in phyloseq.
  if (anyNA(tips) || anyDuplicated(tips) ||
      !all(tips %in% tree$tip.label)) {
    return(NULL)
  }
  drop <- setdiff(tree$tip.label, tips)
  if (length(drop)) {
    .mp_require("ape", "Pruning the row tree")
    tree <- ape::drop.tip(tree, drop)
  }
  tree$tip.label <- row_ids[match(tree$tip.label, tips)]
  tree
}

#' Merge a metadata file or data frame into a phyloseq object
#' @noRd
.mp_attach_metadata <- function(ps, metadata) {
  md <- if (is.character(metadata) && length(metadata) == 1) {
    utils::read.delim(metadata, header = TRUE, row.names = 1,
                      check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    as.data.frame(metadata, stringsAsFactors = FALSE)
  }
  phyloseq::sample_data(ps) <- phyloseq::sample_data(md)
  ps
}
