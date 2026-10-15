# Differential-abundance engine (mp_run_da) and parsing (mp_diff) ---------

test_that("mp_run_da (wilcox) feeds mp_diff", {
  ps <- make_test_ps()
  res <- mp_run_da(ps, group.by = "Group", method = "wilcox")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("feature", "feature_id", "contrast", "method", "log2FC",
                    "lfcSE", "pvalue", "padj", "mean", "robust",
                    "structural_zero") %in% names(res)))
  expect_equal(unique(res$contrast), "Treatment vs Control")
  expect_s3_class(mp_diff(res, plot.type = "volcano"), "ggplot")
})

test_that("planted signal has the right sign and is found", {
  ps <- make_da_ps()
  for (m in c("wilcox", "lm")) {
    res <- mp_run_da(ps, group.by = "Group", method = m)
    up <- res[res$feature_id %in% paste0("OTU", 1:4), ]
    dn <- res[res$feature_id %in% paste0("OTU", 5:6), ]
    expect_true(all(up$log2FC > 1), info = m)
    expect_true(all(dn$log2FC < -1), info = m)
    expect_true(all(up$padj < 0.05), info = m)
  }
})

test_that("multi-group designs compare every level with the reference", {
  ps <- make_da_ps(levels = c("A", "B", "C"), n_per = 6)
  res <- mp_run_da(ps, group.by = "Group")
  expect_setequal(unique(res$contrast), c("B vs A", "C vs A"))
  expect_equal(attr(res, "mp_da")$reference, "A")

  res_ref <- mp_run_da(ps, group.by = "Group", reference = "C")
  expect_setequal(unique(res_ref$contrast), c("A vs C", "B vs C"))

  res_pair <- mp_run_da(ps, group.by = "Group", pairwise = TRUE)
  expect_setequal(unique(res_pair$contrast),
                  c("B vs A", "C vs A", "C vs B"))

  two <- mp_run_da(ps, group.by = "Group", groups = c("B", "A"))
  expect_equal(unique(two$contrast), "A vs B")
})

test_that("design errors are classed", {
  ps <- make_da_ps()
  expect_error(mp_run_da(ps, group.by = "Group", method = "nope"),
               class = "micropubr_bad_method")
  expect_error(mp_run_da(ps, group.by = "Group", groups = c("A", "Z")),
               class = "micropubr_bad_groups")
  expect_error(mp_run_da(ps, group.by = "Group", reference = "Z"),
               class = "micropubr_bad_groups")
  expect_error(mp_run_da(ps, group.by = "Group", adjust = "Group"),
               class = "micropubr_bad_var")
  expect_error(mp_run_da(ps, group.by = "Group", method = c("wilcox", "lm"),
                         engine.args = list(n_cl = 2)),
               class = "micropubr_bad_engine_args")
  one <- phyloseq::prune_samples(
    phyloseq::sample_data(ps)$Group == "A", ps)
  expect_error(mp_run_da(one, group.by = "Group"),
               class = "micropubr_need_groups")
})

test_that("a covariate confounded with the group is rejected by lm", {
  ps <- make_da_ps()
  phyloseq::sample_data(ps)$Dup <- phyloseq::sample_data(ps)$Group
  expect_error(mp_run_da(ps, group.by = "Group", method = "lm",
                         adjust = "Dup"),
               class = "micropubr_da_design")
})

test_that("the Wilcoxon test warns that it ignores covariates", {
  ps <- make_da_ps()
  expect_warning(
    res <- mp_run_da(ps, group.by = "Group", method = "wilcox",
                     adjust = "Age"),
    class = "micropubr_guardrail")
  expect_match(attr(res, "mp_da")$notes, "unadjusted", all = FALSE)
})

test_that("several methods stack into one long table", {
  ps <- make_da_ps()
  res <- suppressWarnings(
    mp_run_da(ps, group.by = "Group", method = c("wilcox", "lm"),
              adjust = "Age", tax.level = "Genus"))
  expect_setequal(unique(res$method), c("wilcox", "lm"))
  expect_equal(sum(res$method == "lm"), sum(res$method == "wilcox"))
  expect_true(all(grepl("^Genus", res$feature)))
  expect_s3_class(mp_diff(res, plot.type = "consensus"), "ggplot")
  expect_s3_class(mp_diff(res, plot.type = "concordance"), "ggplot")
  tab <- mp_diff(res, plot.type = "consensus", return.data = TRUE)
  expect_true(all(c("n_significant", "consensus", "agreement") %in%
                    names(tab)))
})

test_that("the linear mixed model runs with a random intercept", {
  skip_if_not_installed("lmerTest")
  ps <- make_da_ps()
  res <- mp_run_da(ps, group.by = "Group", method = "lm", random = "Subject")
  expect_true(all(is.finite(res$log2FC)))
  expect_true(all(res$padj[res$feature_id %in% paste0("OTU", 1:4)] < 0.05))
})

test_that("ANCOM-BC2, DESeq2 and ALDEx2 engines agree on the planted signal", {
  ps <- make_da_ps(n_per = 10)
  for (m in c("ancombc2", "deseq2", "aldex2")) {
    pkg <- c(ancombc2 = "ANCOMBC", deseq2 = "DESeq2", aldex2 = "ALDEx2")[[m]]
    if (!requireNamespace(pkg, quietly = TRUE)) next
    res <- suppressWarnings(mp_run_da(ps, group.by = "Group", method = m))
    up <- res[res$feature_id %in% paste0("OTU", 1:4) & !res$structural_zero, ]
    expect_true(all(up$log2FC > 1), info = m)
    expect_true(mean(up$padj < 0.05) >= 0.75, info = m)
    if (m == "ancombc2") expect_false(all(is.na(res$robust)))
  }
})

test_that("ancombc is an alias of ancombc2", {
  expect_equal(micropubr:::.mp_da_methods(c("ancombc", "wilcox")),
               c("ancombc2", "wilcox"))
})

# --- Parsing real ANCOM-BC output (fixtures saved from ANCOMBC 2.12.1) -----

test_that("ANCOM-BC2 output: group terms only, no intercept, log2 scale", {
  out <- readRDS(test_path("fixtures", "ancombc2_result.rds"))
  da <- micropubr:::.mp_parse_da(out)
  expect_false(any(grepl("Intercept|depth", da$contrast)))
  expect_setequal(unique(da$contrast),
                  c("Mock vs Feces", "Skin vs Feces", "Tongue vs Feces"))
  tested <- da[da$contrast == "Mock vs Feces" & !da$structural_zero, ]
  ref <- out$res$lfc_SampleTypeMock[match(tested$feature, out$res$taxon)]
  expect_equal(tested$log2FC, ref / log(2))
  expect_equal(tested$robust,
               out$res$passed_ss_SampleTypeMock[match(tested$feature,
                                                      out$res$taxon)])
  # Structural zeros are left out of `res` by ANCOM-BC2 and added back.
  expect_true(any(da$structural_zero))
  expect_true(all(is.na(da$log2FC[da$structural_zero])))
  expect_s3_class(mp_diff(out, plot.type = "zeros"), "ggplot")
  expect_s3_class(mp_diff(out, plot.type = "heatmap", fc.threshold = 0),
                  "ggplot")
})

test_that("ANCOM-BC (v1) output with list-of-data-frames `res` is parsed", {
  out <- readRDS(test_path("fixtures", "ancombc_result.rds"))
  da <- micropubr:::.mp_parse_da(out)
  expect_false(any(grepl("Intercept", da$contrast)))
  tested <- da[da$contrast == "Skin vs Feces" & !da$structural_zero, ]
  ref <- out$res$lfc$SampleTypeSkin[match(tested$feature,
                                          out$res$lfc$taxon)]
  expect_equal(tested$log2FC, ref / log(2))
})

test_that("ALDEx2's diff.btw wins over its unitless effect size", {
  ax <- data.frame(diff.btw = c(2, -3), effect = c(0.5, -0.4),
                   we.ep = c(0.01, 0.02), we.eBH = c(0.02, 0.03),
                   rab.all = c(1, 2), row.names = c("a", "b"))
  da <- micropubr:::.mp_parse_da(ax)
  expect_equal(da$log2FC, c(2, -3))
  expect_equal(unique(da$method), "aldex2")
})

test_that("a named list of results is stacked by method", {
  a <- make_da_df(seed = 1)
  b <- make_da_df(seed = 1)
  b$log2FC <- b$log2FC * 0.9
  da <- micropubr:::.mp_parse_da(list(One = a, Two = b))
  expect_setequal(unique(da$method), c("One", "Two"))
  expect_s3_class(mp_diff(list(One = a, Two = b), plot.type = "consensus"),
                  "ggplot")
  expect_error(mp_diff(a, plot.type = "consensus"),
               class = "micropubr_da_one_method")
  expect_error(micropubr:::.mp_parse_da(list(a, b)),
               class = "micropubr_da_parse")
})

test_that("multi-contrast results facet, select and draw heatmaps", {
  df <- make_da_long()
  expect_s3_class(mp_diff(df, plot.type = "volcano"), "ggplot")
  expect_s3_class(mp_diff(df, plot.type = "heatmap"), "ggplot")
  one <- mp_diff(df, contrast = "B vs A", return.data = TRUE)
  expect_equal(unique(one$contrast), "B vs A")
  expect_error(mp_diff(df, contrast = "nope"),
               class = "micropubr_da_no_contrast")
  expect_error(mp_diff(df, plot.type = "zeros"),
               class = "micropubr_da_no_zeros")
})

test_that("robust.only demotes results that fail the sensitivity analysis", {
  df <- make_da_df()
  df$robust <- TRUE
  df$robust[1:3] <- FALSE
  all_sig <- mp_diff(df, robust.only = FALSE, return.data = TRUE)
  robust <- mp_diff(df, return.data = TRUE)
  expect_true(all(robust$direction[1:3] == "NS"))
  expect_gt(sum(all_sig$direction != "NS"), sum(robust$direction != "NS"))
})

test_that("a covariate collinear with the group raises a guardrail", {
  ps <- make_da_ps(n_per = 10)
  g <- phyloseq::sample_data(ps)$Group
  batch <- g
  batch[c(1, 11)] <- rev(batch[c(1, 11)])   # 90% overlap with the group
  phyloseq::sample_data(ps)$Batch <- batch
  expect_warning(
    res <- mp_run_da(ps, group.by = "Group", method = "lm", adjust = "Batch"),
    "strongly associated", class = "micropubr_guardrail")
  phyloseq::sample_data(ps)$Batch <- rep(c("x", "y"), 10)
  expect_no_warning(
    mp_run_da(ps, group.by = "Group", method = "lm", adjust = "Batch"))
})

test_that("taxa unassigned at the tested rank are labelled by a higher rank", {
  ps <- make_test_ps(ntax = 8)
  tt <- as(phyloseq::tax_table(ps), "matrix")
  tt[c("OTU1", "OTU2"), "Genus"] <- NA
  tt[c("OTU1", "OTU2"), "Phylum"] <- "Firmicutes"
  phyloseq::tax_table(ps) <- phyloseq::tax_table(tt)
  res <- mp_run_da(ps, group.by = "Group", tax.level = "Genus",
                   prevalence = 0)
  # tax_glom() merges the two (now identical) lineages into one taxon.
  labs <- unique(res$feature[res$feature_id %in% c("OTU1", "OTU2")])
  expect_identical(labs, "Unassigned Firmicutes")
})
