# Microbial genomics modules ---------------------------------------------

make_cov_df <- function(contig = "chr1", n = 300, sample = "s1", seed = 1) {
  set.seed(seed)
  data.frame(
    contig = contig,
    position = seq_len(n),
    depth = pmax(0, round(50 + 20 * sin(seq(0, 6 * pi, length.out = n)) +
                            stats::rnorm(n, 0, 5))),
    sample = sample,
    stringsAsFactors = FALSE
  )
}

test_that("mp_coverage plots from a data frame", {
  cov <- make_cov_df()
  expect_s3_class(mp_coverage(cov), "ggplot")
  expect_s3_class(mp_coverage(cov, plot.type = "area"), "ggplot")
  expect_s3_class(mp_coverage(cov, log.y = TRUE), "ggplot")
})

test_that("mp_coverage windows and returns data", {
  cov <- make_cov_df(n = 500)
  binned <- mp_coverage(cov, window = 50, return.data = TRUE)
  expect_s3_class(binned, "data.frame")
  expect_lt(nrow(binned), 500)
  expect_true(all(c("contig", "position", "depth", "sample") %in%
                    names(binned)))
})

test_that("mp_coverage overlays multiple samples and contigs", {
  cov <- rbind(make_cov_df(contig = "chr1", sample = "s1", seed = 1),
               make_cov_df(contig = "chr2", sample = "s2", seed = 2))
  expect_s3_class(mp_coverage(cov), "ggplot")
})

test_that("mp_coverage accepts highlight and annotation tracks", {
  cov <- make_cov_df()
  hl <- data.frame(start = 50, end = 100)
  ann <- data.frame(start = c(20, 150), end = c(80, 220),
                    strand = c("+", "-"), label = c("geneA", "geneB"))
  expect_s3_class(mp_coverage(cov, highlight = hl, annotations = ann),
                  "ggplot")
})

test_that("mp_coverage reads a samtools depth file", {
  cov <- make_cov_df(n = 100)
  tf <- tempfile(fileext = ".txt")
  utils::write.table(cov[, c("contig", "position", "depth")], tf,
                     sep = "\t", row.names = FALSE, col.names = FALSE,
                     quote = FALSE)
  expect_s3_class(mp_coverage(tf, format = "depth"), "ggplot")
})

test_that("mp_coverage rejects malformed data frames", {
  expect_error(mp_coverage(data.frame(a = 1, b = 2)),
               class = "micropubr_bad_coverage_cols")
})

test_that("mp_circos draws when circlize is available", {
  skip_if_not_installed("circlize")
  genome <- data.frame(contig = c("c1", "c2"), length = c(2000, 1500))
  cov <- rbind(
    data.frame(contig = "c1", position = seq(1, 2000, 40),
               value = runif(50, 10, 100)),
    data.frame(contig = "c2", position = seq(1, 1500, 40),
               value = runif(38, 10, 100)))
  pf <- tempfile(fileext = ".pdf")
  grDevices::pdf(pf)
  on.exit({ grDevices::dev.off(); unlink(pf) }, add = TRUE)
  res <- mp_circos(genome, coverage = cov)
  expect_type(res, "list")
  expect_true("genome" %in% names(res))
})
