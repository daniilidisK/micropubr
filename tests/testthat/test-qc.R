# Sequence quality-control profiles --------------------------------------

# A deterministic per-cycle quality summary for two files, declining with
# cycle (and lower for the second file), in the shape mp_qc_profile() reads
# from FASTQ via ShortRead.
make_qc_df <- function(files = c("R1.fastq.gz", "R2.fastq.gz"),
                       cycles = 1:80, seed = 1) {
  set.seed(seed)
  grid <- expand.grid(file = files, cycle = cycles, score = 2:40,
                      stringsAsFactors = FALSE)
  mu <- 38 - 0.07 * grid$cycle - ifelse(grid$file == files[2], 5, 0)
  grid$count <- round(3000 * stats::dnorm(grid$score, mu, 2.5))
  grid[grid$count > 0, ]
}

test_that("mp_qc_profile plots from a per-cycle data frame", {
  qc <- make_qc_df()
  expect_s3_class(mp_qc_profile(qc), "ggplot")
  expect_s3_class(mp_qc_profile(qc, plot.type = "line"), "ggplot")
  expect_s3_class(mp_qc_profile(qc, aggregate = TRUE), "ggplot")
})

test_that("mp_qc_profile returns a tidy per-cycle data frame", {
  qc <- make_qc_df()
  out <- mp_qc_profile(qc, return.data = TRUE)
  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("file", "cycle", "score", "count"))
  expect_true(all(out$count > 0))
  # The returned frame round-trips back into the function.
  expect_s3_class(mp_qc_profile(out), "ggplot")
})

test_that("aggregate collapses files into a single panel", {
  qc <- make_qc_df()
  agg <- mp_qc_profile(qc, aggregate = TRUE, return.data = TRUE)
  expect_equal(as.character(unique(agg$file)), "All reads")
})

test_that("show controls which summary overlays are drawn", {
  qc <- make_qc_df()
  expect_s3_class(mp_qc_profile(qc, show = "mean"), "ggplot")
  expect_s3_class(mp_qc_profile(qc, show = "quartiles", plot.type = "line"),
                  "ggplot")
})

test_that("weighted quantiles are monotone and bracket the mean", {
  # A single cycle with a known score distribution.
  d <- data.frame(score = 2:40, count = round(1000 * dnorm(2:40, 30, 3)))
  q <- .mp_qc_wquantile(d$score, d$count, c(0.25, 0.5, 0.75))
  expect_true(q[1] <= q[2] && q[2] <= q[3])
  expect_true(q[1] <= 30 && q[3] >= 30)
})

test_that("quality.line = NULL omits the reference line", {
  qc <- make_qc_df()
  p <- mp_qc_profile(qc, quality.line = NULL)
  has_hline <- any(vapply(p$layers, function(l)
    inherits(l$geom, "GeomHline"), logical(1)))
  expect_false(has_hline)
})

test_that("a bad data frame is rejected with an informative error", {
  bad <- data.frame(foo = 1:3, bar = 4:6)
  expect_error(mp_qc_profile(bad), class = "micropubr_bad_qc_cols")
})

test_that("missing FASTAQ files raise a clear error", {
  expect_error(
    mp_qc_profile(c("does_not_exist_1.fastq", "nope_2.fastq.gz")),
    class = "micropubr_no_fastq")
})

test_that("interactive = TRUE returns a plotly widget", {
  skip_if_not_installed("plotly")
  qc <- make_qc_df()
  expect_s3_class(mp_qc_profile(qc, interactive = TRUE), "plotly")
})
