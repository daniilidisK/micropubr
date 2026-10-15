# Provenance, methods text, alt text and guardrails ------------------------

test_that("mp_methods describes alpha and ordination figures", {
  ps <- make_da_ps(levels = c("A", "B", "C"), n_per = 6)
  p1 <- mp_alpha(ps, measure = "Shannon", group.by = "Group",
                 test = "kruskal")
  p2 <- mp_ordination(ps, group.by = "Group", permutations = 99)
  txt <- suppressWarnings(mp_methods(p1, p2))
  expect_s3_class(txt, "mp_methods")
  expect_match(txt, "Kruskal-Wallis")
  expect_match(txt, "Dunn's tests")
  expect_match(txt, "PERMANOVA")
  expect_match(txt, "Anderson MJ \\(2001\\)")
  expect_match(txt, "micropubr")
  expect_output(print(txt), "References")
})

test_that("format = 'list' returns the components", {
  p <- mp_alpha(make_test_ps(), measure = "Shannon", group.by = "Group")
  out <- mp_methods(p, format = "list")
  expect_named(out, c("text", "references", "notes", "software", "versions",
                      "citations"))
  expect_type(out$references, "character")
  # The software sentence is kept apart, so several figures can share one.
  expect_false(grepl("Figures were produced", out$text))
  expect_match(out$software, "^Figures were produced in R")
  expect_true(all(c("R", "micropubr", "ggplot2") %in% names(out$versions)))
  expect_true(all(c("R", "micropubr", "ggplot2") %in% names(out$citations)))
  md <- mp_methods(p, format = "markdown")
  expect_match(md, "**References**", fixed = TRUE)
})

test_that("differential-abundance methods text reflects the design", {
  ps <- make_da_ps()
  res <- mp_run_da(ps, group.by = "Group", method = "lm", adjust = "Age")
  p <- mp_diff(res)
  out <- mp_methods(p, format = "list")
  expect_match(out$text, "linear models on CLR-transformed counts")
  expect_match(out$text, "adjusted for Age")
  expect_match(out$text, "B vs A")
})

test_that("figures without provenance are reported", {
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
    ggplot2::geom_point()
  expect_error(mp_methods(p), class = "micropubr_no_provenance")
})

test_that("provenance survives composition with mp_combine", {
  ps <- make_test_ps()
  p1 <- mp_alpha(ps, measure = "Shannon", group.by = "Group")
  p2 <- mp_taxa(ps, tax.level = "Phylum")
  comb <- mp_combine(p1, p2)
  provs <- micropubr:::.mp_collect_provenance(comb)
  expect_setequal(vapply(provs, `[[`, "", "fun"), c("mp_alpha", "mp_taxa"))
})

test_that("mp_alt_text describes chart type and result", {
  ps <- make_test_ps()
  p <- mp_alpha(ps, measure = "Shannon", group.by = "Group")
  alt <- mp_alt_text(p)
  expect_type(alt, "character")
  expect_match(alt, "^Box plots of Shannon")
  expect_match(alt, "Median Shannon is highest in")
  o <- mp_ordination(ps, group.by = "Group", permutations = 99)
  expect_match(mp_alt_text(o), "^Scatter plot")
  expect_error(mp_alt_text("not a plot"), class = "micropubr_bad_input")
})

test_that("guardrails warn, annotate the caption, or stay silent", {
  on.exit(mp_reset_defaults(), add = TRUE)
  ps <- make_da_ps(levels = c("A", "B", "C"), n_per = 3)
  expect_warning(mp_alpha(ps, measure = "Shannon", group.by = "Group"),
                 class = "micropubr_guardrail")

  mp_set_defaults(guardrails = "off")
  expect_no_warning(mp_alpha(ps, measure = "Shannon", group.by = "Group"))

  mp_set_defaults(guardrails = "caption")
  p <- mp_ordination(ps, group.by = "Group", permutations = 99)
  expect_match(p$labels$caption, "Note: Group\\(s\\)")
  expect_error(mp_set_defaults(guardrails = "loud"),
               class = "micropubr_bad_default")
})

test_that("guardrail notes are kept for mp_methods even when silenced", {
  on.exit(mp_reset_defaults(), add = TRUE)
  mp_set_defaults(guardrails = "off")
  ps <- make_da_ps(levels = c("A", "B", "C"), n_per = 3)
  p <- mp_alpha(ps, measure = "Shannon", group.by = "Group")
  expect_match(mp_methods(p, format = "list")$notes, "fewer than 5",
               all = FALSE)
})

test_that("DA methods text survives results without recorded versions", {
  ps <- make_da_ps()
  res <- mp_run_da(ps, group.by = "Group", method = c("wilcox", "lm"))
  info <- attr(res, "mp_da")
  info$method <- c("wilcox", "lm", "ancombc2")
  info$versions <- info$versions["micropubr"]
  attr(res, "mp_da") <- info
  txt <- mp_methods(mp_diff(res), format = "list")$text
  expect_match(txt, "^Differential abundance among Group groups")
  expect_match(txt, "ANCOM-BC2")
})
