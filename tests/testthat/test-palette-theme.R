test_that("mp_palette returns valid hex colours of the right length", {
  expect_length(mp_palette("taxa10"), 10)
  expect_length(mp_palette("taxa40"), 40)
  expect_length(mp_palette("diverging", n = 7), 7)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}", mp_palette("taxa20"))))
})

test_that("mp_palette interpolates sequential/diverging ramps to n", {
  expect_length(mp_palette("sequential", n = 50), 50)
  expect_length(mp_palette("viridis", n = 3), 3)
})

test_that("qualitative palette beyond its size falls back with a message", {
  expect_message(out <- mp_palette("taxa10", n = 25),
                 class = "micropubr_palette_fallback")
  expect_length(out, 25)
})

test_that("unknown palette name errors informatively", {
  expect_error(mp_palette("nope"), class = "micropubr_bad_palette")
})

test_that("mp_palettes lists palettes grouped by type", {
  pl <- mp_palettes()
  expect_type(pl, "list")
  expect_true(all(c("qualitative", "sequential", "diverging") %in% names(pl)))
})

test_that("themes are ggplot theme objects", {
  for (thm in list(mp_theme_publication(), mp_theme_nature(),
                   mp_theme_cell(), mp_theme_dark())) {
    expect_s3_class(thm, "theme")
    expect_s3_class(thm, "gg")
  }
})

test_that("mp_theme_publication honours grid argument", {
  expect_s3_class(mp_theme_publication(grid = "both"), "theme")
})

test_that("mp_darken recolours without undoing the figure's layout", {
  ps <- make_test_ps()
  light <- mp_taxa(ps, tax.level = "Phylum", group.by = "Group")
  dark <- mp_darken(light)
  th <- ggplot2::theme_get() + dark$theme
  # Rotated sample labels survive; no grid lines appear.
  expect_equal(ggplot2::calc_element("axis.text.x.bottom", th)$angle,
               ggplot2::calc_element("axis.text.x.bottom",
                                     ggplot2::theme_get() + light$theme)$angle)
  expect_s3_class(ggplot2::calc_element("panel.grid.major.x", th),
                  "element_blank")
  expect_equal(ggplot2::calc_element("plot.background", th)$fill, "#1B1B1B")

  # A fixed black label layer becomes light, and the original is untouched.
  df <- make_da_df()
  p <- mp_diff(df, plot.type = "volcano")
  labels <- function(x) Filter(function(l) inherits(l$geom, "GeomTextRepel"),
                               x$layers)[[1]]
  expect_equal(labels(mp_darken(p))$aes_params$colour, "grey90")
  expect_equal(labels(p)$aes_params$colour, "black")
  expect_s3_class(mp_darken(mp_alpha(ps, measure = c("Shannon", "Observed"),
                                     group.by = "Group")), "patchwork")
  expect_error(mp_darken("x"), class = "micropubr_bad_input")
})
