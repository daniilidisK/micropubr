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
