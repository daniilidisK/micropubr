# Built-in colour palettes ------------------------------------------------

# Qualitative core: Okabe-Ito, the de-facto colourblind-safe 8-colour set,
# reordered so the first colours stay maximally distinct for small n.
.mp_okabe_ito <- c(
  "#0072B2", "#E69F00", "#009E73", "#D55E00",
  "#CC79A7", "#56B4E9", "#F0E442", "#999999"
)

# Paul Tol "muted" -- extends the colourblind-friendly range to ~10.
.mp_tol_muted <- c(
  "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
  "#DDCC77", "#CC6677", "#882255", "#AA4499", "#661100"
)

# A curated 40-colour qualitative palette for taxonomic bars. Distinctness is
# prioritised; near-white and near-black are avoided so fills read cleanly.
.mp_taxa40 <- c(
  "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9",
  "#F0E442", "#332288", "#44AA99", "#117733", "#999933", "#882255",
  "#AA4499", "#661100", "#6699CC", "#888844", "#BB5566", "#004488",
  "#DDAA33", "#BB5566", "#000080", "#F6222E", "#FE00FA", "#16FF32",
  "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF",
  "#AA0DFE", "#C4451C", "#85660D", "#1CBE4F", "#FA0087", "#C075A6",
  "#782AB6", "#822E1C", "#1C7F93", "#683B79"
)

# Registry mapping palette names to (type, colours).
.mp_palettes <- list(
  taxa10     = list(type = "qualitative", colors = .mp_taxa40[1:10]),
  taxa20     = list(type = "qualitative", colors = .mp_taxa40[1:20]),
  taxa40     = list(type = "qualitative", colors = .mp_taxa40),
  okabe      = list(type = "qualitative", colors = .mp_okabe_ito),
  groups2    = list(type = "qualitative", colors = .mp_okabe_ito[1:2]),
  groups3    = list(type = "qualitative", colors = .mp_okabe_ito[1:3]),
  groups4    = list(type = "qualitative", colors = .mp_okabe_ito[1:4]),
  groups5    = list(type = "qualitative", colors = .mp_okabe_ito[1:5]),
  groups6    = list(type = "qualitative", colors = .mp_okabe_ito[1:6]),
  groups7    = list(type = "qualitative", colors = .mp_okabe_ito[1:7]),
  groups8    = list(type = "qualitative", colors = .mp_okabe_ito[1:8]),
  sequential = list(type = "sequential",
                    colors = c("#FFFFD9", "#41B6C4", "#225EA8", "#081D58")),
  viridis    = list(type = "sequential",
                    colors = c("#440154", "#3B528B", "#21908C",
                               "#5DC863", "#FDE725")),
  diverging  = list(type = "diverging",
                    colors = c("#2166AC", "#67A9CF", "#F7F7F7",
                               "#EF8A62", "#B2182B"))
)

#' micropubr colour palettes
#'
#' Colourblind-tolerant palettes designed for microbiome figures: qualitative
#' palettes for taxa and groups, plus sequential and diverging ramps for
#' heatmaps and fold-changes.
#'
#' @param name Palette name. One of `"taxa10"`, `"taxa20"`, `"taxa40"`,
#'   `"okabe"`, `"groups2"`..`"groups8"`, `"sequential"`, `"viridis"`,
#'   `"diverging"`.
#' @param n Number of colours to return. For qualitative palettes, if `n`
#'   exceeds the palette size a distinct fallback ramp is generated (with a
#'   message). For sequential/diverging palettes the colours are interpolated
#'   to exactly `n`. Defaults to the full palette.
#' @param reverse Reverse the colour order.
#'
#' @return A character vector of hex colours.
#' @export
#' @examples
#' mp_palette("taxa10")
#' mp_palette("diverging", n = 7)
mp_palette <- function(name = "taxa10", n = NULL, reverse = FALSE) {
  if (!name %in% names(.mp_palettes)) {
    rlang::abort(
      c(sprintf("Unknown palette \"%s\".", name),
        i = sprintf("Available: %s.",
                    paste(names(.mp_palettes), collapse = ", "))),
      class = "micropubr_bad_palette"
    )
  }
  spec <- .mp_palettes[[name]]
  cols <- spec$colors
  if (!is.null(n)) {
    cols <- .mp_expand_palette(cols, n, type = spec$type)
  }
  if (reverse) cols <- rev(cols)
  cols
}

#' List available palettes
#' @return A character vector of palette names, invisibly printed by type.
#' @export
mp_palettes <- function() {
  types <- vapply(.mp_palettes, `[[`, character(1), "type")
  out <- split(names(.mp_palettes), types)
  out
}

#' Expand or interpolate a colour vector to length `n`
#' @noRd
.mp_expand_palette <- function(cols, n, type = "qualitative") {
  if (n <= length(cols) && type == "qualitative") {
    return(cols[seq_len(n)])
  }
  if (type == "qualitative") {
    # Too many categories for the curated set: fall back to a distinct HCL ramp.
    rlang::inform(
      sprintf(paste0("Requested %d qualitative colours but the palette has ",
                     "%d; generating a distinct fallback ramp."),
              n, length(cols)),
      class = "micropubr_palette_fallback"
    )
    return(grDevices::hcl.colors(n, palette = "Dark 3"))
  }
  grDevices::colorRampPalette(cols)(n)
}

#' Resolve a user `palette` argument to a colour vector of length `n`
#'
#' Accepts `NULL` (pick a sensible default for `type`), a registered palette
#' name, or an explicit vector of colours (recycled/interpolated to `n`).
#' @noRd
.mp_resolve_palette <- function(palette, n, type = "qualitative") {
  default_name <- switch(type,
                         qualitative = if (n <= 8) "groups8" else "taxa40",
                         sequential = "viridis",
                         diverging = "diverging",
                         "taxa40")
  # Fall back to the package-wide default palette (set via mp_set_defaults())
  # before the built-in per-type default.
  if (is.null(palette)) palette <- getOption("micropubr.palette", NULL)
  if (is.null(palette)) {
    return(mp_palette(default_name, n = n))
  }
  if (length(palette) == 1 && palette %in% names(.mp_palettes)) {
    return(mp_palette(palette, n = n))
  }
  # Treat as an explicit colour vector.
  .mp_expand_palette(palette, n, type = type)
}
