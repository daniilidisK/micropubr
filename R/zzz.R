# Internal package state and load hooks ----------------------------------

# Lightweight environment for caching defaults (e.g. global theme overrides).
.micropubr_state <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  # Allow users to override the default base font size globally via option.
  if (is.null(getOption("micropubr.font.size"))) {
    options(micropubr.font.size = 11)
  }
  if (is.null(getOption("micropubr.verbose"))) {
    options(micropubr.verbose = FALSE)
  }
  invisible()
}

.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion(pkgname)
  packageStartupMessage(
    sprintf("micropubr %s \u2014 publication-ready microbiome figures.", v),
    "\nSee `vignette(\"getting-started\", package = \"micropubr\")` to begin."
  )
}
