# Build the pkgdown site with internal-only files hidden from the root so
# pkgdown does not publish them. pkgdown renders every root-level `.md` it finds
# (see `pkgdown:::package_mds()`), so agent and assistant instruction files must
# be moved aside for the duration of the build and restored afterwards.
internal_mds <- c("AGENTS.md", "CLAUDE.md")
moved <- character(0)

for (name in internal_mds) {
  src <- file.path(getwd(), name)
  if (file.exists(src)) {
    bak <- file.path(tempdir(), paste0(name, ".bak"))
    file.copy(src, bak, overwrite = TRUE)
    unlink(src)
    moved <- c(moved, name)
    message("Temporarily removed ", name, " from the root for pkgdown build.")
  }
}

# Restore unconditionally. `finally` runs even if the build errors; a top-level
# `on.exit()` would not fire reliably under Rscript, which would leave these
# files moved out of the root.
restore <- function() {
  for (name in moved) {
    bak <- file.path(tempdir(), paste0(name, ".bak"))
    if (file.exists(bak)) {
      file.copy(bak, file.path(getwd(), name), overwrite = TRUE)
      unlink(bak)
      message("Restored ", name, " to the project root.")
    }
  }
}

tryCatch(
  {
    pkgdown::clean_site()
    pkgdown::build_site(preview = FALSE)
  },
  finally = restore()
)
