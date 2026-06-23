skip_if_no_netcal <- function() {
  testthat::skip_if_not_installed("reticulate")

  if (!nzchar(Sys.getenv("RETICULATE_USE_MANAGED_VENV"))) {
    Sys.setenv(RETICULATE_USE_MANAGED_VENV = "no")
  }

  if (!reticulate::py_available(initialize = FALSE)) {
    testthat::skip("Python is not configured for reticulate.")
  }

  if (!reticulate::py_module_available("netcal")) {
    testthat::skip("Python module netcal is not available.")
  }

  invisible(TRUE)
}

import_netcal <- function(module) {
  skip_if_no_netcal()

  tryCatch(
    reticulate::import(module, convert = TRUE),
    error = function(e) {
      testthat::skip(paste("Could not import", module, "from netcal."))
    }
  )
}
