# setup.R - Convenience loader for all R library files
#
# Sources all 6 R library files in dependency order.
# For use in notebooks and interactive R sessions.
#
# Usage:
#   source("R/setup.R")                           # from project root
#   Sys.setenv(NEMASCAN_R_LIB = "/path/to/R")     # or set env var
#   source("setup.R")

# Skip if library is already loaded (e.g., sourced by test harness)
if (!exists("write_marker_set", mode = "function") ||
    !exists("safe_log10p", mode = "function")) {

  # Determine library directory:
  #   1. NEMASCAN_R_LIB env var (explicit override)
  #   2. Auto-detect from this file's location (when source()'d)
  r_lib_dir <- Sys.getenv("NEMASCAN_R_LIB", unset = "")

  if (r_lib_dir == "" || !dir.exists(r_lib_dir)) {
    # Auto-detect: this file lives in R/, so the library dir is its parent
    this_file <- tryCatch(
      normalizePath(sys.frame(1)$ofile),
      error = function(e) NULL
    )
    if (!is.null(this_file)) {
      r_lib_dir <- dirname(this_file)
    }
  }

  if (r_lib_dir == "" || !dir.exists(r_lib_dir)) {
    stop("Cannot find R library directory. Set NEMASCAN_R_LIB or source() from project root.")
  }

  # Source files in dependency order:
  #   utils.R   - no dependencies (logging, parsing)
  #   io.R      - depends on utils.R (log_msg)
  #   database.R - depends on utils.R (log_msg, validation)
  #   queries.R  - depends on database.R (config, connections)
  #   analysis.R - depends on utils.R, database.R
  #   qtl_database.R - depends on database.R
  source(file.path(r_lib_dir, "utils.R"))
  source(file.path(r_lib_dir, "io.R"))
  source(file.path(r_lib_dir, "database.R"))
  source(file.path(r_lib_dir, "queries.R"))
  source(file.path(r_lib_dir, "analysis.R"))
  source(file.path(r_lib_dir, "qtl_database.R"))
}
