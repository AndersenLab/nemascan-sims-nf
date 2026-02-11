# Source R library from project R/ directory
r_dir <- Sys.getenv("R_SOURCE_DIR", unset = "")
if (r_dir == "" || !dir.exists(r_dir)) {
  # Detect project root relative to tests/testthat/
  r_dir <- normalizePath(file.path(dirname(getwd()), "R"), mustWork = FALSE)
  if (!dir.exists(r_dir)) r_dir <- normalizePath("R", mustWork = FALSE)
}
for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

# Helper: path to fixtures directory
fixture_path <- function(filename) {
  normalizePath(file.path("tests", "fixtures", filename), mustWork = FALSE)
}

# Helper: create a temporary database directory (cleaned up automatically)
create_temp_db <- function() {
  tmp <- tempfile(pattern = "testdb_")
  dir.create(tmp, recursive = TRUE)
  withr::defer_parent(unlink(tmp, recursive = TRUE))
  tmp
}
