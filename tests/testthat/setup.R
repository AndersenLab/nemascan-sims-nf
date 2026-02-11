# Detect project root (works whether run from project root or tests/testthat/)
.project_root <- getwd()
if (basename(.project_root) == "testthat") {
  .project_root <- dirname(dirname(.project_root))  # tests/testthat -> project root
} else if (basename(.project_root) == "tests") {
  .project_root <- dirname(.project_root)
}

# Source R library from project R/ directory
r_dir <- Sys.getenv("R_SOURCE_DIR", unset = "")
if (r_dir == "" || !dir.exists(r_dir)) {
  r_dir <- file.path(.project_root, "R")
}
for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

# Helper: path to fixtures directory
fixture_path <- function(filename) {
  file.path(.project_root, "tests", "fixtures", filename)
}

# Helper: create a temporary database directory (cleaned up automatically)
create_temp_db <- function() {
  tmp <- tempfile(pattern = "testdb_")
  dir.create(tmp, recursive = TRUE)
  withr::defer_parent(unlink(tmp, recursive = TRUE))
  tmp
}
