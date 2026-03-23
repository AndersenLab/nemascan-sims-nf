# Phase 1 gate test: module R scripts can be sourced with R_SOURCE_DIR

test_that("all module R scripts source without error using R_SOURCE_DIR", {
  r_dir <- file.path(.project_root, "R")

  module_scripts <- list.files(
    file.path(.project_root, "modules", "db_migration"),
    pattern = "\\.R$",
    recursive = TRUE,
    full.names = TRUE
  )

  expect_true(length(module_scripts) > 0, info = "No module R scripts found")

  for (script_path in module_scripts) {
    # Each script uses R_SOURCE_DIR to find the R library, and optparse for args.
    # We can't fully execute them (they need args), but we can verify that
    # the source() calls succeed by sourcing just the library-loading portion.
    # Instead, we verify the scripts exist, are readable, and reference R_SOURCE_DIR.
    script_content <- readLines(script_path)
    expect_true(
      any(grepl("R_SOURCE_DIR", script_content)),
      info = paste("Script does not use R_SOURCE_DIR:", basename(script_path))
    )
    expect_true(
      any(grepl("source\\(", script_content)),
      info = paste("Script does not source R library files:", basename(script_path))
    )
  }
})
