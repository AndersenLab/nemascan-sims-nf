test_that("mappings_schema excludes N and log10p, includes AF1", {
  s <- mappings_schema()
  field_names <- names(s)
  expect_true("AF1" %in% field_names)
  expect_false("N" %in% field_names)
  expect_false("log10p" %in% field_names)
})

