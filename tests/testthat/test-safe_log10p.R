test_that("safe_log10p matches known values", {
  # Fixture-independent ground truth
  expect_equal(safe_log10p(1e-5), 5)
  expect_equal(safe_log10p(1e-100), 100)
  expect_equal(safe_log10p(0), 300)
})

test_that("safe_log10p handles normal p-values", {
  expect_equal(safe_log10p(0.05), -log10(0.05), tolerance = 1e-10)
  expect_equal(safe_log10p(1), 0)
})

test_that("safe_log10p caps infinite values at LOG10P_MAX", {
  expect_equal(safe_log10p(0), LOG10P_MAX)
  expect_equal(safe_log10p(1e-300), 300)
})

test_that("safe_log10p handles vectors", {
  result <- safe_log10p(c(0.05, 0, 1e-300, 1))
  expect_length(result, 4)
  expect_true(all(is.finite(result)))
})
