test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

# testing that bandpass_emg returns correct numeric vector length
test_that("bandpass_emg returns numeric vector of same length", {
  signal <- rnorm(1000)
  result <- bandpass_emg(signal, fs = 1000)

  expect_type(result, "double")
  expect_equal(length(result), length(signal))
})

# failure test (incorrect input) for above testing
test_that("bandpass_emg rejects invalid input", {
  expect_error(bandpass_emg(NULL, fs = 1000))
  expect_error(bandpass_emg(c(1, NA, 3), fs = 1000))
})
