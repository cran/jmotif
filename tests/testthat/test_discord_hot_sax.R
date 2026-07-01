test_that("find discord with HOT SAX", {
  # As of 1.3.0 HOT-SAX computes the NN distance on Z-NORMED windows (the SAX
  # premise of shape similarity); previously the _znorm calls were commented out
  # and it reported raw-amplitude discords (e.g. pos 411 / 1.5046). The z-normed
  # discords below match the brute-force reference and the saxpy / Java
  # implementations exactly.
  discords <- find_discords_hotsax(ecg0606, 100, 4, 4, 0.01, 4)
  expect_equal(discords[discords$position == 430, ]$nn_distance, 5.279080, tolerance = 1e-6)
  expect_equal(discords[discords$position == 318, ]$nn_distance, 4.175756, tolerance = 1e-6)
  expect_equal(discords[discords$position == 2080, ]$nn_distance, 2.392998, tolerance = 1e-6)
  expect_equal(discords[discords$position == 25, ]$nn_distance, 2.375541, tolerance = 1e-6)
})
