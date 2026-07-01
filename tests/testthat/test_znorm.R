X <- c(-1, -2, -1, 0, 2, 1, 1, 0)
Xscaled <- as.numeric(X) / 100  # ensure it's a numeric vector

test_that("testing ZNorm computation", {

  # Basic ZNorm: mean 0, POPULATION sd 1 (as of 1.3.0 znorm divides by n, not
  # n-1). R's sd() is the sample sd, so check the population sd directly.
  z <- znorm(X)
  expect_equal(sqrt(mean(z^2)), 1, tolerance = 1e-8)
  expect_equal(mean(z), 0, tolerance = 1e-8)

  # ZNorm with small threshold returns original vector
  expect_equal(znorm(Xscaled, 0.1), Xscaled, tolerance = 1e-8)

  # Invalid input triggers an error
  # expect_error(znorm(list("c", "d")))
})

test_that("concrete znorm examples match expected values", {

  # First example
  dat1 <- c(2.02, 2.33, 2.99, 6.85, 9.2, 8.8, 7.5, 6, 5.85, 3.85, 4.85, 3.85, 2.22, 1.45, 1.34)
  # Expected values use POPULATION std (divide by n), the 1.3.0 convention.
  dat1_znorm <- c(-1.014066, -0.8925349, -0.6337913, 0.8794671, 1.800751, 1.643937,
                  1.13429, 0.5462366, 0.4874312, -0.2966404, 0.09539539, -0.2966404,
                  -0.9356589, -1.237526, -1.28065)
  expect_equal(znorm(dat1), dat1_znorm, tolerance = 1e-6)

  # Test cloning behavior: threshold higher than stdev
  expect_equal(znorm(dat1, 3.0), dat1, tolerance = 1e-6)

  # Second example
  dat2 <- c(0.5, 1.29, 2.58, 3.83, 3.25, 4.25, 3.83, 5.63, 6.44, 6.25, 8.75, 8.83, 3.25, 0.75, 0.72)
  dat2_znorm <- c(-1.33469, -1.03429, -0.5437627, -0.06844565, -0.2889927, 0.09126087,
                  -0.06844565, 0.6160108, 0.9240163, 0.8517681, 1.802402, 1.832822,
                  -0.2889927, -1.239627, -1.251034)
  expect_equal(znorm(dat2), dat2_znorm, tolerance = 1e-6)

})

