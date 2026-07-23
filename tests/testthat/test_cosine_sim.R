a = c(2, 1, 0, 2, 0, 1, 1, 1)
b <- c(2, 1, 1, 1, 1, 0, 1, 1)

test_that("cosine similarity test", {

  expect_equal(0.1784162, as.numeric(cosine_dist(rbind(a,b))), tolerance = 1e-6)

})

test_that("cosine_dist handles zero-norm rows", {

  m <- rbind(c(0, 0, 0), c(1, 1, 1))
  d <- as.numeric(cosine_dist(m))
  expect_false(any(is.na(d)))
  expect_equal(d, 1)

  m2 <- rbind(c(0, 0, 0), c(0, 0, 0))
  expect_equal(as.numeric(cosine_dist(m2)), 1)

})
