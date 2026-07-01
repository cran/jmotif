# As of 1.3.0 the brute-force discord search computes the NN distance on
# Z-NORMED windows (it gained an n_threshold parameter and now z-norms each
# window), matching HOT-SAX and the saxpy / Java implementations. Previously it
# operated on RAW subsequences and reported pos 411 / 1.5046 on the full series.
#
# On ecg0606[1:600], win=100, the z-normed top discords are:
#   #0 at 430, NN distance 5.329944
#   #1 at 173, NN distance 5.114571

test_that("find discord brute force", {
  discords <- find_discords_brute_force(ecg0606[1:600], 100, 2, 0.01)
  expect_equal(discords[discords$position == 430, ]$nn_distance, 5.329944, tolerance = 1e-6)
})
