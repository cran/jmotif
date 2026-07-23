test_that("find discord with RRA", {
  # RRA z-normalizes subsequences before the distance (matching HOT-SAX and the
  # canonical jMotif GrammarViz RRA), so the top discord is the variable-length
  # region beginning at the ecg0606 anomaly near position 430 -- the same window
  # HOT-SAX reports -- rather than the raw-distance region near 363 of earlier
  # versions.
  discords <- find_discords_rra(ecg0606, 100, 4, 4, "none", 0.01, 4)
  expect_true(discords[1,]$start <= 430)
  expect_true(discords[1,]$end > 430)
})

test_that("RRA tier-B exact span w100 seed=0 (rra_ecg_top_region)", {
  discords <- find_discords_rra(ecg0606, 100, 4, 4, "none", 0.01, 1, seed = 0L)
  expect_equal(discords[1, ]$start, 430L)
  expect_equal(discords[1, ]$end, 531L)
  expect_equal(discords[1, ]$nn_distance, 0.054235681140201275, tolerance = 1e-12)
})

test_that("RRA tier-B exact span w120 seed=0 (rra_ecg_top_region_w120)", {
  discords <- find_discords_rra(ecg0606, 120, 4, 4, "none", 0.01, 1, seed = 0L)
  expect_equal(discords[1, ]$start, 430L)
  expect_equal(discords[1, ]$end, 551L)
  expect_equal(discords[1, ]$nn_distance, 0.04856148607272408, tolerance = 1e-12)
})

test_that("RRA w150/p7/a4 does not emit a zero-gap top discord", {
  discords <- find_discords_rra(ecg0606, 150, 7, 4, "none", 0.01, 1, seed = 0L)
  expect_equal(discords[1, ]$start, 1006L)
  expect_equal(discords[1, ]$end, 1175L)
  expect_false(discords[1, ]$rule_id == -1L && discords[1, ]$length == 1L)
  expect_true(discords[1, ]$length >= 7L)
  expect_true(discords[1, ]$nn_distance > 0)
})
