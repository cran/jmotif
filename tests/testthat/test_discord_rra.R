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

