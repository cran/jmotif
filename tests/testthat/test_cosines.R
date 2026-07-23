data("CBF")
w <- 60; p <- 6; a <- 6

bag1 <- manyseries_to_wordbag(CBF[["data_train"]][CBF[["labels_train"]] == 1,],
                              w, p, a, "exact", 0.01)
bag2 <- manyseries_to_wordbag(CBF[["data_train"]][CBF[["labels_train"]] == 2,],
                              w, p, a, "exact", 0.01)
bag3 <- manyseries_to_wordbag(CBF[["data_train"]][CBF[["labels_train"]] == 3,],
                              w, p, a, "exact", 0.01)

sample <- (CBF[["data_test"]][CBF[["labels_test"]] == 3,])[1,]

test_that("test cosine sim #1", {

  tfidf <- bags_to_tfidf( list("cylinder" = bag1, "bell" = bag2, "funnel" = bag3) )
  bag <- series_to_wordbag(sample, w, p, a, "exact", 0.01)
  cosines <- cosine_sim(list("bag" = bag, "tfidf" = tfidf))
  expect_equal(3, which(cosines$cosines == max(cosines$cosines)))

})

test_that("cosine_sim returns 0 for zero-norm bag", {

  tfidf <- data.frame(
    words = c("w1", "w2"),
    classA = c(0, 0),
    classB = c(1, 2),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  bag <- data.frame(words = character(0), counts = integer(0), stringsAsFactors = FALSE)
  cosines <- cosine_sim(list(bag = bag, tfidf = tfidf))
  expect_false(any(is.na(cosines$cosines)))
  expect_equal(cosines$cosines, c(0, 0))

})
