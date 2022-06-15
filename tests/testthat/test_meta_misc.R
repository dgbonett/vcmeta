library(vcmeta)

test_that("ci.fisher returns valid CI vector", {
  res <- ci.fisher(0.05, 0.50, .10)
  
  testthat::expect_equal(class(res), "numeric")
  testthat::expect_equal(length(res), 2)
})


test_that("cor.from.t returns valid vector", {
  res <- cor.from.t(9.4, 9.8, 1.26, 1.40, 0.6708, 10)
  
  testthat::expect_equal(class(res), "numeric")
  testthat::expect_equal(length(res), 1)
})

test_that("meta.chitest returns valid vector", {
  colnames_expected <- c("Q", "df", "p")
  
  est <- c(.297, .324, .281, .149) 
  se <- c(.082, .051, .047, .094)
  res <- meta.chitest(est, se)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})
