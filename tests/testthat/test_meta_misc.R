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
