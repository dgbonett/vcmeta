library(vcmeta)

# test_that("meta.ave.fisher returns valid CI vector", {
#   colnames_expected <- c("Estimate",        "LL",        "UL")
#   
#   res <- meta.ave.fisher(0.05, 0.376, .054)
#   
#   testthat::expect_equal(class(res), c("matrix", "array"))
#   testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
#   testthat::expect_equal(colnames(res), colnames_expected)
# })


test_that("cor.from.t returns valid vector", {
  colnames_expected <- c("Estimate")
  
  res <- cor.from.t(9.4, 9.8, 1.26, 1.40, 2.27, 30)

  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(length(res), 1)
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
  
})

test_that("meta.chitest returns valid vector", {
  colnames_expected <- c("Q", "df", "p")
  
  est <- c(.297, .324, .281, .149) 
  se <- c(.082, .051, .047, .094)
  res <- meta.chitest(est, se)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("stdmean2.from.t returns valid matrix", {
  colnames_expected <- c("Estimate", "SE")

  res <- stdmean2.from.t(3.27, 25, 25)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(1, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("table.from.odds returns valid matrix", {
  colnames_expected <- c("cell 11",    "cell 12",    "cell 21",    "cell 22")
  
  res <- table.from.odds(.17, .5, 3.18, 100)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})


test_that("table.from.phi returns valid matrix", {
  colnames_expected <- c("cell 11",    "cell 12",    "cell 21",    "cell 22")
  
  res <- table.from.phi(.28, .64, .38, 200)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(2, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  testthat::expect_snapshot(res)
})
