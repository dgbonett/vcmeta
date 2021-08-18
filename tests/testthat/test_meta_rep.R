library(vcmeta)

test_that("replicate.mean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
  
  res <- replicate.mean2(
    .05, 21.9, 16.1, 3.82, 3.21, 40, 40, 25.2, 19.1, 3.98, 3.79, 75, 75
  )
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("replicate.mean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
  
  res <- replicate.mean.ps(
    .05, 
    86.22, 
    70.93, 
    14.89, 
    12.32, 
    .765, 
    20, 
    84.81, 
    77.24,
    15.68, 
    16.95, 
    .702, 
    75
  )
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("replicate.stdmean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  res <- replicate.stdmean2(
    .05, 21.9, 16.1, 3.82, 3.21, 40, 40, 25.2, 19.1, 3.98, 3.79, 75, 75
  )

  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("replicate.stdmean.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  res <- replicate.stdmean.ps(
    alpha = .05, 
    m11 = 86.22, 
    m12 = 70.93, 
    sd11 = 14.89, 
    sd12 = 12.32, 
    cor1 = .765, 
    n1 = 20, 
    m21 = 84.81, 
    m22 = 77.24, 
    sd21 = 15.68, 
    sd22 = 16.95, 
    cor2 = .702, 
    n2 = 75
  )
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("replicate.cor returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")
  
  res <- replicate.cor(.05, .598, 80, .324, 200, 0)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("replicate.gen returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")
  
  res <- replicate.gen(.05, .782, .210, .650, .154)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})


test_that("replicate.gen returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")
  
  res <- replicate.gen(.05, .782, .210, .650, .154)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
})