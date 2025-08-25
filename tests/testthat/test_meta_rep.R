library(vcmeta)

test_that("replicate.mean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")
  
  res <- replicate.mean2(
    .05, 21.9, 16.1, 3.82, 3.21, 40, 40, 25.2, 19.1, 3.98, 3.79, 75, 75
  )
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
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
  
  testthat::expect_snapshot(res)
})


test_that("replicate.stdmean2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")
  
  res <- replicate.stdmean2(
    .05, 21.9, 16.1, 3.82, 3.21, 40, 40, 25.2, 19.1, 3.98, 3.79, 75, 75, 9
  )

  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
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
    n2 = 75,
    stdzr = 0
  )
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})


test_that("replicate.cor returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")
  
  res <- replicate.cor(.05, .598, 80, .324, 200, 0)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})


test_that("replicate.cor.gen returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")

  res <- replicate.cor.gen(.05, .454, .170, .318, .098)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})


test_that("replicate.gen returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")
  
  res <- replicate.gen(.05, .782, .210, .650, .154)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})


test_that("replicate.gen returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")
  
  res <- replicate.gen(.05, .782, .210, .650, .154)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})


test_that("replicate.oddsratio returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")
  
  res <- replicate.prop2(.05, 21, 16, 40, 40, 19, 13, 60, 60)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})

test_that("replicate.prop2 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "exp(LL)", "exp(UL)")  
  
  res <- replicate.oddsratio(.05, 1.39, .302, 1.48, .206)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})

test_that("replicate.slope returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "t", "p", "LL", "UL", "df")  
  
  res <- replicate.slope(.05, 23.4, 5.16, 50, 18.5, 4.48, 90, 4)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})



test_that("replicate.spear returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")  

  res <- replicate.spear(.05, .598, 80, .324, 200)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})


test_that("replicate.mean1 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL", "df")  
  
  res <- replicate.mean1(.05, 21.9, 3.82, 40, 25.2, 3.98, 75)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})


test_that("replicate.prop1 returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")  
  
  res <- replicate.prop1(.05, 21, 300, 35, 400)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})


test_that("replicate.ratio.prop2 returns valid matrix", {
  colnames_expected <- c("Estimate", "LL", "UL")  
  
  res <- replicate.ratio.prop2(.05, 21, 16, 40, 40, 19, 13, 60, 60)
  
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})


test_that("replicate.prop.ps returns valid matrix", {
  colnames_expected <- c("Estimate", "SE", "z", "p", "LL", "UL")  
  
  f1 <- c(42, 2, 15, 61)
  f2 <- c(69, 5, 31, 145)
  res <- replicate.prop.ps(.05, f1, f2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)  
  testthat::expect_snapshot(res)
})


test_that("replicate.agree example", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")  
  
  res <- replicate.agree(.05, 85, 100, 160, 200, 2)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})


test_that("replicate.cronbach example", {
  colnames_expected <- c("Estimate", "SE", "LL", "UL")  
  
  res <- replicate.cronbach(.05, .883, 100, .869, 200, 6)
  
  testthat::expect_equal(class(res), c("matrix", "array"))
  testthat::expect_equal(dim(res), c(4, length(colnames_expected)))
  testthat::expect_equal(colnames(res), colnames_expected)
  
  testthat::expect_snapshot(res)
})